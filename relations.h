#pragma once
#include <algorithm>
#include <atomic>
#include <chrono>
#include <deque>
#include <iostream>
#include <random>
#include <shared_mutex>
#include <thread>
#include "berlekamp.h"
#include "matrix.h"
#include "polynomial.h"
#include "print.h"
#include "xrelation.h"

class RelationGenerator {
private:
   size_t nbThreads;
   Latex* latex;

public:
   vector<HFormula> names;
   vector<Fraction<Univariate>> rational_fractions;

   vector<Univariate> polynomials;
   vector<Univariate> polynomial_basis;

   void addPolynomial(Univariate poly, int index = 0);
   void addFraction(HFormula& name, Fraction<Univariate> frac);

   void printRelation(const vector<Rational>& relation, const vector<size_t>& iCol_in_rows);
   void printRelations();

   void prepareBasis(void);
   void shuffleBasis(void);

   RelationGenerator(Latex* _latex) {
      latex = _latex;

      char* nbThreads_string = getenv("NB_THREADS");

      nbThreads = std::thread::hardware_concurrency();
      if(nbThreads_string != NULL) {
         nbThreads = stoi(string(nbThreads_string));
      }
   }
};

void RelationGenerator::addFraction(HFormula& name, Fraction<Univariate> frac) {
   names.push_back(name);
   rational_fractions.push_back(frac);

   polynomials.push_back(frac.getNumerator());
   polynomials.push_back(frac.getDenominator());
}

vector<pair<size_t, Rational>> decompose(Univariate poly, const vector<Univariate>& basis) {
   vector<pair<size_t, Rational>> decomposition;
   for(size_t iFactor = 0;iFactor < basis.size();iFactor++) {
      int nb = 0;
      while(poly.size() > 1 && isMultipleOf(poly, basis[iFactor])) {
         poly = poly / basis[iFactor];
         nb++;
      }

      if(nb != 0) {
         decomposition.push_back({iFactor, nb});
      }
   }

   assert(poly.size() <= 1);

   return decomposition;
}

void factorisation_worker(
	std::shared_mutex* mtx,
	deque<pair<size_t, Univariate>>* waiting_queue,
	mutex* waiting_queue_mtx,
	deque<atomic<Univariate*>>* basis,
	atomic<size_t>* basis_size
) {
	/*
	 * Note: having `mtx` and `waiting_queue_mtx` does not seem to increase performances,
	 *       but at least the future reader will not think there is a hidden dependency.
	 */
	while(true) {
		waiting_queue_mtx->lock();

		if(waiting_queue->empty()) {
			waiting_queue_mtx->unlock();
			return;
		}

		size_t iElement = waiting_queue->back().first;
		Univariate poly = waiting_queue->back().second;
		waiting_queue->pop_back();

		waiting_queue_mtx->unlock();

		while(poly.size() > 1) {
			if(iElement == *basis_size) {
				mtx->lock();

				if(iElement == *basis_size) {
					Univariate* ptr = new Univariate();
					*ptr = poly;
					basis->emplace_back(ptr);
					(*basis_size)++;

					mtx->unlock();
					break;
				}

				mtx->unlock();
			}

			while(true) {
				mtx->lock_shared();
				/* If we do not take the lock, the element we want to read
				 * from `basis` might be deleted in the middle of the read */
				Univariate element = *((*basis)[iElement]);
				mtx->unlock_shared();

				Univariate pgcd = gcd(poly, element);

				if(pgcd.size() <= 1) break;

				bool simplify_poly = true;
				if(pgcd.size() != element.size()) {
					simplify_poly = false;
					mtx->lock();
					if(element == *((*basis)[iElement])) {
						Univariate* ptr = new Univariate();
						*ptr = pgcd;

						Univariate* oldElement = (*basis)[iElement];
						(*basis)[iElement] = ptr;
						delete oldElement;

						/* Strange: computing simplified is faster here than before the lock */
						Univariate simplified = element;
						while(isMultipleOf(simplified, pgcd)) {
							simplified = simplified / pgcd;
						}
						if(simplified.size() > 1) {
							waiting_queue_mtx->lock();
							waiting_queue->push_back({iElement+1, simplified});
							waiting_queue_mtx->unlock();
						}
						simplify_poly = true;
					}
					mtx->unlock();
				}

				if (simplify_poly) {
					while(isMultipleOf(poly, pgcd)) {
						poly = poly / pgcd;
					}
				}
			}

			iElement++;
		}
	}
}

void RelationGenerator::prepareBasis(void) {
   vector<thread> threads(nbThreads);
   std::shared_mutex mtx;

   deque<pair<size_t, Univariate>> waiting_queue;
   mutex waiting_queue_mtx;
   for(Univariate poly : polynomials) {
   	waiting_queue.push_front({0, poly});
   }

   deque<atomic<Univariate*>> basis;
   atomic<size_t> basis_size = 0;

   for(auto& thread_i: threads) {
      thread_i = thread(
         factorisation_worker,
         &mtx, &waiting_queue, &waiting_queue_mtx, &basis, &basis_size
      );
   }

#if 0
   /* Progress-bar */
   size_t initial_queue_size = waiting_queue.size();
   size_t curr_size = initial_queue_size;
   while(curr_size > 0) {
      waiting_queue_mtx.lock();
      curr_size = waiting_queue.size();
      waiting_queue_mtx.unlock();
      cerr << "\x1b[2K\x1b[100D[" << curr_size << "/" << initial_queue_size << "]";
      cerr.flush();
      std::this_thread::sleep_for(std::chrono::milliseconds(100));
   }
   cerr << endl;
#endif

   for(auto& thread_i: threads) {
      thread_i.join();
   }

   for(Univariate* poly : basis) {
      polynomial_basis.push_back(*poly);
      delete poly;
   }

#if 0
   cout << KCYN "BASIS: size:" << polynomial_basis.size() << KRST << endl;
   for(auto poly : polynomial_basis) {
       cout << toString(poly, "x") << endl;
   }
   cout << KCYN "============" KRST << endl;
#endif
}

/* This function is used to stress-test our setup, e.g. relations size should be independent of this */
void RelationGenerator::shuffleBasis(void) {
   auto rng = std::default_random_engine {};
   std::shuffle(std::begin(polynomial_basis), std::end(polynomial_basis), rng);
}

void decomposition_worker(
	mutex* mtx,
	deque<Fraction<Univariate>>* waiting_queue,
	const vector<Univariate>* basis,
	Matrix<Rational>* decompositions) {

	while(true) {
		mtx->lock();

		if(waiting_queue->empty()) {
			mtx->unlock();
			return;
		}

		Fraction<Univariate> fraction = waiting_queue->back();
		waiting_queue->pop_back();
		size_t id = waiting_queue->size();

		mtx->unlock();

		MatrixRow<Rational> numerator = decompose(fraction.getNumerator(), *basis);
		MatrixRow<Rational> denominator = decompose(fraction.getDenominator(), *basis);

		MatrixRow<Rational> decomposition = numerator - denominator;

		mtx->lock();
		decompositions->coeffs[id] = decomposition;
		mtx->unlock();
	}
}

void RelationGenerator::printRelations() {
   auto t1 = std::chrono::high_resolution_clock::now();
   Matrix<Rational> decompositions(rational_fractions.size(), 0);

   vector<thread> threads(nbThreads);
   mutex mtx;
   deque<Fraction<Univariate>> waiting_queue(rational_fractions.begin(), rational_fractions.end());

   for(auto& thread_i: threads) {
      thread_i = thread(
         decomposition_worker,
         &mtx, &waiting_queue, &polynomial_basis, &decompositions
      );
   }

   for(auto& thread_i: threads) {
      thread_i.join();
   }

   decompositions.actualizeNCols();
   auto t2 = std::chrono::high_resolution_clock::now();

   std::chrono::duration<float> e21 = t2 - t1;
   cerr << "Factored " << rational_fractions.size() << " fractions"
        << KGRY << " (" << e21.count() << "s)" KRST << endl;

   auto t3 = std::chrono::high_resolution_clock::now();
   //Matrix<Rational> relations = row_echelon_form(kernel_basis(decompositions));
   
   decompositions = prepare_matrix(decompositions);
   Matrix<Rational> rows = kernel_basis(decompositions);
   auto t4 = std::chrono::high_resolution_clock::now();

   std::chrono::duration<float> e43 = t4 - t3;
   cerr << "Relations computed. " << "Size: " << rows.nbRows() << " * " << rows.nbCols()
        << KGRY << " (" << e43.count() << "s)" KRST << endl;

   cerr << "Simplifying.." << endl;
   //Matrix<Rational> relations_matrix = LLL(rows, Rational(3) / Rational(4));
   //Matrix<Rational> relations_matrix = row_echelon_form(rows);
   Matrix<Rational> relations_matrix = rows;

   vector<Relation> relations;

   for(const auto& relation_row : relations_matrix.coeffs) {
      relations.push_back(Relation(relation_row.coeffs, names));
   }

   auto t5 = std::chrono::high_resolution_clock::now();
   for(auto& relation: relations) {
       relation.classify();
   }
   std::sort(relations.begin(), relations.end());
   auto t6 = std::chrono::high_resolution_clock::now();

   std::chrono::duration<float> e65 = t6 - t5;
   cerr << "Classified " << relations.size() << " relations"
        << KGRY << " (" << e65.count() << "s)" KRST << endl;
    #if DEBUG_TIME_CLASSIFY
    size_t debug_idx = 0;
    for (const auto& duration: relation_time_classify_debug) {
        cerr << "RTCD " << debug_idx++ << " : " << duration.count() << "s" << endl;
    }
    #endif /* DEBUG_TIME_CLASSIFY */


   for(auto& relation: relations) {
       cout << relation << endl;
       relation.print(latex->stream, 1);
   }
}
