#pragma once
#include <algorithm>
#include <atomic>
#include <chrono>
#include <deque>
#include <iostream>
#include <mutex>
#include <random>
#include <thread>
#include "berlekamp.h"
#include "matrix.h"
#include "polynomial.h"
#include "xrelation.h"

class RelationGenerator {
private:
   size_t nbThreads;

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

   RelationGenerator() {
      char* nbThreads_string = getenv("NB_THREADS");

      nbThreads = std::thread::hardware_concurrency();
      nbThreads = 1; /* Temporary, default to 1 until we fix concurrency issues */
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

vector<pair<int, int>> decompose(Univariate poly, const vector<Univariate>& basis) {
   vector<pair<int, int>> decomposition;
   for(size_t iFactor = 0;iFactor < basis.size();iFactor++) {
      int nb = 0;
      while(poly.size() > 1 && poly % basis[iFactor] == Univariate(0)) {
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
	mutex* mtx,
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
				mtx->lock();
				/* If we do not take the lock, the element we want to read
				 * from `basis` might be deleted in the middle of the read */
				Univariate element = *((*basis)[iElement]);
				mtx->unlock();

				Univariate pgcd = gcd(poly, element);

				if(pgcd.size() <= 1) break;

				Univariate simplified = element;
				while(simplified % pgcd == Univariate(0)) {
					simplified = simplified / pgcd;
				}

				mtx->lock();
				if(element == *((*basis)[iElement])) {
					if(pgcd.size() != element.size()) {
						Univariate* ptr = new Univariate();
						*ptr = pgcd;

						Univariate* oldElement = (*basis)[iElement];
						(*basis)[iElement] = ptr;
						delete oldElement;
					}

					if(simplified.size() > 1) {
						waiting_queue_mtx->lock();
						waiting_queue->push_back({iElement, simplified});
						waiting_queue_mtx->unlock();
					}

					mtx->unlock();

					while(poly % pgcd == Univariate(0))
						poly = poly / pgcd;

					continue;
				}
				mtx->unlock();
			}

			iElement++;
		}
	}
}

void RelationGenerator::prepareBasis(void) {
   vector<thread> threads(nbThreads);
   mutex mtx;

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

		vector<pair<int, int>> numerator = decompose(fraction.getNumerator(), *basis);
      vector<pair<int, int>> denominator = decompose(fraction.getDenominator(), *basis);

      vector<Rational> decomposition(basis->size(), Rational(0));

      for(pair<int, int> poly : numerator) {
         decomposition[poly.first] = decomposition[poly.first] + Rational(poly.second);
      }

      for(pair<int, int> poly : denominator) {
         decomposition[poly.first] = decomposition[poly.first] - Rational(poly.second);
      }

		mtx->lock();
		decompositions->coeffs[id] = decomposition;
		mtx->unlock();
	}
}

void RelationGenerator::printRelations() {
   auto t1 = std::chrono::high_resolution_clock::now();
   Matrix<Rational> decompositions(rational_fractions.size(), 0);

   cerr << "Factoring fractions : " << rational_fractions.size() << endl;

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

   cerr << "Done" << endl;

   auto t2 = std::chrono::high_resolution_clock::now();
   decompositions.actualizeNCols();

   //Matrix<Rational> relations = row_echelon_form(kernel_basis(decompositions));
   Matrix<Rational> rows = kernel_basis(decompositions);

   auto t3 = std::chrono::high_resolution_clock::now();

   std::chrono::duration<float> e21 = t2 - t1;
   std::chrono::duration<float> e32 = t3 - t2;
   cerr << "Relations computed (" << e21.count() << "s+ " << e32.count() << "s). "
        << "Size : " << rows.nbRows() << " * " << rows.nbCols() << endl;
   cerr << "Simplifying.." << endl;


   //Matrix<Rational> relations_matrix = LLL(rows, Rational(3) / Rational(4));
   //Matrix<Rational> relations_matrix = row_echelon_form(rows);
   Matrix<Rational> relations_matrix = rows;

   vector<Relation> relations;

   for(const auto& relation_row : relations_matrix.coeffs) {
      relations.push_back(Relation(relation_row.coeffs, names));
   }

   for(auto& relation: relations) {
       relation.classify();
   }
   std::sort(relations.begin(), relations.end());

   for(auto& relation: relations) {
       cout << relation << endl;
       relation.print(latex, 1);
   }
}
