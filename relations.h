#pragma once
#include <chrono>
#include <iostream>
#include <thread>
#include <mutex>
#include "berlekamp.h"
#include "matrix.h"
#include "polynomial.h"
#include "xrelation.h"

struct RelationGenerator {
   vector<const FormulaName*> names;
   vector<Fraction<Univariate>> rational_fractions;

   vector<Univariate> polynomial_basis;

   void addPolynomial(Univariate poly, int index = 0);
   void addFraction(const FormulaName *name, Fraction<Univariate> frac);

   void printRelation(const vector<Rational>& relation, const vector<size_t>& iCol_in_rows);
   void printRelations();
};

void RelationGenerator::addFraction(const FormulaName *name, Fraction<Univariate> frac) {
   names.push_back(name);
   rational_fractions.push_back(frac);

   addPolynomial(frac.getNumerator());
   addPolynomial(frac.getDenominator());
}

void RelationGenerator::addPolynomial(Univariate poly, int index) {
   for(int iElement = index;iElement < (int)polynomial_basis.size();iElement++) {
   	if(poly.size() <= 1) {
   		return;
   	}

      Univariate element = polynomial_basis[iElement];

      Univariate pgcd = gcd(poly, element);

      if(pgcd.size() <= 1) continue;

      polynomial_basis[iElement] = pgcd;

      while(element % pgcd == Univariate(0)) {
      	Univariate save = element;
      	element = element / pgcd;
      }

      addPolynomial(element, iElement);

      while(poly % pgcd == Univariate(0))
      	poly = poly / pgcd;
      iElement--;
   }

   if(poly.size() > 1) {
   	polynomial_basis.push_back(poly);
   }
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
   return decomposition;
}

void decomposition_worker(
	mutex* mtx, 
	vector<Fraction<Univariate>>* waiting_queue,
	vector<Univariate>* basis,
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

	char* nbThreads_string = getenv("NB_THREADS");
	
	size_t nbThreads = 4;
	if(nbThreads_string != NULL) {
		nbThreads = stoi(string(nbThreads_string));
	}
	
   cerr << "Factoring fractions : " << rational_fractions.size() << endl;

   vector<thread> threads(nbThreads);
   mutex mtx;
   vector<Fraction<Univariate>> waiting_queue = rational_fractions;
   
   for(size_t iThread = 0;iThread < nbThreads;iThread++) {
   	threads[iThread] = thread(
   		decomposition_worker, 
   		&mtx, &waiting_queue, &polynomial_basis, &decompositions
   	);
   }
   
   for(size_t iThread = 0;iThread < nbThreads;iThread++) {
   	threads[iThread].join();
   }
   
   cerr << "Done" << endl;

   auto t2 = std::chrono::high_resolution_clock::now();
   //Matrix<Rational> relations = row_echelon_form(kernel_basis(decompositions));
   Matrix<Rational> rows = kernel_basis(decompositions);

   auto t3 = std::chrono::high_resolution_clock::now();

   std::chrono::duration<float> e21 = t2 - t1;
   std::chrono::duration<float> e32 = t3 - t2;
   cerr << "Relations computed (" << e21.count() << "s+ " << e32.count() << "s). Size : "
        << rows.nbRows() << " * " << rows.nbCols() << endl;
   cerr << "Simplifying.." << endl;

   /* On vire les colonnes inutiles */
   auto t4 = std::chrono::high_resolution_clock::now();
   Matrix<Rational> cleaned_rows(rows.nbRows(), 0);

   vector<size_t> iCol_in_rows;

   for(size_t iCol = 0;iCol < rows.nbCols();iCol++) {
   	bool isNull = true;
   	for(size_t iRow = 0;iRow < rows.nbRows();iRow++) {
   		if(!(rows.coeffs[iRow][iCol] == Rational(0))) {
   			isNull = false;
   		}
   	}

   	if(!isNull) {
   		iCol_in_rows.push_back(iCol);

   		for(size_t iRow = 0;iRow < rows.nbRows();iRow++) {
   			cleaned_rows.coeffs[iRow].push_back(rows.coeffs[iRow][iCol]);
			}
   	}
   }
   auto t5 = std::chrono::high_resolution_clock::now();
   std::chrono::duration<float> e54 = t5 - t4;

   cerr << "Matrix simplified (" << e54.count() << "s). Size : "
        << cleaned_rows.nbRows() << " * " << cleaned_rows.nbCols() << endl;
   cerr << "Simplifying.." << endl;

   //Matrix<Rational> relations_matrix = LLL(cleaned_rows, Rational(3) / Rational(4));
   //Matrix<Rational> relations_matrix = row_echelon_form(cleaned_rows);
   Matrix<Rational> relations_matrix = cleaned_rows;

   vector<Relation> relations;

   for(const vector<Rational>& relation_row : relations_matrix.coeffs) {
      relations.push_back(Relation(relation_row, names, iCol_in_rows));
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
