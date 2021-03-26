#pragma once
#include <chrono>
#include <iostream>
#include "berlekamp.h"
#include "matrix.h"
#include "polynomial.h"
#include "xrelation.h"

struct RelationGenerator {
   vector<const FormulaName*> names;
   vector<Fraction<Univariate>> rational_fractions;

   vector<Univariate> polynomial_basis;

   void addPolynomial(Univariate poly, int index = 0);
   vector<pair<int, int>> decompose(Univariate poly);
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

vector<pair<int, int>> RelationGenerator::decompose(Univariate poly) {
   vector<pair<int, int>> decomposition;
   for(int iFactor = 0;iFactor < (int)polynomial_basis.size();iFactor++) {
      int nb = 0;
      while(poly.size() > 1 && poly % polynomial_basis[iFactor] == Univariate(0)) {
         poly = poly / polynomial_basis[iFactor];
         nb++;
      }

      if(nb != 0) {
         decomposition.push_back({iFactor, nb});
      }
   }
   return decomposition;
}

void RelationGenerator::printRelations() {
   auto t1 = std::chrono::high_resolution_clock::now();
   Matrix<Rational> decompositions(0, 0);

   cerr << "Factoring fractions : " << rational_fractions.size() << endl;
   for(auto& fraction: rational_fractions) {
      vector<pair<int, int>> numerator = decompose(fraction.getNumerator());
      vector<pair<int, int>> denominator = decompose(fraction.getDenominator());

      vector<Rational> decomposition(polynomial_basis.size(), Rational(0));

      for(pair<int, int> poly : numerator) {
         decomposition[poly.first] = decomposition[poly.first] + Rational(poly.second);
      }

      for(pair<int, int> poly : denominator) {
         decomposition[poly.first] = decomposition[poly.first] - Rational(poly.second);
      }

      decompositions.coeffs.push_back(decomposition);
      cout << decompositions.coeffs[decompositions.coeffs.size() - 1].size() << endl;
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

   for (auto& row: rows.coeffs) {
      cout << row.coeffs.size() << " -> ";
      for (auto& pos : row.coeffs) {
         cout << pos.first << "-" << pos.second << " ";
      }
      cerr << endl;
   }

   /* On vire les colonnes inutiles */
   auto t4 = std::chrono::high_resolution_clock::now();
   Matrix<Rational> cleaned_rows(rows.nbRows(), 0);

   vector<size_t> iCol_in_rows;

   size_t pos = 0;
   for(size_t iCol = 0;iCol < rows.nbCols();iCol++) {
   	bool isNull = true;
   	for(size_t iRow = 0;iRow < rows.nbRows();iRow++) {
   		if(!(rows.coeffs[iRow].getCoeff(iCol) == Rational(0))) {
   			isNull = false;
   		}
   	}

   	if(!isNull) {
   		iCol_in_rows.push_back(iCol);

   		for(size_t iRow = 0;iRow < rows.nbRows();iRow++) {
   			cleaned_rows.coeffs[iRow].setCoeff(pos, rows.coeffs[iRow].getCoeff(iCol));
			}
         pos++;
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
