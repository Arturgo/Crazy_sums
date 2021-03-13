#pragma once
#include <chrono>
#include "matrix.h"
#include "polynomial.h"
#include "berlekamp.h"
#include <iostream>

struct RelationGenerator {
   vector<string> names;
   vector<Fraction<Univariate>> rational_fractions;

   vector<Univariate> polynomial_basis;

   void addPolynomial(Univariate poly, int index = 0);
   vector<pair<int, int>> decompose(Univariate poly);
   void addFraction(string name, Fraction<Univariate> frac);

   void printRelations();
};

void readRelations(RelationGenerator& manager) {
   int nbFractions;
   cin >> nbFractions;
   
   for(int iFraction = 0;iFraction < nbFractions;iFraction++) {
      string name;
      cin >> name;
      
      int degHaut, degBas;
      cin >> degHaut >> degBas;
      
      Univariate haut(0), bas(0);
      
      for(int iDeg = 0;iDeg <= degHaut;iDeg++) {
         int coeff;
         cin >> coeff;
         haut.setCoeff(iDeg, coeff);
      }
      
      for(int iDeg = 0;iDeg <= degBas;iDeg++) {
         int coeff;
         cin >> coeff;
         bas.setCoeff(iDeg, coeff);
      }
      
      manager.addFraction(name, Fraction<Univariate>(haut, bas));
   }
}

void RelationGenerator::addFraction(string name, Fraction<Univariate> frac) {
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

   //Matrix<Rational> relations = LLL(cleaned_rows, Rational(3) / Rational(4));
   //Matrix<Rational> relations = row_echelon_form(cleaned_rows);
   Matrix<Rational> relations = cleaned_rows;

   for(const vector<Rational>& relation : relations.coeffs) {
      for(size_t iCol = 0;iCol < relation.size();iCol++) {
         if(!(relation[iCol] == Rational(0))) {
            cout << names[iCol_in_rows[iCol]] << "^" << toString(relation[iCol]) << " ";
         }
      }
      cout << "= 1" << endl;
   }
}

