#include "matrix.h"
#include "berlekamp.h"

struct RelationGenerator {
   vector<string> names;
   vector<Fraction<Univariate>> rational_fractions;
   vector<Fraction<Univariate>> rational_fractions_BI;

   vector<Univariate> polynomial_basis;
   vector<IntegralsP> polynomial_basis_BI;
   vector<ModularP> polynomial_basis_M;
   
   void addPolynomial(Univariate poly, int index = 0);
   void addPolynomial2(ModularP poly);
   void addPolynomial2(IntegralsP poly, BigInt modulo, BigInt invertTable[]);
   vector<pair<int, int>> decompose(Univariate poly);
   vector<pair<int, int>> decompose(ModularP poly);
   vector<pair<int, int>> decompose(IntegralsP poly, BigInt p, BigInt invertTable[]);
   void addFraction(string name, Fraction<Univariate> frac);
   void addFraction2(string name, Fraction<Univariate> frac);
   void addFraction2(string name, Fraction<Univariate> frac, BigInt modulo, BigInt invertTable[]);
   void polynomial_basis_add(IntegralsP poly);
   void polynomial_basis_add(ModularP poly);
   void printRelations();
   void printRelations2();
   void printRelations(BigInt p, BigInt invertTable[]);
};

void RelationGenerator::addFraction(string name, Fraction<Univariate> frac) {
   names.push_back(name);
   rational_fractions.push_back(frac);
   
   addPolynomial(frac.getNumerator());
   addPolynomial(frac.getDenominator());
}

void RelationGenerator::addFraction2(string name, Fraction<Univariate> frac, BigInt p, BigInt invertTable[]) {
   names.push_back(name);
   rational_fractions.push_back(frac);
   
   addPolynomial2(toBigInt(frac.getNumerator()), p, invertTable);
   addPolynomial2(toBigInt(frac.getDenominator()), p, invertTable);
}

void RelationGenerator::addFraction2(string name, Fraction<Univariate> frac) {
   names.push_back(name);
   rational_fractions.push_back(frac);
   
   //cout << "num" << endl;
   addPolynomial2(toModular(frac.getNumerator()));
   //cout << "denom" << endl;
   addPolynomial2(toModular(frac.getDenominator()));
   //cout << "next" << endl;
}

void RelationGenerator::addPolynomial(Univariate poly, int index) {
   for(int iElement = index;iElement < (int)polynomial_basis.size();iElement++) {
   	if(poly.size() <= 1) {
   		return;
   	}
   	
      Univariate element = polynomial_basis[iElement];
      
      Univariate pgcd = univariate_gcd(poly, element);
      
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

void RelationGenerator::polynomial_basis_add(IntegralsP poly) {
   for(int iElement = 0;iElement < (int)polynomial_basis_BI.size();iElement++) {
      if (polynomial_basis_BI[iElement] == poly) {
         return;
      }
   }
   polynomial_basis_BI.push_back(poly);
}

void RelationGenerator::polynomial_basis_add(ModularP poly) {
   for(int iElement = 0;iElement < (int)polynomial_basis_M.size();iElement++) {
      if (polynomial_basis_M[iElement] == poly) {
         return;
      }
   }
   polynomial_basis_M.push_back(poly);
}

void RelationGenerator::addPolynomial2(ModularP poly) {
   while (poly.size() > 1) {
      //cout << "begins with " << toString(poly, "X") << endl;
      ModularP div = berlekamp(poly);
      //cout << "berlekamped " << toString(div, "X") << endl;
      div.reduce();
      poly.reduce();
      if (div.size() == poly.size()) {
         //cout << "ended" << endl;
         return polynomial_basis_add(poly);
      } else {
         bool t = true;
         //cout << "begin" << endl;
         while (t) {
            //cout << toString(div, "X") << " | " << toString(poly, "X") << endl;
            ModularP q = quotient(poly, div);
            q.reduce();
            if (q.size() == 0) {
               t = false;
            } else {
               poly = q;
            }
         }
         //cout << "end" << endl;
         addPolynomial2(div);
         //cout << "end2" << endl;
      }
   }
}

void RelationGenerator::addPolynomial2(IntegralsP poly, BigInt p, BigInt invertTable[]) {
   while (poly.size() > 1) {
      //cout << "begins with " << toString(poly, "X") << endl;
      IntegralsP div = berlekamp(poly, p, invertTable);
      //cout << "berlekamped" << endl;
      for (int iCoeff = 0; iCoeff < (int)div.size(); iCoeff++) {
         div.setCoeff(iCoeff, div.getCoeff(iCoeff)%p);
      }
      div.reduce();
      poly.reduce();
      if (div.size() == poly.size()) {
         return polynomial_basis_add(poly);
      } else {
         bool t = true;
         //cout << "begin" << endl;
         while (t) {
            //cout << toString(div, "X") << " | " << toString(poly, "X") << endl;
            IntegralsP quotient = division(poly, div, p, invertTable);
            quotient.reduce();
            if (quotient.size() == 0) {
               t = false;
            } else {
               poly = quotient;
            }
         }
         //cout << "end" << endl;
         addPolynomial2(div, p, invertTable);
         //cout << "end2" << endl;
      }
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

vector<pair<int, int>> RelationGenerator::decompose(ModularP poly) {
   vector<pair<int, int>> decomposition;
   for(int iFactor = 0;iFactor < (int)polynomial_basis_M.size();iFactor++) {
      int nb = 0;
      bool t = true;
      while(poly.size() > 1 && t) {
         ModularP q = quotient(poly, polynomial_basis_M[iFactor]);
         if (q.size() == 0) {
            t = false;
         } else {
            nb++;
            poly = q;
            //cout << nb << " " << toString(q, "X") << endl;
         }
      }
      
      if(nb != 0) {
         decomposition.push_back({iFactor, nb});
      }
   }
   return decomposition;
}


vector<pair<int, int>> RelationGenerator::decompose(IntegralsP poly, BigInt modulo, BigInt invertTable[]) {
   vector<pair<int, int>> decomposition;
   for(int iFactor = 0;iFactor < (int)polynomial_basis_BI.size();iFactor++) {
      int nb = 0;
      bool t = true;
      while(poly.size() > 1 && t) {
         //cout << "enter" << endl;
         IntegralsP quotient = division(poly, polynomial_basis_BI[iFactor], modulo, invertTable);
         //cout << "left" << endl;
         if (quotient.size() == 0) {
            t = false;
         } else {
            nb++;
            poly = quotient;
         }
      }
      
      if(nb != 0) {
         decomposition.push_back({iFactor, nb});
      }
   }
   return decomposition;
}

void RelationGenerator::printRelations() {
   Matrix<Rational> decompositions(0, 0);
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

   cout << "size : "  << decompositions.coeffs.size() << " " << decompositions.coeffs[0].size() << endl;

   //Matrix<Rational> relations = row_echelon_form(kernel_basis(decompositions));
   Matrix<Rational> rows = kernel_basis(decompositions);
   
   cerr << "Relations computed. Size : " << endl;
   cerr << rows.nbRows() << " " << rows.nbCols() << endl;
   cerr << "Simplifying.." << endl;
   
   // On vire les colonnes inutiles
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
   
   cerr << "Matrix simplified. Size : " << endl;
   cerr << cleaned_rows.nbRows() << " " << cleaned_rows.nbCols() << endl;
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

void RelationGenerator::printRelations2() {
   Matrix<Rational> decompositions(0, 0);
   for(auto& fraction: rational_fractions) {
      vector<pair<int, int>> numerator = decompose(toModular(fraction.getNumerator()));
      vector<pair<int, int>> denominator = decompose(toModular(fraction.getDenominator()));
      
      vector<Rational> decomposition(polynomial_basis_M.size(), Rational(0));
      
      //cout << toString(fraction.getNumerator(), "X") << " = ";
      for(pair<int, int> poly : numerator) {
         decomposition[poly.first] = decomposition[poly.first] + Rational(poly.second);
         //cout << toString(polynomial_basis_M[poly.first], "X") << "^" << poly.second << " ";
      }
      //cout << endl;
      
      //cout << toString(fraction.getDenominator(), "X") << " = ";
      for(pair<int, int> poly : denominator) {
         decomposition[poly.first] = decomposition[poly.first] - Rational(poly.second);
         //cout << toString(polynomial_basis_M[poly.first], "X") << "^-" << poly.second << " ";
      }
      //cout << endl;

      
      
      decompositions.coeffs.push_back(decomposition);
   }

   cout << "size : "  << decompositions.coeffs.size() << " " << decompositions.coeffs[0].size() << endl;

   //Matrix<Rational> relations = row_echelon_form(kernel_basis(decompositions));
   Matrix<Rational> rows = kernel_basis(decompositions);
   
   cerr << "Relations computed. Size : " << endl;
   cerr << rows.nbRows() << " " << rows.nbCols() << endl;
   cerr << "Simplifying.." << endl;
   
   // On vire les colonnes inutiles
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
   
   cerr << "Matrix simplified. Size : " << endl;
   cerr << cleaned_rows.nbRows() << " " << cleaned_rows.nbCols() << endl;
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


void RelationGenerator::printRelations(BigInt modulo, BigInt invertTable[]) {
   Matrix<Rational> decompositions(0, 0);
   for(auto& fraction: rational_fractions) {
      //cout << "begin decompose" << endl;
      vector<pair<int, int>> numerator = decompose(toBigInt(fraction.getNumerator()), modulo, invertTable);
      //cout << "mid decompose" << endl;
      vector<pair<int, int>> denominator = decompose(toBigInt(fraction.getDenominator()), modulo, invertTable);
      //cout << "end decompose" << endl;
      
      vector<Rational> decomposition(polynomial_basis_BI.size(), Rational(0));
      
      for(pair<int, int> poly : numerator) {
         decomposition[poly.first] = decomposition[poly.first] + Rational(poly.second);
      }
      
      //cout << "end decomposebis" << endl;
      for(pair<int, int> poly : denominator) {
         decomposition[poly.first] = decomposition[poly.first] - Rational(poly.second);
      }
      
      decompositions.coeffs.push_back(decomposition);
   }
   
   cout << "size : "  << decompositions.coeffs.size() << " " << decompositions.coeffs[0].size() << endl;
   //Matrix<Rational> relations = row_echelon_form(kernel_basis(decompositions));
   Matrix<Rational> rows = kernel_basis(decompositions);
   
   cerr << "Relations computed. Size : " << endl;
   cerr << rows.nbRows() << " " << rows.nbCols() << endl;
   cerr << "Simplifying.." << endl;
   
   // On vire les colonnes inutiles
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
   
   cerr << "Matrix simplified. Size : " << endl;
   cerr << cleaned_rows.nbRows() << " " << cleaned_rows.nbCols() << endl;
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

