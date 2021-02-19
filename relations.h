#include "matrix.h"
#include "polynomial.h"

struct RelationGenerator {
   vector<string> names;
   vector<Fraction<Univariate>> rational_fractions;

   vector<pair<Univariate, bool>> polynomial_basis;
   vector<int> cleaned_basis;
   
   void addPolynomial(Univariate poly);
   vector<pair<int, int>> decompose(Univariate poly);
   void addFraction(string name, Fraction<Univariate> frac);
   
   void printRelations();
};

void RelationGenerator::addFraction(string name, Fraction<Univariate> frac) {
   names.push_back(name);
   rational_fractions.push_back(frac);
   
   addPolynomial(frac.getNumerator());
   addPolynomial(frac.getDenominator());
}

void RelationGenerator::addPolynomial(Univariate poly) {
   if(poly.size() <= 1) {
      return;
   }
   for(int iElement = 0;iElement < (int)polynomial_basis.size();iElement++) {
      pair<Univariate, bool> element = polynomial_basis[iElement];
      
      if(!element.second) continue;
      Univariate pgcd = univariate_gcd(poly, element.first);
      
      if(pgcd.size() == poly.size() && pgcd.size() == element.first.size()) {
         return;
      }
      
      if(pgcd.size() <= 1) continue;
      
      polynomial_basis[iElement].second = false;
      
      addPolynomial(pgcd);
      addPolynomial(poly / pgcd);
      addPolynomial(element.first / pgcd);
      return;
   }
   
   polynomial_basis.push_back({poly, true});
}

vector<pair<int, int>> RelationGenerator::decompose(Univariate poly) {
   vector<pair<int, int>> decomposition;
   int iPoly = 0;
   for(int id : cleaned_basis) {
      int nb = 0;
      while(poly.size() > 1 && poly % polynomial_basis[id].first == Univariate(0)) {
         poly = poly / polynomial_basis[id].first;
         nb++;
      }
      
      if(nb != 0) {
         decomposition.push_back({iPoly, nb});
      }
      iPoly++;
   }
   return decomposition;
}

void RelationGenerator::printRelations() {
   for(size_t iElement = 0;iElement < polynomial_basis.size();iElement++) {
      if(polynomial_basis[iElement].second) {
         cleaned_basis.push_back(iElement);
      }
   }
   
   Matrix<Rational> decompositions(0, 0);
   for(auto& fraction: rational_fractions) {
      vector<pair<int, int>> numerator = decompose(fraction.getNumerator());
      vector<pair<int, int>> denominator = decompose(fraction.getDenominator());
      
      vector<Rational> decomposition(cleaned_basis.size(), Rational(0));
      
      for(pair<int, int> poly : numerator) {
         decomposition[poly.first] = decomposition[poly.first] + Rational(poly.second);
      }
      
      for(pair<int, int> poly : denominator) {
         decomposition[poly.first] = decomposition[poly.first] - Rational(poly.second);
      }
      
      decompositions.coeffs.push_back(decomposition);
   }
   
   //Matrix<Rational> relations = row_echelon_form(kernel_basis(decompositions));
   Matrix<Rational> rows = kernel_basis(decompositions);
   
   cerr << "Relations computed. Size : " << endl;
   cerr << rows.nbRows() << " " << rows.nbCols() << endl;
   cerr << "Simplifying.." << endl;
   
   Matrix<Rational> relations = LLL(rows, Rational(3) / Rational(4));
   
   for(const vector<Rational>& relation : relations.coeffs) {
      for(size_t iCol = 0;iCol < relation.size();iCol++) {
         if(!(relation[iCol] == Rational(0))) {
            cout << names[iCol] << "^" << toString(relation[iCol]) << " ";
         }
      }
      cout << "= 1" << endl;
   }
}

