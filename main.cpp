#include <iostream>
#include "relations.h"
using namespace std;

#define PRIME 53

Univariate X, U;

Bivariate x, y, u;

Fraction<Univariate> get_fraction(int phi, int sigma, int s) {
   Fraction<Univariate> fraction_p(U, U);
   
   Bivariate poly_k = u;
   
   for(int iPhi = 0;iPhi < phi;iPhi++) {
      fraction_p = fraction_p * Fraction<Univariate>(X - U, X);
      poly_k = poly_k * y;
   }
   
   for(int iSigma = 0;iSigma < sigma;iSigma++) {
      fraction_p = fraction_p * Fraction<Univariate>(U, X - U);
      poly_k = poly_k * (x * y - u);
   }
   
   Fraction<Univariate> fraction_k(0);
   for(int iDeg = 0;iDeg < (int)poly_k.size();iDeg++) {
      if(poly_k.getCoeff(iDeg).size() == 0) {
         continue;
      }
      
      Fraction<Univariate> ajout = Fraction<Univariate>(poly_k.getCoeff(iDeg), U) * Fraction<Univariate>(
         U, 
         (U << (s - iDeg)) - U
      );
      
      fraction_k = fraction_k + ajout;
   }
   
   return Fraction<Univariate>(U, U) + fraction_p * fraction_k;
}

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

int main() {
   X.setCoeff(1, 1);
   U.setCoeff(0, 1);
   
   u.setCoeff(0, U);
   x.setCoeff(0, X);
   y.setCoeff(1, U);

   modulo.value = PRIME;
   inverseTable.push_back(BigInt(0));
   for (int i = 1; i < PRIME; i++) {
      int j = 1;
      while ((i*j)%PRIME != 1) {
         j++;
      }
      inverseTable.push_back(BigInt(j));
   }

   BigInt p = PRIME;
   BigInt inv[PRIME];
   for (int i = 1; i < PRIME; i++) {
      int j = 1;
      while ((i*j)%PRIME != 1) {
         j++;
      }
      inv[i] = j;
   }
   RelationGenerator manager;
   for(int phi = 0;phi <= 2;phi++) {
      for(int sigma = 0;sigma <= 2;sigma++) {
         for(int s = phi + sigma + 2;s <= phi + sigma + 5;s++) {
            Fraction<Univariate> frac = get_fraction(phi, sigma, s);
            string name = "L(Phi^" + to_string(phi) + " * Sigma^" + to_string(sigma)
             + ", " + to_string(s) + ")";
            //decompose(frac.getNumerator(), p, inv);
            //decompose(frac.getDenominator(), p, inv);
            cout << name << endl;
            cout << toString(frac.getNumerator(), "X") << "/" << toString(frac.getDenominator(), "X") << endl;
            
            //manager.addFraction(name, frac);            
            manager.addFraction2(name, frac);            
         }
      }
   }
   
   cerr << "Data generated" << endl;
   
   manager.printRelations2();
   return 0;
}
