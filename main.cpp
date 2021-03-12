#include <iostream>
#include "relations.h"
#include "arith_f.h"
using namespace std;


#define PRIME 53

Bivariate x, y, u;

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
  	
   RelationGenerator manager;
   for(int i_phi = 0;i_phi <= 2;i_phi++) {
      for(int i_sigma_1 = 0;i_sigma_1 <= 2;i_sigma_1++) {
      	for(int i_sigma_2 = 0;i_sigma_2 <= 2;i_sigma_2++) {
      		int sum = i_phi + i_sigma_1 + 2 * i_sigma_2;
		      for(int s = sum + 2;s <= sum + 6;s++) {
		         Fraction<Univariate> frac = 
		         	(pow(phi(), i_phi)
		          * pow(sigma_k(1), i_sigma_1)
		          * pow(sigma_k(2), i_sigma_2)
		          * pow(inv_id(), s)).get_fraction();
		         
		         string name;
		         if(i_phi != 0)
		         	name += "Phi^" + to_string(i_phi) + " ";
		         if(i_sigma_1 != 0)
		         	name += "Sigma_1^" + to_string(i_sigma_1) + " ";
		         if(i_sigma_2 != 0)
		         	name += "Sigma_2^" + to_string(i_sigma_2) + " ";
		         
		         if(name == "")
		         	name = "1 ";
		         name = "L(" + name + ", " + to_string(s) + ")";
		         
		         cout << name << endl;
		         cout << toString(frac.getNumerator(), "X") << "/" << toString(frac.getDenominator(), "X") << endl;
		         
		         manager.addFraction2(name, frac);            
		      }
			  }
      }
   }
   
   cerr << "Data generated" << endl;
   
   manager.printRelations2();
   return 0;
}
