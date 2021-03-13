#include <iostream>
#include "relations.h"
#include "arith_f.h"
using namespace std;

const int PRIME = 53;

int main() {
   X.setCoeff(1, 1);
   U.setCoeff(0, 1);
   
   precomputeInverses(PRIME);
   
   RelationGenerator manager;
   for(int i_phi = 0;i_phi <= 2;i_phi++) {
      for(int i_sigma_1 = 0;i_sigma_1 <= 2;i_sigma_1++) {
      	for(int i_sigma_2 = 0;i_sigma_2 <= 0;i_sigma_2++) {
      		int sum = i_phi + i_sigma_1 + 2 * i_sigma_2;
		      for(int s = sum + 2;s <= sum + 8;s++) {
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
		         
		         manager.addFraction(name, frac);            
		      }
			}
      }
   }
   
   cerr << "Data generated" << endl;
   
   manager.printRelations();
   return 0;
}
