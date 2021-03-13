#include <iostream>
#include "relations.h"
#include "arith_f.h"
using namespace std;

/*Bivariate x, y, u;

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
*/
int main(int argc, char *argv[]) {
   X.setCoeff(1, 1);
   U.setCoeff(0, 1);
   
  // u.setCoeff(0, U);
   //x.setCoeff(0, X);
   //y.setCoeff(1, U);
   int maxi_phi = 2;
   int maxi_sigma_1 = 2;
   int maxi_sigma_2 = 2;
   int maxi_sum = 6;
   if (argc > 1) {
      if (argc == 5) {
         maxi_phi = atoi(argv[1]);
         maxi_sigma_1 = atoi(argv[2]);
         maxi_sigma_2 = atoi(argv[3]);
         maxi_sum = atoi(argv[4]);
      } else {
         cerr << "wrong number of argument you should give 0 (for default values) or 4 (phi, sigma_1, sigma_2, s) not " << argc - 1 << endl;
         exit(-1);
      }
   }
  	
   RelationGenerator manager;
   for(int i_phi = 0;i_phi <= maxi_phi;i_phi++) {
      for(int i_sigma_1 = 0;i_sigma_1 <= maxi_sigma_1;i_sigma_1++) {
      	for(int i_sigma_2 = 0;i_sigma_2 <= maxi_sigma_2;i_sigma_2++) {
      		int sum = i_phi + i_sigma_1 + 2 * i_sigma_2;
		      for(int s = sum + 2;s <= sum + maxi_sum;s++) {
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
