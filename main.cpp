#include <iostream>
#include "arith_f.h"
#include "print.h"
#include "relations.h"
using namespace std;

static void name_append_component(std::string &name, std::string component, int power)
{
   if (power != 0) {
      if (name != "") {
         name += " ";
      }
      name += component;
      if (power != 1) {
         name += "^" + std::to_string(power);
      }
   }
}

static void name_make_lfunc(std::string &name, int exponent)
{
   if (name == "") {
      name = KCYN "ζ(" + to_string(exponent) + ")" KRST;
   } else {
      name = KRED "L(" + name + ", " + to_string(exponent) + ")" KRST;
   }
}

static void add_relation(RelationGenerator &manager, int i_phi,
                         int i_sigma_1, int i_sigma_2, int s)
{
   auto t1 = std::chrono::high_resolution_clock::now();
   Fraction<Univariate> frac =
      (pow(phi(), i_phi)
    * pow(sigma_k(1), i_sigma_1)
    * pow(sigma_k(2), i_sigma_2)
    * pow(inv_id(), s)).get_fraction();
    
   auto t2 = std::chrono::high_resolution_clock::now();

   string name;
   name_append_component(name, "ϕ", i_phi);
   name_append_component(name, "σ_1", i_sigma_1);
   name_append_component(name, "σ_2", i_sigma_2);

   name_make_lfunc(name, s);
   
   auto t3 = std::chrono::high_resolution_clock::now();
   manager.addFraction(name, frac);
   auto t4 = std::chrono::high_resolution_clock::now();

   std::chrono::duration<float> e21 = t2 - t1;
   std::chrono::duration<float> e43 = t4 - t3;
   cout << KBLD + name + KRST
        << KGRY "   (" << e21.count() << "s + " << e43.count() << "s)" KRST << endl;
   cout << toString(frac.getNumerator(), "x") << KGRN "/" KRST
        << toString(frac.getDenominator(), "x") << endl;
}
int main(int argc, char *argv[]) {
   precomputeInverses(53);
   X.setCoeff(1, 1);
   U.setCoeff(0, 1);

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
            for(int s = sum + 2;s <= sum + 2 + maxi_sum;s++) {
               add_relation(manager, i_phi, i_sigma_1, i_sigma_2, s);
            }
         }
      }
   }

   cerr << "Data generated" << endl;

   manager.printRelations();
   return 0;
}
