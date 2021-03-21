#include <fstream>
#include <iostream>
#include "arith_f.h"
#include "print.h"
#include "relations.h"
using namespace std;

static FormulaName* name_append_component(const FormulaName *name, FormulaName::LeafType component, int power)
{
    return new FormulaNameProduct(name, new FormulaNamePower(new FormulaNameLeaf(component), power));
}

static FormulaName* name_append_component(const FormulaName *name, FormulaName::LeafType component,
                                          int extra, int power)
{
    return new FormulaNameProduct(name, new FormulaNamePower(
               new FormulaNameLeaf(component, (FormulaName::LeafExtraArg){.k = extra, .l = 0}), power));
}

static FormulaName* name_make_lfunc(const FormulaName *name, int exponent)
{
    return new FormulaNameLFunction(name, exponent);
}

static void add_relation(RelationGenerator &manager, int i_phi,
                         int i_sigma_1, int i_sigma_2, int i_sigma_3,
                         int i_mu, int i_theta,
                         int s)
{
   auto t1 = std::chrono::high_resolution_clock::now();

   assert(i_mu <= 2);

   Fraction<Univariate> frac =
      (pow(phi(), i_phi)
    * pow(theta(), i_theta)
    * pow(sigma_k(1), i_sigma_1)
    * pow(sigma_k(2), i_sigma_2)
    * pow(sigma_k(3), i_sigma_3)
    * pow(mobius(), i_mu)
    * pow(inv_id(), s)).get_fraction();

   auto t2 = std::chrono::high_resolution_clock::now();

   FormulaName* name = new FormulaName();
   name = name_append_component(name, FormulaName::LEAF_JORDAN_T, 1, i_phi);
   name = name_append_component(name, FormulaName::LEAF_THETA, i_theta);
   name = name_append_component(name, FormulaName::LEAF_SIGMA, 1, i_sigma_1);
   name = name_append_component(name, FormulaName::LEAF_SIGMA, 2, i_sigma_2);
   name = name_append_component(name, FormulaName::LEAF_SIGMA, 3, i_sigma_3);
   name = name_append_component(name, FormulaName::LEAF_MU, i_mu);

   name = name_make_lfunc(name, s);

   auto t3 = std::chrono::high_resolution_clock::now();
   manager.addFraction(name, frac);
   auto t4 = std::chrono::high_resolution_clock::now();

   std::chrono::duration<float> e21 = t2 - t1;
   std::chrono::duration<float> e43 = t4 - t3;

   cout << KBLD << *name << KRST
        << KGRY "   (" << e21.count() << "s + " << e43.count() << "s)" KRST << endl;
#if 0
   cout << toString(frac.getNumerator(), "x") << KGRN "/" KRST
        << toString(frac.getDenominator(), "x") << endl;
#endif
   //name->print_full(latex, 1);
}
int main(int argc, char *argv[]) {
   precomputeInverses(53);
   X.setCoeff(1, 1);
   U.setCoeff(0, 1);
   latex_init();

   int maxi_phi = 2;
   int maxi_theta = 0;
   int maxi_sigma_1 = 2;
   int maxi_sigma_2 = 1;
   int maxi_sigma_3 = 1;
   int maxi_mu = 2;
   int maxi_sum = 8;
   if (argc > 1) {
      if (argc == 5) {
         maxi_phi = atoi(argv[1]);
         maxi_sigma_1 = atoi(argv[2]);
         maxi_sigma_2 = atoi(argv[3]);
         maxi_sigma_3 = atoi(argv[4]);
         maxi_mu = atoi(argv[5]);
         maxi_theta = atoi(argv[6]);
         maxi_sum = atoi(argv[7]);
      } else {
         cerr << "wrong number of argument you should give 0 (for default values) or 4 (phi, sigma_1, sigma_2, s) not " << argc - 1 << endl;
         exit(-1);
      }
   }

   RelationGenerator manager;

    for(int s = 2;s <= 2+maxi_sum*3;s++) {
        add_relation(manager, 0, 0, 0, 0, 0, 0, s);
    }

   for(int i_phi = 0;i_phi <= maxi_phi;i_phi++) {
      for(int i_theta = 0;i_theta <= maxi_theta;i_theta++) {
         for(int i_sigma_1 = 0;i_sigma_1 <= maxi_sigma_1;i_sigma_1++) {
            for(int i_sigma_2 = 0;i_sigma_2 <= maxi_sigma_2;i_sigma_2++) {
               for(int i_sigma_3 = 0;i_sigma_3 <= maxi_sigma_3;i_sigma_3++) {
                  for(int i_mu = 0;i_mu <= maxi_mu;i_mu++) {
                     int sum = i_phi + 0*i_theta + i_sigma_1 + 2*i_sigma_2 + 3*i_sigma_3 + 0*i_mu;
                     for(int s = sum + 2;s <= sum + 2 + max(0, maxi_sum);s++) {
                        if(i_phi+i_theta+i_sigma_1+i_sigma_2+i_sigma_3+i_mu > 0) {
                           add_relation(manager, i_phi, i_sigma_1, i_sigma_2, i_sigma_3, i_mu, i_theta, s);
                        }
                     }
                  }
               }
            }
         }
      }
   }

   cerr << "Data generated" << endl;

   manager.printRelations();
   latex_end();
   return 0;
}
