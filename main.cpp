#include <fstream>
#include <iostream>
#include "arith_f.h"
#include "print.h"
#include "relations.h"
using namespace std;

static FormulaName* name_append_component(const FormulaName *name, std::string component, int power)
{
    return new FormulaNameProduct(name, new FormulaNamePower(new FormulaNameLeaf(component), power));
}

static FormulaName* name_make_lfunc(const FormulaName *name, int exponent)
{
    return new FormulaNameLFunction(name, new FormulaNameLeaf(exponent));
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

   FormulaName* name = new FormulaName();
   name = name_append_component(name, "\\phi{}", 0);
   name = name_append_component(name, "\\phi{}", i_phi);
   name = name_append_component(name, "\\sigma{}_1", i_sigma_1);
   name = name_append_component(name, "\\sigma{}_2", i_sigma_2);

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
   latex_end();
   return 0;
}
