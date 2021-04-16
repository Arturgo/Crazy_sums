constexpr int PRIME_MODULO = 211;

#include <fstream>
#include <iostream>
#include "arith_f.h"
#include "print.h"
#include "relations.h"
using namespace std;

static HFormula name_append_component(const HFormula& name, FormulaNode::LeafType component, int power)
{
    return HFormulaProduct(name, HFormulaPower(HFormulaLeaf(component), power));
}

static HFormula name_append_component(const HFormula& name, FormulaNode::LeafType component,
                                      int extra, int power)
{
    return HFormulaProduct(name, HFormulaPower(
               HFormulaLeaf(component, (FormulaNode::LeafExtraArg){.k = extra, .l = 0}), power));
}

static HFormula name_make_lfunc(const HFormula& name, int exponent)
{
    return HFormulaLFunction(name, exponent);
}

static void add_relation(RelationGenerator &manager, Latex& latex,
                         int i_lambda, int i_tau, int i_theta, int i_phi, int i_J_2,
                         int i_sigma_1, int i_sigma_2, int i_sigma_3, int i_mu,
                         int i_zeta1,
                         int min_s, int max_s)
{
   auto t1 = std::chrono::high_resolution_clock::now();

   assert(i_mu <= 2);
   assert(i_lambda <= 1);

   FArith formula =
      pow(liouville(), i_lambda)
    * pow(nb_divisors(), i_tau)
    * pow(theta(), i_theta)
    * pow(phi(), i_phi)
    * pow(jordan_totient(2), i_J_2)
    * pow(sigma_k(1), i_sigma_1)
    * pow(sigma_k(2), i_sigma_2)
    * pow(sigma_k(3), i_sigma_3)
    * pow(mobius(), i_mu)
    * pow(zeta_1(), i_zeta1)
   ;

   auto t2 = std::chrono::high_resolution_clock::now();
   std::chrono::duration<float> elapsed = t2 - t1;

   HFormula name = HFormulaOne(); /* https://youtu.be/i8knduidWCw */
   name = name_append_component(name, FormulaNode::LEAF_LIOUVILLE, i_lambda);
   name = name_append_component(name, FormulaNode::LEAF_NBDIVISORS, i_tau);
   name = name_append_component(name, FormulaNode::LEAF_THETA, i_theta);
   name = name_append_component(name, FormulaNode::LEAF_JORDAN_T, 1, i_phi);
   name = name_append_component(name, FormulaNode::LEAF_JORDAN_T, 2, i_J_2);
   name = name_append_component(name, FormulaNode::LEAF_SIGMA, 1, i_sigma_1);
   name = name_append_component(name, FormulaNode::LEAF_SIGMA, 2, i_sigma_2);
   name = name_append_component(name, FormulaNode::LEAF_SIGMA, 3, i_sigma_3);
   name = name_append_component(name, FormulaNode::LEAF_MU_K, 1, i_mu);
   name = name_append_component(name, FormulaNode::LEAF_ZETAK, i_zeta1, ((i_zeta1 > 0) ? 1 : 0));

   for (int s=min_s; s<=max_s; s++) {
      auto t3 = std::chrono::high_resolution_clock::now();
      FArith fformula = formula * pow(inv_id(), s);
      Fraction<Univariate> frac = fformula.get_fraction();
      auto t4 = std::chrono::high_resolution_clock::now();
      HFormula fname = name_make_lfunc(name, s);

      manager.addFraction(fname, frac);

      if (1) {
         elapsed += t4 - t3;
         cout << KBLD << fname << KRST
              << KGRY "   [" << fformula.A.nbCols() << "]  (" << elapsed.count() << "s)" KRST << endl;
         elapsed -= elapsed;
      }
      if (0) {
         cout << frac << endl;
         fname.get()->print_full(latex.stream, 1);
      }
   }
}

int main([[maybe_unused]] int argc, [[maybe_unused]] char *argv[]) {
   precomputeInverses();
   X.setCoeff(1, 1);
   U.setCoeff(0, 1);
   x = Fraction<Univariate>(X);
   u = Fraction<Univariate>(U);
   z = Fraction<Univariate>(Z);
   Latex latex;

   int maxi_lambda = 1;
   int maxi_tau = 1;
   int maxi_theta = 0;
   int maxi_phi = 1;
   int maxi_J_2 = 0;
   int maxi_sigma_1 = 2;
   int maxi_sigma_2 = 1;
   int maxi_sigma_3 = 1;
   int maxi_mu = 0;
   int maxi_zeta1 = 0;
   int maxi_sum = 8;

   auto t1 = std::chrono::high_resolution_clock::now();

   RelationGenerator manager(&latex);

   add_relation(manager, latex, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2+maxi_sum*3);

   for(int i_lambda = 0;i_lambda <= maxi_lambda;i_lambda++) {
   for(int i_tau = 0;i_tau <= maxi_tau;i_tau++) {
   for(int i_theta = 0;i_theta <= maxi_theta;i_theta++) {
   for(int i_phi = 0;i_phi <= maxi_phi;i_phi++) {
   for(int i_J_2 = 0;i_J_2 <= maxi_J_2;i_J_2++) {
   for(int i_sigma_1 = 0;i_sigma_1 <= maxi_sigma_1;i_sigma_1++) {
   for(int i_sigma_2 = 0;i_sigma_2 <= maxi_sigma_2;i_sigma_2++) {
   for(int i_sigma_3 = 0;i_sigma_3 <= maxi_sigma_3;i_sigma_3++) {
   for(int i_mu = 0;i_mu <= maxi_mu;i_mu++) {
   for(int i_zeta1 = 0;i_zeta1 <= maxi_zeta1;i_zeta1++) {
      int sum = 0*i_lambda + 1*i_tau + 0*i_theta + 1*i_phi + 2*i_J_2
                + 1*i_sigma_1 + 2*i_sigma_2 + 3*i_sigma_3 + 0*i_mu
                + 1*i_zeta1;
      if ((i_phi > 0) && (i_sigma_1 > 0) && (i_mu > 0)) {
         continue; /* Avoid generating C-8: J_2µ == φσµ */
      }
      if(i_lambda+i_tau+i_theta+i_phi+i_J_2+i_sigma_1+i_sigma_2+i_sigma_3+i_mu+i_zeta1 > 0) {
         int min_s = sum + 2;
         int max_s = min_s + max(0, maxi_sum);
         add_relation(manager, latex, i_lambda, i_tau, i_theta, i_phi, i_J_2,
                      i_sigma_1, i_sigma_2, i_sigma_3, i_mu, i_zeta1, min_s, max_s);
      }
   }}}}}}}}}}

   auto t2 = std::chrono::high_resolution_clock::now();
   std::chrono::duration<float> e21 = t2 - t1;
   cout.flush();
   cerr << "Data generated ("<< manager.rational_fractions.size() << " fractions)"
        << KGRY << " (" << e21.count() << "s)" KRST << endl;

   auto t3 = std::chrono::high_resolution_clock::now();
   manager.prepareBasis();
   auto t4 = std::chrono::high_resolution_clock::now();
   std::chrono::duration<float> e43 = t4 - t3;
   cerr << "Basis prepared ("<< manager.polynomial_basis.size() << " polynomials)"
        << KGRY << " (" << e43.count() << "s)" KRST << endl;

   manager.printRelations();

   return 0;
}
