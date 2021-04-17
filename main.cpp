constexpr int PRIME_MODULO = 997;

#include <fstream>
#include <iostream>
#include "generation.h"
#include "print.h"
#include "relations.h"
using namespace std;

/*
 * Edit me as you wish to change the generated L-functions
 */
static constexpr GenerationConstraintLine generation_constraints_lines[] = {
    /* leaf_type                  , min_exp, max_exp, [extra min/max k, [min/max l]] */
    {FormulaNode::LEAF_LIOUVILLE  ,       0, 1,       {}},
    {FormulaNode::LEAF_TAUK       ,       0, 2,       {2, 2}},
    {FormulaNode::LEAF_THETA      ,       0, 1,       {}},
    {FormulaNode::LEAF_JORDAN_T   ,       0, 1,       {1, 2}},
    {FormulaNode::LEAF_BETA       ,       0, 1,       {1, 3}},
    {FormulaNode::LEAF_SIGMA      ,       0, 2,       {1, 3}},
    {FormulaNode::LEAF_SIGMA_PRIME,       0, 2,       {1, 3}},
    {FormulaNode::LEAF_KSI        ,       0, 1,       {2, 3}},
    {FormulaNode::LEAF_MU         ,       0, 0,       {1, 2}},
    {FormulaNode::LEAF_NU         ,       0, 0,       {2, 2}},
    {FormulaNode::LEAF_ZETAK      ,       0, 0,       {1, 3}},
};

/*
 * Edit me as you wish to change the generated L-functions
 */
static constexpr GenerationConstraint generation_constraints = {
    .lines = generation_constraints_lines,
    .lines_count = sizeof(generation_constraints_lines)/sizeof(generation_constraints_lines[0]),
    .min_sum = 0,
    .max_sum = 8,
    .max_score = 6,
};

int main([[maybe_unused]] int argc, [[maybe_unused]] char *argv[]) {
    precomputeInverses();
    X.setCoeff(1, 1);
    U.setCoeff(0, 1);
    x = Fraction<Univariate>(X);
    u = Fraction<Univariate>(U);
    z = Fraction<Univariate>(Z);
    Latex latex;

    auto t1 = std::chrono::high_resolution_clock::now();
    RelationGenerator manager(&latex);
    add_relations(manager, latex, generation_constraints);

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
