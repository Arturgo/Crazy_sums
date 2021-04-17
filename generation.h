#pragma once

#include <unordered_map>
#include "arith_f.h"
#include "relations.h"

typedef struct GenerationFacts {
    int i_phi = 0;
    int i_sigma_1 = 0;
    int i_mu = 0;
} GenerationFacts;

typedef struct GenerationConstraintExtraArg {
    int min_k = 0;
    int max_k = 0;
    int min_l = 0;
    int max_l = 0;
} GenerationConstraintExtraArg;

typedef struct GenerationConstraintLine {
    FormulaNode::LeafType leaf_type;
    int min_exp;
    int max_exp;
    GenerationConstraintExtraArg extra_constraint;
} GenerationConstraintLine;

typedef struct GenerationConstraint {
    const GenerationConstraintLine* lines;
    size_t lines_count;
    int min_sum;
    int max_sum;
    int max_score;
} GenerationConstraint;

static bool bad_formula(GenerationFacts facts)
{
    if ((facts.i_phi > 0) && (facts.i_sigma_1 > 0) && (facts.i_mu > 0)) {
        return true; /* Avoid generating C-8: J_2µ == φσµ */
    }
    return false;
}

static HFormula name_append_component(const HFormula& name, FormulaNode::LeafType component,
                                      int extra_k, int extra_l, int power)
{
    return HFormulaProduct(name, HFormulaPower(
               HFormulaLeaf(component, (FormulaNode::LeafExtraArg){.k = extra_k, .l = extra_l}), power));
}

static void add_relations(RelationGenerator &manager, Latex& latex,
                          const GenerationConstraint& generation_constraints, size_t constraint_idx,
                          const FArith& formula, const HFormula& name, int sum, int score,
                          GenerationFacts facts)
{
    if (score > generation_constraints.max_score) {
        return;
    }

    if (constraint_idx < generation_constraints.lines_count) {
        const GenerationConstraintLine* gc = &generation_constraints.lines[constraint_idx];
        if (gc->min_exp == 0) {
            add_relations(manager, latex, generation_constraints, constraint_idx+1,
                          formula, name, sum, score, facts);
        }

        for (int extra_k=gc->extra_constraint.min_k; extra_k<=gc->extra_constraint.max_k; extra_k++) {
        for (int extra_l=gc->extra_constraint.min_l; extra_l<=gc->extra_constraint.max_l; extra_l++) {
        for (int exp=gc->min_exp; exp<=gc->max_exp; exp++) {
            if (exp == 0) {
                continue;
            }
            //TODO: account for time spent here (w.r.t `elapsed`)
            GenerationFacts ffacts = facts;
            FArith fformula = one();
            int sum_extra;
            switch (gc->leaf_type) {
                case FormulaNode::LEAF_LIOUVILLE:
                    fformula = liouville();
                    sum_extra = 0;
                    assert(exp <= 1);
                    break;
                case FormulaNode::LEAF_NBDIVISORS:
                    fformula = nb_divisors();
                    sum_extra = 1;
                    break;
                case FormulaNode::LEAF_THETA:
                    fformula = theta();
                    sum_extra = 0;
                    break;
                case FormulaNode::LEAF_JORDAN_T:
                    fformula = jordan_totient(extra_k);
                    sum_extra = extra_k;
                    if (extra_k == 1) {
                        ffacts.i_phi += exp;
                    }
                    break;
                case FormulaNode::LEAF_SIGMA:
                    fformula = sigma_k(extra_k);
                    sum_extra = extra_k;
                    if (extra_k == 1) {
                        ffacts.i_sigma_1 += exp;
                    }
                    break;
                case FormulaNode::LEAF_MU_K:
                    fformula = mobius_k(extra_k);
                    sum_extra = 0;
                    ffacts.i_mu += exp;
                    assert(exp <= 2);
                    break;
                case FormulaNode::LEAF_ZETAK:
                    fformula = zeta_1();
                    sum_extra = 1;
                    assert(exp <= 1);
                    break;
                default:
                    assert(false); /* Unknown leaf */
            }
            fformula = formula * pow(fformula, exp);
            HFormula nname = name_append_component(name, gc->leaf_type, extra_k, extra_l, exp);
            int ssum = sum + (sum_extra * exp);
            int sscore = score + extra_k + extra_l + exp;
            add_relations(manager, latex, generation_constraints, constraint_idx+1,
                          fformula, nname, ssum, sscore, ffacts);
        }}}
    } else {
        if ((score == 0) || bad_formula(facts)) {
            return;
        }
        int min_s = 2 + sum + generation_constraints.min_sum;
        int max_s = min_s + generation_constraints.max_sum;
        for (int s=min_s; s<=max_s; s++) {
            auto t3 = std::chrono::high_resolution_clock::now();
            FArith fformula = formula * pow(inv_id(), s);
            Fraction<Univariate> frac = fformula.get_fraction();
            auto t4 = std::chrono::high_resolution_clock::now();
            HFormula fname = HFormulaLFunction(name, s);

            manager.addFraction(fname, frac);

            if (1) {
                std::chrono::duration<float> elapsed = t4 - t3;
                cout << KBLD << fname << KRST
                     << KGRY "   [" << fformula.A.nbCols() << "]  (" << elapsed.count() << "s)" KRST << endl;
            }
            if (0) {
                cout << frac << endl;
                fname.get()->print_full(latex.stream, 1);
            }
        }
    }
}

static void add_relations(RelationGenerator &manager, Latex& latex,
                          const GenerationConstraint& generation_constraints)
{
    FArith formula = one();
    HFormula name = HFormulaOne(); /* https://youtu.be/i8knduidWCw */
    int sum = 0;
    int score = 0;
    GenerationFacts facts;

    /* Add Zetas */
    GenerationConstraint zeta_constraints = {
        .lines = NULL,
        .lines_count = 0,
        .min_sum = 0,
        .max_sum = 10*generation_constraints.max_sum,
        .max_score = 1,
    };
    add_relations(manager, latex, zeta_constraints, 0,
                  formula, name, sum, 1, facts);

    /* Add all other relations */
    add_relations(manager, latex, generation_constraints, 0,
                  formula, name, sum, score, facts);
}
