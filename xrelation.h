#pragma once
#include <iostream>
#include <unordered_map>
#include <set>
#include <sstream>
#include <vector>
#include "fraction.h"
#include "hformula.h"
using namespace std;

#define DEBUG_TIME_CLASSIFY 0
#if DEBUG_TIME_CLASSIFY
    vector<std::chrono::duration<float>> relation_time_classify_debug;
#endif /* DEBUG_TIME_CLASSIFY */

class Relation {
private:
    typedef std::pair<HFormula, Rational> Element;
    std::vector<Element> elements;
    bool known = false;
    std::string known_formula;

    void normalize() {
        /* To be called once at the end of the constructor. Improves printing and formula detection */
        /* Try to normalize the L-side. We want the shortest product on bottom, then the bigger exponent on bottom. */
        size_t l_top = 0, l_bottom = 0;
        const FormulaNode* top_n = NULL;
        const FormulaNode* bottom_n = NULL;
        for (const Element& e: elements) {
            const FormulaNode* inner = e.first.get();
            if (inner->isLFunc() && !inner->isZeta()) {
                if (is_positive(e.second)) {
                    l_top++;
                    top_n = inner;
                } else {
                    l_bottom++;
                    bottom_n = inner;
                }
            }
        }
        if ((l_top < l_bottom) ||
            ((l_top == 1) && (l_bottom == 1) &&
             (!top_n->getLFuncExponentS().is_symbolic()) &&
             (!bottom_n->getLFuncExponentS().is_symbolic()) &&
             (top_n->getLFuncExponent() > bottom_n->getLFuncExponent())
            )
           ) {
            for (Element& e: elements) {
                e.second = -e.second;
            }
        }
    }

    std::ostream& print_element(std::ostream& out, bool latex, const Element& element, bool neg_power) const {
        Rational power = element.second;
        if (neg_power) {
            power = -power;
        }

        out << (latex ? "{" : "");
        element.first.get()->print_inside(out, latex);
        out << (latex ? "}" : "");
        if (power != Rational(1)) {
            out << (latex ? "^{" : "^");
            out << toString(power);
            out << (latex ? "}" : "");
        }
        return out;
    }

    typedef std::function<bool (const Element&)> print_condition;
    size_t print_frac_side(std::ostream& out, bool latex, bool neg_power, print_condition condition) const {
        out << (latex ? "{" : "");
        size_t side_size = 0;
        for(auto element: elements) {
            if (condition(element)) {
                if (side_size>0) {
                    out << (latex ? "\\cdot" : " ");
                }
                print_element(out, latex, element, neg_power);
                side_size++;
            }
        }
        if (side_size == 0) {
            out << "1";
        }
        out << (latex ? "}" : "");
        return side_size;
    }

    void print_frac(std::ostream& out, bool latex, bool neg_power, print_condition condition, bool& is_simple) const {
        std::ostringstream num;
        std::ostringstream denom;
        size_t num_size = print_frac_side(num, latex, neg_power,
            [condition, neg_power](const Element& e){ return condition(e) && (neg_power ^  is_positive(e.second));}
        );
        size_t denom_size = print_frac_side(denom, latex, neg_power ^ true,
            [condition, neg_power](const Element& e){ return condition(e) && (neg_power ^ !is_positive(e.second));}
        );

        if (denom_size > 0) {
            out << (latex ? "\\frac" : "");
        }

        if ((num_size > 1) && (denom_size > 0)) {
            out << (latex ? "" : KYLW "(" KRST);
        }
        out << num.str();
        if ((num_size > 1) && (denom_size > 0)) {
            out << (latex ? "" : KYLW ")" KRST);
        }

        if (denom_size > 0) {
            out << (latex ? "" : KGRN " / " KRST);
            if (denom_size > 1) {
                out << (latex ? "" : KYLW "(" KRST);
            }
            out << denom.str();
            if (denom_size > 1) {
                out << (latex ? "" : KYLW ")" KRST);
            }
        }

        is_simple = false;
        if ((num_size == 1) && (denom_size == 0)) {
            is_simple = true;
        }
    }

    size_t utf8_helper(const string& str, std::string::size_type idx) const {
        size_t cplen = 1;
        if ((str[idx] & 0xf8) == 0xf0) {
            cplen = 4;
        } else if ((str[idx] & 0xf0) == 0xe0) {
            cplen = 3;
        } else if ((str[idx] & 0xe0) == 0xc0) {
            cplen = 2;
        }
        if ((idx + cplen) > str.length()) {
            cplen = 1;
        }
        return cplen;
    }

    bool nat_cmp(const string& a, const string& b) const {
        std::string::size_type ia = 0, ib = 0;
        bool in_nb = false;

        while((ia < a.size()) && (ib < b.size())) {
            if (in_nb) {
                uint64_t nata = 0, natb = 0;
                while ((ia < a.size()) && (utf8_helper(a, ia) == 1) && isdigit(a[ia])) {
                    nata = (10*nata) + (a[ia] - '0');
                    ia++;
                }
                while ((ib < b.size()) && (utf8_helper(b, ib) == 1) && isdigit(b[ib])) {
                    natb = (10*natb) + (b[ib] - '0');
                    ib++;
                }
                if (nata != natb) {
                    return nata < natb;
                }
                in_nb = false;
            } else {
                size_t wa = utf8_helper(a, ia);
                size_t wb = utf8_helper(b, ib);

                if ((wa == 1) && isdigit(a[ia]) && (wb == 1) && isdigit(b[ib])) {
                    in_nb = true;
                    continue;
                }

                if (wa != wb) {
                    return wa < wb;
                }
                for(size_t i=0; i<min(wa,wb); i++) {
                    if (a[ia+i] != b[ib+i]) {
                        return a[ia+i] < b[ib+i];
                    }
                }

                ia += wa;
                ib += wb;
            }
        }

        if (ia >= a.size()) {
            return false;
        } else if (ib >= b.size()) {
            return true;
        }
        return false; /* Equal */
    }

    class RelationSummary {
    public:
        Rational zeta = Rational(0);
        Rational mu = Rational(0);
        Rational sigma = Rational(0);
        Rational sigma_prime = Rational(0);
        Rational theta = Rational(0);
        Rational jordan_t = Rational(0);
        Rational liouville = Rational(0);
        Rational tauk = Rational(0);

        RelationSummary(const std::vector<Element>& elements) {
            for (auto& element: elements) {
                const FormulaNode* name = element.first.get();
                if (name->isZeta()) {
                    zeta += abs(element.second);
                } else if (name->isLFuncNonZeta()) {
                    HFormula h_product = name->getLFuncProduct();
                    const FormulaNode* product = h_product.get();
                    for (size_t idx=0; idx<product->getProductSize(); idx++) {
                        const FormulaNode* formula = product->getProductElem(idx).get();
                        size_t n_times = 1;
                        if (formula->isPower()) {
                            n_times = formula->getPower();
                            formula = formula->getPowerInner().get();
                        }
                        assert(formula->isLeaf());
                        Rational to_add = Rational(n_times) * abs(element.second);
                        if (formula->isMu()) {
                            mu += to_add;
                        } else if (formula->isSigma()) {
                            sigma += to_add;
                        } else if (formula->isSigmaPrime()) {
                            sigma_prime += to_add;
                        } else if (formula->isTheta()) {
                            theta += to_add;
                        } else if (formula->isJordanT()) {
                            jordan_t += to_add;
                        } else if (formula->isLiouville()) {
                            liouville += to_add;
                        } else if (formula->isTauK()) {
                            tauk += to_add;
                        }
                    }
                }
            }
        }

        friend std::ostream& operator << (std::ostream& out, const RelationSummary& s) {
            out << "[RL: ζ:" << s.zeta << " µ:" << s.mu
                << " σ:" << s.sigma << " σ':" << s.sigma_prime << " θ:" << s.theta
                << " J:" << s.jordan_t << " λ:" << s.liouville << " τ:" << s.tauk
                << "]";
            return out;
        }

        typedef std::function<bool (const RelationSummary&)> instance_early_bailout;
        static bool no_early_bailout([[maybe_unused]] const RelationSummary&) {
            return true;
        }
    };

    class SymbolicInstantiation {
    public:
        typedef std::unordered_map<std::string, SomeInt> assignment;
        assignment variables;
        std::set<const Element*> relation_elements;
        std::set<const Element*> formula_elements;

        SymbolicInstantiation(const std::vector<Element>& relation, const std::vector<Element>& formula) {
            for (auto& element: relation) {
                bool is_in = relation_elements.find(&element) != relation_elements.end();
                assert(!is_in); /* Do we need multiset? */
                relation_elements.insert(&element);
            }
            for (auto& element: formula) {
                bool is_in = formula_elements.find(&element) != formula_elements.end();
                assert(!is_in); /* Do we need multiset? */
                formula_elements.insert(&element);
            }
        }
    };

public:
    Relation(const vector<pair<size_t, Rational>>& relation_row, vector<HFormula>& names) {
        for(const pair<size_t, Rational>& coeff: relation_row) {
            if(!(coeff.second == Rational(0))) {
                elements.push_back(make_pair(names[coeff.first], coeff.second));
            }
        }
        normalize();
    }

    Relation(const vector<pair<HFormula, Rational>>& relation_elements) {
        for(auto relation_element: relation_elements) {
            elements.push_back(relation_element);
        }
        normalize();
    }

    std::ostream& print(std::ostream& out, bool latex) const {
        std::ostringstream left;
        bool is_simple;
        print_frac(left, latex, false, [](const Element& e){ return !e.first.get()->isZeta(); }, is_simple);

        if (known) {
            if (latex) {
                std::string color = (known_formula[0] == 'C') ? "DeepSkyBlue" : "DarkGreen";
                out << "\\textcolor{" << color << "}{\\texttt{[" << known_formula << "]}}";
            } else {
                std::ostringstream ss;
                std::string color = (known_formula[0] == 'C') ? KCYN : KGRN;
                out << color << "[" << std::left << std::setw(4) << known_formula << "]" KRST;
            }
        } else {
            if (is_simple) {
                if (latex) {
                    out << "\\textcolor{FireBrick}{\\texttt{[!!!!]}}";
                } else {
                    out << KBLD KRED "[!!!!]" KRST;
                }
            } else {
                if (latex) {
                    out << "\\textcolor{Purple}{\\texttt{[****]}}";
                } else {
                    out << KMAG "[****]" KRST;
                }
            }

        }
        out << (latex ? "~~ " : " ");

        if (latex) {
            out << "$";
        }

        /* Print all non-zeta */
        out << left.str();

        out << (latex ? "=" : KRED " = " KRST);

        /* Then print all zeta */
        print_frac(out, latex, true , [](const Element& e){ return  e.first.get()->isZeta(); }, is_simple);

        if (latex) {
            out << "$" << endl << endl;
        }
        out.flush(); /* Better for debug purposes */
        return out;
    }

    friend std::ostream& operator << (std::ostream& out, const Relation &r) {
        return r.print(out, false);
    }

    bool operator < (const Relation& other) const {
        if (known != other.known) {
            return known > other.known;
        }
        if (known && (known_formula != other.known_formula)) {
            if (known_formula[0] != other.known_formula[0]) {
                return known_formula[0] > other.known_formula[0]; /* Hack */
            }
            return nat_cmp(known_formula, other.known_formula);
        }
        std::ostringstream ss, other_ss;
        ss << *this;
        other_ss << other;
        return nat_cmp(ss.str(), other_ss.str());
    }

private:
    void classify_raw(std::string formula) {
        known_formula = formula;
        known = true;
    }

    vector<string> instantiate_split_helper(const string& str, const string& delimiter) const {
        vector<string> tokens;
        size_t prev = 0;
        size_t pos = 0;
        do
        {
            pos = str.find(delimiter, prev);
            if (pos == string::npos) {
                pos = str.length();
            }
            string token = str.substr(prev, pos-prev);
            tokens.push_back(token);
            prev = pos + delimiter.length();
        }
        while ((pos < str.length()) && (prev < str.length()));
        return tokens;
    }

    bool instantiate_eval_helper(FormulaNode::Symbolic symbolic,
                                 const SymbolicInstantiation& instantiation,
                                 SomeInt* out_sum,
                                 unsigned int* out_uninstanticated_var_cnt,
                                 string* out_uninstanticated_var_name,
                                 SomeInt* out_uninstanticated_var_times,
                                 int debug) const {
        vector<string> n_symbols = instantiate_split_helper(symbolic.str, "+");
        unsigned int uninstanticated_var_cnt = 0;
        string uninstanticated_var_name;
        SomeInt uninstanticated_var_times = 0;
        SomeInt sum = 0;
        for (auto n_symbol: n_symbols) {
            vector<string> tmp = instantiate_split_helper(n_symbol, "*");
            SomeInt times;
            string name;
            if (tmp.size() == 1) {
                if (tmp[0][0] == '-') {
                    times = -1;
                    name = tmp[0].substr(1);
                } else {
                    times = 1;
                    name = tmp[0];
                }
            } else {
                assert(tmp.size() == 2);
                char c = tmp[0][0];
                if ((c >= '0' && c <= '9') || c == '-') {
                    /* Case: -4*x */
                    times = std::stoi(tmp[0]);
                } else {
                    /* Case: k*x */
                    auto it = instantiation.variables.find(tmp[0]);
                    if (it != instantiation.variables.end()) {
                        times = (it->second);
                    } else {
                        /* We do not handle when the first part has not been instantiated. Do better someday? */
                        if (debug>=0) { cerr << string(debug, ' ') << __func__ << " I " KGRY "UPR" KRST << endl; }
                        return false;
                    }
                }
                name = tmp[1];
            }

            /* Safety: check the name contains a single alpha character */
            assert((name.length() == 1) && (isalpha(name[0])));

            auto it = instantiation.variables.find(name);
            if (it != instantiation.variables.end()) {
                sum += times * (it->second);
            } else {
                if (uninstanticated_var_cnt > 0) {
                    /* We cannot instantiate more than one, bail out for now. Do better someday? */
                    if (debug>=0) { cerr << string(debug, ' ') << __func__ << " I " KGRY "TMV" KRST << endl; }
                    return false;
                }
                uninstanticated_var_name = name;
                uninstanticated_var_times = times;
                uninstanticated_var_cnt++;
            }
        }

        if (out_sum != NULL) {
            *out_sum = sum;
        }
        if (out_uninstanticated_var_cnt != NULL) {
            *out_uninstanticated_var_cnt = uninstanticated_var_cnt;
        }
        if (out_uninstanticated_var_name != NULL) {
            *out_uninstanticated_var_name = uninstanticated_var_name;
        }
        if (out_uninstanticated_var_times != NULL) {
            *out_uninstanticated_var_times = uninstanticated_var_times;
        }
        return true;
    }


    int iincr(int debug) const {
        if (debug >= 0) {
            return debug+1;
        }
        return debug;
    }

    bool try_instantiate(FormulaNode::MaybeSymbolic should_be_value, FormulaNode::MaybeSymbolic maybe_symbolic,
                         SymbolicInstantiation& instantiation, int debug) const {
        if (debug>=0) { cerr << string(debug, ' ') << __func__ << " I " << should_be_value << " wrt " << maybe_symbolic << endl; }
        SomeInt value;
        if (should_be_value.is_symbolic()) {
            /* We only accept values on the left side. Fully-instantiated symbolic is OK. */
            FormulaNode::Symbolic symbolic = should_be_value.extract_symbol();
            unsigned int uninstanticated_var_cnt = 0;
            SomeInt sum = 0;
            bool good = instantiate_eval_helper(symbolic, instantiation, &sum, &uninstanticated_var_cnt,
                                                NULL, NULL, iincr(debug));
            if (!good || (uninstanticated_var_cnt > 0)) {
                return false;
            }
            value = sum;
        } else {
            value = should_be_value.extract_value();
        }
        if (!maybe_symbolic.is_symbolic()) {
            return value == maybe_symbolic.extract_value();
        }
        FormulaNode::Symbolic symbolic = maybe_symbolic.extract_symbol();
        unsigned int uninstanticated_var_cnt = 0;
        string uninstanticated_var_name;
        SomeInt uninstanticated_var_times = 0;
        SomeInt sum = 0;
        bool good = instantiate_eval_helper(symbolic, instantiation, &sum, &uninstanticated_var_cnt,
                                            &uninstanticated_var_name, &uninstanticated_var_times, iincr(debug));
        if (!good) {
            return false;
        }
        if (uninstanticated_var_cnt > 0) {
            assert(uninstanticated_var_cnt == 1);
            /* See if we can match the existing value */
            SomeInt remaining = value - sum;
            if (remaining % uninstanticated_var_times == SomeInt(0)) {
                SomeInt newvar_value = remaining / uninstanticated_var_times;
                std::pair<std::string, SomeInt> new_elem (uninstanticated_var_name, newvar_value);
                instantiation.variables.insert(new_elem);
                if (debug>=0) {
                    cerr << string(debug, ' ') << __func__ << " I "
                         << KBLD << uninstanticated_var_name << ":=" << newvar_value << KRST << endl;
                }
                sum += uninstanticated_var_times*newvar_value;
                assert(sum == value);
            }
        }
        if (sum != value) {
            if (debug>=0) { cerr << string(debug, ' ') << __func__ << " I " KGRY "Badsum: " << sum << KRST << endl; }
            return false;
        }
        return true; /* O RLY? */
    }

    bool try_instantiate_product(const vector<HFormula>& v_rel, const vector<HFormula>& v_form,
                                 SymbolicInstantiation& instantiation, int debug) const {
            if((v_rel.size() == 0) && (v_form.size() == 0)) {
                return true; /* Nothing (a.k.a everything) matches */
            }
            if (debug>=0) { cerr << string(debug, ' ') << __func__ << " TP " << v_rel.size() << " wrt " << v_form.size() << endl; }
            for(size_t i=0; i<v_rel.size(); i++) {
                for(size_t j=0; j<v_form.size(); j++) {
                    SymbolicInstantiation new_instantiation = instantiation;
                    bool good = try_instantiate(v_rel[i], v_form[j], new_instantiation, iincr(debug));
                    if (good) {
                        vector<HFormula> new_v_rel = v_rel;
                        new_v_rel.erase(new_v_rel.begin() + i);
                        vector<HFormula> new_v_form = v_form;
                        new_v_form.erase(new_v_form.begin() + j);
                        good = try_instantiate_product(new_v_rel, new_v_form, new_instantiation, iincr(debug));
                        if (good) {
                            instantiation = new_instantiation;
                            return true;
                        }
                    }
                }
            }
            return false;
    }

    bool try_instantiate(const HFormula& h_rel, const HFormula& h_form, SymbolicInstantiation& instantiation, int debug) const {
        if (debug>=0) { cerr << string(debug, ' ') << __func__ << " TX " << h_rel << " wrt " << h_form << endl; }

        const FormulaNode* rel = h_rel.get();
        const FormulaNode* form = h_form.get();
        if (rel->isLeaf() && form->isLeaf()) {
            if (debug>=0) { cerr << string(debug, ' ') << __func__ << " TX L<->L" << endl; }
            if (!rel->isLeafSameAs(form)) {
                return false;
            }
            SymbolicInstantiation new_instantiation = instantiation;
            bool good = try_instantiate(rel->getLeafKS_dangerous(), form->getLeafKS_dangerous(), new_instantiation, iincr(debug));
            if (good) {
                good = try_instantiate(rel->getLeafLS_dangerous(), form->getLeafLS_dangerous(), new_instantiation, iincr(debug));
                if (good) {
                    instantiation = new_instantiation;
                }
            }
            return good;
        } else if (rel->isLeaf() && form->isPower()) {
            if (debug>=0) { cerr << string(debug, ' ') << __func__ << " TX L<->P" << endl; }
            if ((!form->getPowerS().is_symbolic()) && (form->getPowerS().extract_value() == 1)) {
                return try_instantiate(h_rel, form->getPowerInner(), instantiation, iincr(debug));
            }
        } else if (rel->isPower() && form->isLeaf()) {
            if (debug>=0) { cerr << string(debug, ' ') << __func__ << " TX P<->L" << endl; }
            if (rel->getPower() == 1) {
                return try_instantiate(rel->getPowerInner(), h_form, instantiation, iincr(debug));
            }
        } else if (rel->isPower() && form->isPower()) {
            SymbolicInstantiation new_instantiation = instantiation;
            bool good = try_instantiate(rel->getPower(), form->getPowerS(), new_instantiation, iincr(debug));
            if (good) {
                good = try_instantiate(rel->getPowerInner(), form->getPowerInner(), new_instantiation, iincr(debug));
                if (good) {
                    instantiation = new_instantiation;
                }
            }
            return good;
        } else if (rel->isProduct() && form->isProduct()) {
            if ((rel->getProductSize() == 1) && (form->getProductSize() == 1)) {
                return try_instantiate(rel->getProductElem(0), form->getProductElem(0), instantiation, iincr(debug));
            }
            /* Try all combinations... is there a better way? */
            vector<HFormula> rel_product, form_product;
            for (size_t idx=0; idx<rel->getProductSize(); idx++) {
                HFormula formula = rel->getProductElem(idx);
                size_t n_times = 1;
                if (formula.get()->isPower()) {
                    n_times = formula.get()->getPower();
                    formula = formula.get()->getPowerInner();
                }
                for (size_t times=0; times<n_times; times++) {
                    rel_product.push_back(formula);
                }
            }
            for (size_t idx=0; idx<form->getProductSize(); idx++) {
                HFormula formula = form->getProductElem(idx);
                size_t n_times = 1;
                if (formula.get()->isPower()) {
                    n_times = formula.get()->getPower();
                    formula = formula.get()->getPowerInner();
                }
                for (size_t times=0; times<n_times; times++) {
                    form_product.push_back(formula);
                }
            }
            if(rel_product.size() != form_product.size()) {
                return false; /* Since we flatten the product, it should have the same size */
            }
            return try_instantiate_product(rel_product, form_product, instantiation, iincr(debug));
        } else if (rel->isZeta() && form->isZeta()) {
            return try_instantiate(rel->getZetaExponentS(), form->getZetaExponentS(), instantiation, iincr(debug));
        } else if (rel->isLFuncNonZeta() && form->isLFuncNonZeta()) {
            SymbolicInstantiation new_instantiation = instantiation;
            bool good1 = try_instantiate(rel->getLFuncExponentS(), form->getLFuncExponentS(), new_instantiation, iincr(debug));
            bool good2 = try_instantiate(rel->getLFuncProduct(), form->getLFuncProduct(), new_instantiation, iincr(debug));
            if (good2 && !good1) {
                /* See e.g. C13. Else we always fall in the "more than one non-instantiated var, which is not handled above */
                good1 = try_instantiate(rel->getLFuncExponentS(), form->getLFuncExponentS(), new_instantiation, iincr(debug));
            }
            bool good = good1 && good2;
            if (good) {
                instantiation = new_instantiation;
            }
            return good;
        }

        if (debug>=0) { cerr << string(debug, ' ') << __func__ << " TX " KGRY "fail" KRST << endl; }
        return false; /* We could not do it. This does not mean this is not doable */
    }

    bool try_instantiate_e(Element rel, const Element& form, SymbolicInstantiation& instantiation, Rational& rel_power, int debug) const {
        if (debug>=0) {
            cerr << string(debug, ' ') << __func__ << " TI ";
            print_element(cerr, false, rel, false);
            cerr << " wrt ";
            print_element(cerr, false, form, false);
            cerr << endl;
        }
        if (rel.second == form.second) {
            /* OK */
        } else if (is_positive(rel.second) && is_positive(form.second) && (rel.second > form.second)) {
            /* We might have ζ(4)^2 w.r.t ζ(a)ζ(b), so me must match partial powers on 'rel' side. */
            rel.second = form.second;
        } else if (!is_positive(rel.second) && !is_positive(form.second) && (rel.second < form.second)) {
            /* Same as above but for negative powers */
            rel.second = form.second;
        } else {
            if (debug>=0) { cerr << string(debug, ' ') << __func__ << KGRY " TI DIFF-PWR" KRST << endl; }
            return false;
        }

        bool good = try_instantiate(rel.first, form.first, instantiation, iincr(debug));
        if (good) {
            rel_power = rel.second;
        }
        return good;
    }

    bool is_instance_of(SymbolicInstantiation& instantiation, int debug) const {
        if (debug>=0) cerr << string(debug, ' ') << __func__ << " IIO*" << endl;
        if (instantiation.relation_elements.size() == 0) {
            if (instantiation.formula_elements.size() == 0) {
                return true; /* YA RLY! */
            } else {
                /* It matches if every variable has been instantiated and the remaining product equals 1 */
                std::vector<Element> up, down;
                Rational sum = 0;
                for (const Element* element: instantiation.formula_elements) {
                    sum += element->second;
                    if (is_positive(element->second)) {
                        up.push_back(*element);
                    } else {
                        Element new_elem = *element;
                        new_elem.second = -new_elem.second;
                        down.push_back(new_elem);
                    }
                }
                if ((up.size() == 0) || (down.size() == 0)) {
                    return false;
                }
                if (sum == Rational(0)) {
                    SymbolicInstantiation new_instantiation = SymbolicInstantiation(up, down);
                    new_instantiation.variables = instantiation.variables;
                    if (is_instance_of(new_instantiation, iincr(debug))) {
                        /* Check if no new variable was introduced, i.e. there was no non-instantiated variable */
                        return new_instantiation.variables.size() == instantiation.variables.size();
                    }
                }
                return false;
            }
        }
        assert(instantiation.formula_elements.size() > 0);

        for (auto form_elem: instantiation.formula_elements) {
            for (auto rel_elem: instantiation.relation_elements) {
                SymbolicInstantiation new_instantiation = instantiation;
                Rational rel_power;
                if (try_instantiate_e(*rel_elem, *form_elem, new_instantiation, rel_power, iincr(debug))) {
                    Element new_rel_elem = *rel_elem;
                    new_instantiation.relation_elements.erase(rel_elem);
                    if (rel_power != rel_elem->second) {
                        new_rel_elem.second = new_rel_elem.second - rel_power;
                        new_instantiation.relation_elements.insert(&new_rel_elem);
                    }
                    new_instantiation.formula_elements.erase(form_elem);
                    if (is_instance_of(new_instantiation, iincr(debug))) {
                        instantiation = new_instantiation;
                        return true; /* NO WAI!! */
                    }
                }
            }
        }

        return false; /* Could not find a valid instantiation */
    }

    bool is_instance_of(RelationSummary::instance_early_bailout early_bailout, const RelationSummary& summary,
                        const Relation& formula, SymbolicInstantiation::assignment* assignment, int debug) {
        if (debug>=0) cerr << string(debug, ' ') << __func__ << " IIO " << *this << " wrt "<< formula << endl;
        if (formula.elements.size() < elements.size()) {
            return false;
        }
        if (!early_bailout(summary)) {
            return false;
        }
        SymbolicInstantiation instantiation = SymbolicInstantiation(elements, formula.elements);
        bool good = is_instance_of(instantiation, iincr(debug));
        if (good && (assignment != NULL)) {
#if 0
            cerr << KRED << "Assignment:";
            for (auto it = instantiation.variables.begin(); it != instantiation.variables.end(); it++) {
                cerr << " " << it->first << ':' << it->second << " ";
            }
            cerr << KRST << endl;
#endif
            *assignment = instantiation.variables;
        }
        return good;
    }

    typedef bool(Relation::*relation_classifier)(const RelationSummary&, string&);

    bool check_D3(const RelationSummary& summary, string& out_name) {
        std::string name = "D-3";
        vector<pair<HFormula, Rational>> vect{
            {HFormulaLFunction(HFormulaProduct(HFormulaLeaf(
                FormulaNode::LEAF_TAUK, (FormulaNode::LeafExtraArg){.k = 2, .l = 0})),
                FormulaNode::Symbolic("s")), Rational(1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("s")), Rational(-2)},
        };
        Relation formula = Relation(vect);
        formula.classify_raw(name);
        bool good = is_instance_of(RelationSummary::no_early_bailout, summary, formula, NULL, -1);
        out_name = name;
        return good;
    }

    bool check_D6(const RelationSummary& summary, string& out_name) {
        std::string name = "D-6";
        vector<pair<HFormula, Rational>> vect{
            {HFormulaLFunction(HFormulaProduct(HFormulaLeaf(
                FormulaNode::LEAF_SIGMA, (FormulaNode::LeafExtraArg){.k = FormulaNode::Symbolic("k"), .l = 0})),
                FormulaNode::Symbolic("s")), Rational(1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("s")), Rational(-1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("s+-k")), Rational(-1)},
        };
        Relation formula = Relation(vect);
        formula.classify_raw(name);
        SymbolicInstantiation::assignment assignment;
        bool good = is_instance_of(RelationSummary::no_early_bailout, summary, formula, &assignment, -1);
        if (!good) {
            return false;
        }
        out_name = name;
        if (assignment["k"] == 1) {
            out_name = "D-4";
        }
        return good;
    }

    bool check_D10(const RelationSummary& summary, string& out_name) {
        std::string name = "D-10";
        vector<pair<HFormula, Rational>> vect{
            {HFormulaLFunction(HFormulaProduct(HFormulaPower(
                HFormulaLeaf(FormulaNode::LEAF_MU, (FormulaNode::LeafExtraArg){.k = 1, .l = 0}), 2)), FormulaNode::Symbolic("s")), Rational(1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("2*s")), Rational(1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("s")), Rational(-1)},
        };
        Relation formula = Relation(vect);
        formula.classify_raw(name);
        bool good = is_instance_of(RelationSummary::no_early_bailout, summary, formula, NULL, -1);
        out_name = name;
        return good;
    }

    bool check_D11(const RelationSummary& summary, string& out_name) {
        std::string name = "D-11";
        vector<pair<HFormula, Rational>> vect{
            {HFormulaLFunction(HFormulaProduct(
                HFormulaPower(HFormulaLeaf(FormulaNode::LEAF_TAUK, (FormulaNode::LeafExtraArg){.k = 2, .l = 0}), 2)),
                FormulaNode::Symbolic("s")), Rational(1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("s")), Rational(-4)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("2*s")), Rational(1)},
        };
        Relation formula = Relation(vect);
        formula.classify_raw(name);
        bool good = is_instance_of(RelationSummary::no_early_bailout, summary, formula, NULL, -1);
        out_name = name;
        return good;
    }

    bool check_D12(const RelationSummary& summary, string& out_name) {
        (void) summary;
        (void) out_name;
        // TODO ! (also, rm D3)
        return false;
    }

    bool check_D15(const RelationSummary& summary, string& out_name) {
        std::string name = "D-15";
        vector<pair<HFormula, Rational>> vect{
            {HFormulaLFunction(HFormulaProduct(HFormulaLeaf(
                FormulaNode::LEAF_ZETAK, (FormulaNode::LeafExtraArg){.k = FormulaNode::Symbolic("k"), .l = 0})), FormulaNode::Symbolic("s")), Rational(1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("s+-k")), Rational(-1)},
        };
        Relation formula = Relation(vect);
        formula.classify_raw(name);
        SymbolicInstantiation::assignment assignment;
        bool good = is_instance_of(RelationSummary::no_early_bailout, summary, formula, &assignment, -1);
        if (!good) {
            return false;
        }
        out_name = name;
        if (assignment["k"] == 1) {
            out_name = "D-13";
        }
        return good;
    }

    bool check_D18(const RelationSummary& summary, string& out_name) {
        std::string name = "D-18";
        vector<pair<HFormula, Rational>> vect{
            {HFormulaLFunction(HFormulaProduct(HFormulaLeaf(FormulaNode::LEAF_THETA)), FormulaNode::Symbolic("s")), Rational(1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("2*s")), Rational(1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("s")), Rational(-2)},
        };
        Relation formula = Relation(vect);
        formula.classify_raw(name);
        bool good = is_instance_of(RelationSummary::no_early_bailout, summary, formula, NULL, -1);
        out_name = name;
        return good;
    }

    bool check_D21(const RelationSummary& summary, string& out_name) {
        std::string name = "D-21";
        vector<pair<HFormula, Rational>> vect{
            {HFormulaLFunction(HFormulaProduct(HFormulaLeaf(FormulaNode::LEAF_LIOUVILLE), HFormulaLeaf(FormulaNode::LEAF_THETA)),
                               FormulaNode::Symbolic("s")), Rational(1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("2*s")), Rational(-1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("s")), Rational(2)},
        };
        Relation formula = Relation(vect);
        formula.classify_raw(name);
        bool good = is_instance_of(RelationSummary::no_early_bailout, summary, formula, NULL, -1);
        out_name = name;
        return good;
    }

    bool check_D22(const RelationSummary& summary, string& out_name) {
        std::string name = "D-22";
        vector<pair<HFormula, Rational>> vect{
            {HFormulaLFunction(HFormulaProduct(HFormulaLeaf(
                FormulaNode::LEAF_MU, (FormulaNode::LeafExtraArg){.k = FormulaNode::Symbolic("k"), .l = 0})), FormulaNode::Symbolic("s")), Rational(1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("k*s")), Rational(1)},
        };
        Relation formula = Relation(vect);
        formula.classify_raw(name);
        SymbolicInstantiation::assignment assignment;
        bool good = is_instance_of(RelationSummary::no_early_bailout, summary, formula, &assignment, -1);
        out_name = name;
        if (assignment["k"] == 1) {
            out_name = "D-2";
        }
        return good;
    }

    bool check_D24(const RelationSummary& summary, string& out_name) {
        std::string name = "D-24";
                vector<pair<HFormula, Rational>> vect{
            {HFormulaLFunction(HFormulaProduct(
                HFormulaLeaf(FormulaNode::LEAF_XI, (FormulaNode::LeafExtraArg){.k = FormulaNode::Symbolic("k"), .l = 0})),
                FormulaNode::Symbolic("s")), Rational(1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("s")), Rational(-1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("k*s")), Rational(1)},
        };
        Relation formula = Relation(vect);
        formula.classify_raw(name);
        bool good = is_instance_of(RelationSummary::no_early_bailout, summary, formula, NULL, -1);
        out_name = name;
        return good;
    }

    bool check_D25(const RelationSummary& summary, string& out_name) {
        std::string name = "D-25";
                vector<pair<HFormula, Rational>> vect{
            {HFormulaLFunction(HFormulaProduct(
                HFormulaLeaf(FormulaNode::LEAF_LIOUVILLE),
                HFormulaLeaf(FormulaNode::LEAF_SIGMA_PRIME, (FormulaNode::LeafExtraArg){.k = FormulaNode::Symbolic("k"), .l = 0})),
                FormulaNode::Symbolic("s")), Rational(1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("s")), Rational(1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("s+-k")), Rational(-1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("2*s")), Rational(-1)},
        };
        Relation formula = Relation(vect);
        formula.classify_raw(name);
        bool good = is_instance_of(RelationSummary::no_early_bailout, summary, formula, NULL, -1);
        out_name = name;
        return good;
    }

    bool check_D26(const RelationSummary& summary, string& out_name) {
        std::string name = "D-26";
                vector<pair<HFormula, Rational>> vect{
            {HFormulaLFunction(HFormulaProduct(HFormulaLeaf(
                FormulaNode::LEAF_PSI, (FormulaNode::LeafExtraArg){.k = FormulaNode::Symbolic("k"), .l = 0})), FormulaNode::Symbolic("s")), Rational(1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("2*s")), Rational(1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("s+-k")), Rational(-1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("s")), Rational(-1)},
        };
        Relation formula = Relation(vect);
        formula.classify_raw(name);
        bool good = is_instance_of(RelationSummary::no_early_bailout, summary, formula, NULL, -1);
        out_name = name;
        return good;
    }

    /* Warning: formula in the catalog is wrong! */
    bool check_D27(const RelationSummary& summary, string& out_name) {
        std::string name = "D-27";
        vector<pair<HFormula, Rational>> vect{
            {HFormulaLFunction(HFormulaProduct(
                HFormulaLeaf(FormulaNode::LEAF_LIOUVILLE),
                HFormulaLeaf(FormulaNode::LEAF_JORDAN_T, (FormulaNode::LeafExtraArg){.k = FormulaNode::Symbolic("k"), .l = 0})
                ), FormulaNode::Symbolic("s")), Rational(1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("s")), Rational(-1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("2*s+-2*k")), Rational(-1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("s+-k")), Rational(1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("2*s")), Rational(1)},
        };
        Relation formula = Relation(vect);
        formula.classify_raw(name);
        bool good = is_instance_of(RelationSummary::no_early_bailout, summary, formula, NULL, -1);
        out_name = name;
        return good;
    }

    bool check_D28(const RelationSummary& summary, string& out_name) {
        std::string name = "D-28";
                vector<pair<HFormula, Rational>> vect{
            {HFormulaLFunction(HFormulaProduct(HFormulaLeaf(
                FormulaNode::LEAF_NU, (FormulaNode::LeafExtraArg){.k = FormulaNode::Symbolic("k"), .l = 0})), FormulaNode::Symbolic("s")), Rational(1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("k*s")), Rational(-1)},
        };
        Relation formula = Relation(vect);
        formula.classify_raw(name);
        bool good = is_instance_of(RelationSummary::no_early_bailout, summary, formula, NULL, -1);
        out_name = name;
        return good;
    }


    bool check_D30(const RelationSummary& summary, string& out_name) {
        std::string name = "D-30";
                vector<pair<HFormula, Rational>> vect{
            {HFormulaLFunction(HFormulaProduct(HFormulaLeaf(
                FormulaNode::LEAF_RHO, (FormulaNode::LeafExtraArg){.k = FormulaNode::Symbolic("k"), .l = FormulaNode::Symbolic("t")})), FormulaNode::Symbolic("s")), Rational(1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("s+-k")), Rational(-1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("t*s")), Rational(-1)},
        };
        Relation formula = Relation(vect);
        formula.classify_raw(name);
        bool good = is_instance_of(RelationSummary::no_early_bailout, summary, formula, NULL, -1);
        out_name = name;
        return good;
    }

    bool check_D37(const RelationSummary& summary, string& out_name) {
        std::string name = "D-37";
                vector<pair<HFormula, Rational>> vect{
            {HFormulaLFunction(HFormulaProduct(HFormulaLeaf(
                FormulaNode::LEAF_SIGMA_PRIME, (FormulaNode::LeafExtraArg){.k = FormulaNode::Symbolic("k"), .l = 0})), FormulaNode::Symbolic("s")), Rational(1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("2*s+-2*k")), Rational(-1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("s")), Rational(-1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("s+-k")), Rational(1)},
        };
        Relation formula = Relation(vect);
        formula.classify_raw(name);
        bool good = is_instance_of(RelationSummary::no_early_bailout, summary, formula, NULL, -1);
        out_name = name;
        return good;
    }

    bool check_D38(const RelationSummary& summary, string& out_name) {
        std::string name = "D-38";
                vector<pair<HFormula, Rational>> vect{
            {HFormulaLFunction(HFormulaProduct(
                HFormulaLeaf(FormulaNode::LEAF_LIOUVILLE),
                HFormulaLeaf(FormulaNode::LEAF_SIGMA_PRIME, (FormulaNode::LeafExtraArg){.k = FormulaNode::Symbolic("h"), .l = 0}),
                HFormulaLeaf(FormulaNode::LEAF_SIGMA, (FormulaNode::LeafExtraArg){.k = FormulaNode::Symbolic("k"), .l = 0})),
                FormulaNode::Symbolic("s")), Rational(1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("2*s+-h+-k")), Rational(-1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("2*s+-2*k")), Rational(-1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("s+-h+-k")), Rational(-1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("s+-h")), Rational(-1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("2*s")), Rational(-1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("s")), Rational(1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("s+-k")), Rational(1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("4*s+-2*h+-2*k")), Rational(1)},
        };
        Relation formula = Relation(vect);
        formula.classify_raw(name);

        RelationSummary::instance_early_bailout early_bailout = [](const RelationSummary&r) {
            return (r.liouville == Rational(1)) && (r.sigma_prime == Rational(1)) && (r.sigma == Rational(1))
                && ((r.zeta == Rational(8)) || (r.zeta == Rational(6)));
        };

        bool good = is_instance_of(early_bailout, summary, formula, NULL, -1);
        out_name = name;
        return good;
    }

    bool check_D39(const RelationSummary& summary, string& out_name) {
        std::string name = "D-39";
                vector<pair<HFormula, Rational>> vect{
            {HFormulaLFunction(HFormulaProduct(HFormulaLeaf(
                FormulaNode::LEAF_SIGMA_PRIME, (FormulaNode::LeafExtraArg){.k = FormulaNode::Symbolic("h"), .l = 0}),
                HFormulaLeaf(FormulaNode::LEAF_SIGMA, (FormulaNode::LeafExtraArg){.k = FormulaNode::Symbolic("k"), .l = 0})), FormulaNode::Symbolic("s")), Rational(1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("2*s+-k+-h")), Rational(-1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("2*s+-2*h+-2*k")), Rational(-1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("s+-k")), Rational(-1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("2*s+-2*h")), Rational(-1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("s")), Rational(-1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("s+-h")), Rational(1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("s+-k+-h")), Rational(1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("4*s+-2*h+-2*k")), Rational(1)},
        };
        Relation formula = Relation(vect);
        formula.classify_raw(name);

        RelationSummary::instance_early_bailout early_bailout = [](const RelationSummary&r) {
            return (r.sigma == Rational(1)) && (r.sigma_prime == Rational(1))
                && ((r.zeta == Rational(8)) || (r.zeta == Rational(6)));
        };

        bool good = is_instance_of(early_bailout, summary, formula, NULL, -1);
        out_name = name;
        return good;
    }

    bool check_D40(const RelationSummary& summary, string& out_name) {
        std::string name = "D-40";
                vector<pair<HFormula, Rational>> vect{
            {HFormulaLFunction(HFormulaProduct(
                HFormulaLeaf(FormulaNode::LEAF_LIOUVILLE),
                HFormulaLeaf(FormulaNode::LEAF_SIGMA_PRIME, (FormulaNode::LeafExtraArg){.k = FormulaNode::Symbolic("h"), .l = 0}),
                HFormulaLeaf(FormulaNode::LEAF_SIGMA, (FormulaNode::LeafExtraArg){.k = FormulaNode::Symbolic("k"), .l = 0})),
                FormulaNode::Symbolic("s")), Rational(1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("2*s+-k+-h")), Rational(-1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("s+-h")), Rational(-1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("2*s+-2*k")), Rational(-1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("2*s")), Rational(-1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("s")), Rational(1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("4*s+-2*h+-2*k")), Rational(1)},
        };
        Relation formula = Relation(vect);
        formula.classify_raw(name);

        RelationSummary::instance_early_bailout early_bailout = [](const RelationSummary&r) {
            return (r.liouville == Rational(1)) && (r.sigma_prime == Rational(1)) && (r.sigma == Rational(1))
                && (r.zeta == Rational(6));
        };

        bool good = is_instance_of(early_bailout, summary, formula, NULL, -1);
        out_name = name;
        return good;
    }

    bool check_D41(const RelationSummary& summary, string& out_name) {
        std::string name = "D-41";
                vector<pair<HFormula, Rational>> vect{
            {HFormulaLFunction(HFormulaProduct(
                HFormulaLeaf(FormulaNode::LEAF_SIGMA_PRIME, (FormulaNode::LeafExtraArg){.k = FormulaNode::Symbolic("h"), .l = 0}),
                HFormulaLeaf(FormulaNode::LEAF_SIGMA_PRIME, (FormulaNode::LeafExtraArg){.k = FormulaNode::Symbolic("k"), .l = 0})),
                FormulaNode::Symbolic("s")), Rational(1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("s+-k+-h")), Rational(-1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("2*s+-2*h")), Rational(-1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("2*s+-2*k")), Rational(-1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("s")), Rational(-1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("2*s+-h+-k")), Rational(1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("s+-h")), Rational(1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("s+-k")), Rational(1)},
        };
        Relation formula = Relation(vect);
        formula.classify_raw(name);

        RelationSummary::instance_early_bailout early_bailout = [](const RelationSummary&r) {
            return (r.sigma_prime == Rational(2))
                && ((r.zeta == Rational(7)) || (r.zeta == Rational(5)));
        };

        bool good = is_instance_of(early_bailout, summary, formula, NULL, -1);
        out_name = name;
        return good;
    }

    bool check_D42(const RelationSummary& summary, string& out_name) {
        std::string name = "D-42";
                vector<pair<HFormula, Rational>> vect{
            {HFormulaLFunction(HFormulaProduct(
                HFormulaLeaf(FormulaNode::LEAF_LIOUVILLE),
                HFormulaLeaf(FormulaNode::LEAF_SIGMA_PRIME, (FormulaNode::LeafExtraArg){.k = FormulaNode::Symbolic("h"), .l = 0}),
                HFormulaLeaf(FormulaNode::LEAF_SIGMA_PRIME, (FormulaNode::LeafExtraArg){.k = FormulaNode::Symbolic("k"), .l = 0})),
                FormulaNode::Symbolic("s")), Rational(1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("s+-k")), Rational(-1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("s+-h")), Rational(-1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("2*s+-2*k+-2*h")), Rational(-1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("2*s")), Rational(-1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("2*s+-h+-k")), Rational(1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("s+-h+-k")), Rational(1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("s")), Rational(1)},
        };
        Relation formula = Relation(vect);
        formula.classify_raw(name);

        RelationSummary::instance_early_bailout early_bailout = [](const RelationSummary&r) {
            return (r.liouville == Rational(1)) && (r.sigma_prime == Rational(2))
                && ((r.zeta == Rational(7)) || (r.zeta == Rational(5)));
        };

        bool good = is_instance_of(early_bailout, summary, formula, NULL, -1);
        out_name = name;
        return good;
    }

    bool check_D46(const RelationSummary& summary, string& out_name) {
        std::string name = "D-46";
        vector<pair<HFormula, Rational>> vect{
            {HFormulaLFunction(HFormulaProduct(
                HFormulaLeaf(FormulaNode::LEAF_LIOUVILLE),
                HFormulaLeaf(FormulaNode::LEAF_SIGMA, (FormulaNode::LeafExtraArg){.k = FormulaNode::Symbolic("h"), .l = 0}),
                HFormulaLeaf(FormulaNode::LEAF_SIGMA, (FormulaNode::LeafExtraArg){.k = FormulaNode::Symbolic("k"), .l = 0})
                ), FormulaNode::Symbolic("s")), Rational(1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("2*s+-2*h+-2*k")), Rational(-1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("2*s+-2*k")), Rational(-1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("2*s+-2*h")), Rational(-1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("2*s")), Rational(-1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("s")), Rational(1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("s+-h+-k")), Rational(1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("s+-k")), Rational(1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("s+-h")), Rational(1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("2*s+-h+-k")), Rational(1)},
        };
        Relation formula = Relation(vect);
        formula.classify_raw(name);

        RelationSummary::instance_early_bailout early_bailout = [](const RelationSummary&r) {
            return (r.liouville == Rational(1)) && (r.sigma == Rational(2))
                && ((r.zeta == Rational(9)) || (r.zeta == Rational(7)) || (r.zeta == Rational(5)));
        };

        bool good = is_instance_of(early_bailout, summary, formula, NULL, -1);
        out_name = name;
        return good;
    }

    bool check_D47(const RelationSummary& summary, string& out_name) {
        std::string name = "D-47";
        vector<pair<HFormula, Rational>> vect{
            {HFormulaLFunction(HFormulaProduct(
                HFormulaLeaf(FormulaNode::LEAF_LIOUVILLE),
                HFormulaLeaf(FormulaNode::LEAF_SIGMA, (FormulaNode::LeafExtraArg){.k = FormulaNode::Symbolic("k"), .l = 0})
                ), FormulaNode::Symbolic("s")), Rational(1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("2*s")), Rational(-1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("2*s+-2*k")), Rational(-1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("s")), Rational(1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("s+-k")), Rational(1)},
        };
        Relation formula = Relation(vect);
        formula.classify_raw(name);

        RelationSummary::instance_early_bailout early_bailout = [](const RelationSummary&r) {
            return (r.liouville == Rational(1)) && (r.sigma == Rational(1));
        };

        bool good = is_instance_of(early_bailout, summary, formula, NULL, -1);
        out_name = name;
        return good;
    }

    bool check_D50(const RelationSummary& summary, string& out_name) {
        std::string name = "D-50";
        vector<pair<HFormula, Rational>> vect{
            {HFormulaLFunction(HFormulaProduct(
                HFormulaLeaf(FormulaNode::LEAF_LIOUVILLE),
                HFormulaPower(HFormulaLeaf(FormulaNode::LEAF_TAUK, (FormulaNode::LeafExtraArg){.k = 2, .l = 0}), 2)
                ), FormulaNode::Symbolic("s")), Rational(1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("2*s")), Rational(-3)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("s")), Rational(4)},
        };
        Relation formula = Relation(vect);
        formula.classify_raw(name);
        bool good = is_instance_of(RelationSummary::no_early_bailout, summary, formula, NULL, -1);
        out_name = name;
        return good;
    }

    bool check_D52(const RelationSummary& summary, string& out_name) {
        std::string name = "D-52";
        vector<pair<HFormula, Rational>> vect{
            {HFormulaLFunction(HFormulaProduct(HFormulaLeaf(
                FormulaNode::LEAF_JORDAN_T, (FormulaNode::LeafExtraArg){.k = FormulaNode::Symbolic("k"), .l = 0})), FormulaNode::Symbolic("s")), Rational(1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("s+-k")), Rational(-1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("s")), Rational(1)},
        };
        Relation formula = Relation(vect);
        formula.classify_raw(name);
        SymbolicInstantiation::assignment assignment;
        bool good = is_instance_of(RelationSummary::no_early_bailout, summary, formula, &assignment, -1);
        if (!good) {
            return false;
        }
        out_name = name;
        if (assignment["k"] == 1) {
            out_name = "D-5";
        }
        return good;
    }

    bool check_D53(const RelationSummary& summary, string& out_name) {
        std::string name = "D-53";
        vector<pair<HFormula, Rational>> vect{
            {HFormulaLFunction(HFormulaProduct(
                HFormulaPower(HFormulaLeaf(FormulaNode::LEAF_LIOUVILLE), 1)),
                FormulaNode::Symbolic("s")), Rational(1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("2*s")), Rational(-1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("s")), Rational(1)},
        };
        Relation formula = Relation(vect);
        formula.classify_raw(name);
        bool good = is_instance_of(RelationSummary::no_early_bailout, summary, formula, NULL, -1);
        out_name = name;
        return good;
    }

    bool check_D58(const RelationSummary& summary, string& out_name) {
        std::string name = "D-58";
        vector<pair<HFormula, Rational>> vect{
            {HFormulaLFunction(HFormulaProduct(
                HFormulaLeaf(FormulaNode::LEAF_SIGMA, (FormulaNode::LeafExtraArg){.k = FormulaNode::Symbolic("a"), .l = 0}),
                HFormulaLeaf(FormulaNode::LEAF_SIGMA, (FormulaNode::LeafExtraArg){.k = FormulaNode::Symbolic("b"), .l = 0})
                ), FormulaNode::Symbolic("s")), Rational(1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("s")), Rational(-1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("s+-a")), Rational(-1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("s+-b")), Rational(-1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("s+-a+-b")), Rational(-1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("2*s+-a+-b")), Rational(1)},
        };
        Relation formula = Relation(vect);
        formula.classify_raw(name);

        RelationSummary::instance_early_bailout early_bailout = [](const RelationSummary&r) {
            return (r.sigma == Rational(2)) && (r.zeta == Rational(5));
        };

        bool good = is_instance_of(early_bailout, summary, formula, NULL, -1);
        out_name = name;
        return good;
    }

    bool check_C1(const RelationSummary& summary, string& out_name) {
        std::string name = "C-1";
        vector<pair<HFormula, Rational>> vect{
            {HFormulaLFunction(HFormulaProduct(
                HFormulaLeaf(FormulaNode::LEAF_SIGMA, (FormulaNode::LeafExtraArg){.k = FormulaNode::Symbolic("s"), .l = 0}),
                HFormulaPower(HFormulaLeaf(FormulaNode::LEAF_MU, (FormulaNode::LeafExtraArg){.k = 1, .l = 0}), 2)
                ), FormulaNode::Symbolic("2*s")), Rational(1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("s")), Rational(-1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("3*s")), Rational(1)},
        };
        Relation formula = Relation(vect);
        formula.classify_raw(name);
        bool good = is_instance_of(RelationSummary::no_early_bailout, summary, formula, NULL, -1);
        out_name = name;
        return good;
    }

    bool check_C11(const RelationSummary& summary, string& out_name) {
        std::string name = "C-11";
        vector<pair<HFormula, Rational>> vect{
            {HFormulaLFunction(HFormulaProduct(
                HFormulaLeaf(FormulaNode::LEAF_JORDAN_T, (FormulaNode::LeafExtraArg){.k = FormulaNode::Symbolic("s"), .l = 0}),
                HFormulaLeaf(FormulaNode::LEAF_MU, (FormulaNode::LeafExtraArg){.k = 1, .l = 0})
                ), FormulaNode::Symbolic("2*s")), Rational(1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("2*s")), Rational(-1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("3*s")), Rational(-1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("s")), Rational(1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("6*s")), Rational(1)},
        };
        Relation formula = Relation(vect);
        formula.classify_raw(name);
        bool good = is_instance_of(RelationSummary::no_early_bailout, summary, formula, NULL, -1);
        out_name = name;
        return good;
    }

    bool check_C13(const RelationSummary& summary, string& out_name) {
        std::string name = "C-13";
        vector<pair<HFormula, Rational>> vect{
            {HFormulaLFunction(HFormulaProduct(
                HFormulaLeaf(FormulaNode::LEAF_JORDAN_T, (FormulaNode::LeafExtraArg){.k = FormulaNode::Symbolic("i"), .l = 0}),
                HFormulaLeaf(FormulaNode::LEAF_SIGMA, (FormulaNode::LeafExtraArg){.k = FormulaNode::Symbolic("i"), .l = 0})
                ), FormulaNode::Symbolic("2*i+k")), Rational(1)},
            {HFormulaLFunction(HFormulaProduct(
                HFormulaLeaf(FormulaNode::LEAF_JORDAN_T, (FormulaNode::LeafExtraArg){.k = FormulaNode::Symbolic("2*i"), .l = 0}),
                HFormulaLeaf(FormulaNode::LEAF_SIGMA, (FormulaNode::LeafExtraArg){.k = FormulaNode::Symbolic("i+k"), .l = 0})
                ), FormulaNode::Symbolic("3*i+2*k")), Rational(-1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("2*i+2*k")), Rational(-1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("i+2*k")), Rational(1)},
        };
        Relation formula = Relation(vect);
        formula.classify_raw(name);
        bool good = is_instance_of(RelationSummary::no_early_bailout, summary, formula, NULL, -1);
        out_name = name;
        return good;
    }

    bool check_C14(const RelationSummary& summary, string& out_name) {
        std::string name = "C-14";
        vector<pair<HFormula, Rational>> vect{
            {HFormulaLFunction(HFormulaProduct(
                HFormulaLeaf(FormulaNode::LEAF_THETA),
                HFormulaLeaf(FormulaNode::LEAF_SIGMA, (FormulaNode::LeafExtraArg){.k = FormulaNode::Symbolic("s"), .l = 0})
                ), FormulaNode::Symbolic("2*s")), Rational(1)},
            {HFormulaLFunction(HFormulaProduct(
                HFormulaPower(HFormulaLeaf(FormulaNode::LEAF_JORDAN_T, (FormulaNode::LeafExtraArg){.k = FormulaNode::Symbolic("s"), .l = 0}), 2)
                ), FormulaNode::Symbolic("4*s")), Rational(-1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("s")), Rational(-2)},
        };
        Relation formula = Relation(vect);
        formula.classify_raw(name);
        bool good = is_instance_of(RelationSummary::no_early_bailout, summary, formula, NULL, -1);
        out_name = name;
        return good;
    }

    bool check_C15(const RelationSummary& summary, string& out_name) {
        std::string name = "C-15";
        vector<pair<HFormula, Rational>> vect{
            {HFormulaLFunction(HFormulaProduct(
                HFormulaLeaf(FormulaNode::LEAF_JORDAN_T, (FormulaNode::LeafExtraArg){.k = FormulaNode::Symbolic("s"), .l = 0}),
                HFormulaLeaf(FormulaNode::LEAF_SIGMA, (FormulaNode::LeafExtraArg){.k = FormulaNode::Symbolic("s"), .l = 0})
                ), FormulaNode::Symbolic("3*s")), Rational(1)},
            {HFormulaLFunction(HFormulaProduct(
                HFormulaPower(HFormulaLeaf(FormulaNode::LEAF_JORDAN_T, (FormulaNode::LeafExtraArg){.k = FormulaNode::Symbolic("s"), .l = 0}), 2)
                ), FormulaNode::Symbolic("4*s")), Rational(1)},
            {HFormulaLFunction(HFormulaProduct(
                HFormulaPower(HFormulaLeaf(FormulaNode::LEAF_JORDAN_T, (FormulaNode::LeafExtraArg){.k = FormulaNode::Symbolic("s"), .l = 0}), 2),
                HFormulaPower(HFormulaLeaf(FormulaNode::LEAF_SIGMA, (FormulaNode::LeafExtraArg){.k = FormulaNode::Symbolic("s"), .l = 0}), 2)
                ), FormulaNode::Symbolic("5*s")), Rational(-1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("2*s")), Rational(-1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("3*s")), Rational(1)},
        };
        Relation formula = Relation(vect);
        formula.classify_raw(name);
        bool good = is_instance_of(RelationSummary::no_early_bailout, summary, formula, NULL, -1);
        out_name = name;
        return good;
    }

    bool check_C17(const RelationSummary& summary, string& out_name) {
        std::string name = "C-17";
        vector<pair<HFormula, Rational>> vect{
            {HFormulaLFunction(HFormulaProduct(
                HFormulaLeaf(FormulaNode::LEAF_LIOUVILLE),
                HFormulaLeaf(FormulaNode::LEAF_TAUK, (FormulaNode::LeafExtraArg){.k = 2, .l = 0})),
                FormulaNode::Symbolic("s")), Rational(1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("2*s")), Rational(-2)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("s")), Rational(2)},
        };
        Relation formula = Relation(vect);
        formula.classify_raw(name);
        bool good = is_instance_of(RelationSummary::no_early_bailout, summary, formula, NULL, -1);
        out_name = name;
        return good;
    }

    bool check_C18(const RelationSummary& summary, string& out_name) {
        std::string name = "C-18";
        vector<pair<HFormula, Rational>> vect{
            {HFormulaLFunction(HFormulaProduct(
                HFormulaLeaf(FormulaNode::LEAF_TAUK, (FormulaNode::LeafExtraArg){.k = 2, .l = 0}),
                HFormulaLeaf(FormulaNode::LEAF_SIGMA, (FormulaNode::LeafExtraArg){.k = FormulaNode::Symbolic("k"), .l = 0})
                ), FormulaNode::Symbolic("s")), Rational(1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("s+-k")), Rational(-2)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("s")), Rational(-2)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("2*s+-k")), Rational(1)},
        };
        Relation formula = Relation(vect);
        formula.classify_raw(name);
        bool good = is_instance_of(RelationSummary::no_early_bailout, summary, formula, NULL, -1);
        out_name = name;
        return good;
    }

    bool check_C19(const RelationSummary& summary, string& out_name) {
        std::string name = "C-19";
        vector<pair<HFormula, Rational>> vect{
            {HFormulaLFunction(HFormulaProduct(
                HFormulaLeaf(FormulaNode::LEAF_LIOUVILLE),
                HFormulaLeaf(FormulaNode::LEAF_TAUK, (FormulaNode::LeafExtraArg){.k = 2, .l = 0}),
                HFormulaLeaf(FormulaNode::LEAF_SIGMA, (FormulaNode::LeafExtraArg){.k = FormulaNode::Symbolic("k"), .l = 0})
                ), FormulaNode::Symbolic("s")), Rational(1)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("2*s+-2*k")), Rational(-2)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("2*s")), Rational(-2)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("s+-k")), Rational(2)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("s")), Rational(2)},
            {HFormulaLFunction(HFormulaOne(), FormulaNode::Symbolic("2*s+-k")), Rational(1)},
        };
        Relation formula = Relation(vect);
        formula.classify_raw(name);
        bool good = is_instance_of(RelationSummary::no_early_bailout, summary, formula, NULL, -1);
        out_name = name;
        return good;
    }


    vector<relation_classifier> classifiers{
        /*
            * D formulae,
            * from Gould, H. W., & Shonhiwa, T. (2008). A catalog of interesting Dirichlet series.
            * Missouri Journal of Mathematical Sciences, 20(1), 2-18.
            */

        /* D-1 we won't find */
        /* Check for D-2: handled in D-22 */
        &Relation::check_D3,
        /* Check for D-4: handled in D-6 */
        /* Check for D-5: handled in D-52 */
        &Relation::check_D6,
        /* D-7 & D-8 we won't find */
        /* Check for D-9: this is D-18 */
        &Relation::check_D10,
        &Relation::check_D11,
        &Relation::check_D12,                     // TODO
        /* Check for D-13: handled in D-15 */
        /* D-14 we won't find */
        &Relation::check_D15,
        /* D-16 & D-17 we won't find */
        &Relation::check_D18,
        /* D-19 & D-20 we won't find */
        &Relation::check_D21,
        &Relation::check_D22,
        /* D-23 PHI_K not implemented */          // TODO
        &Relation::check_D24,
        &Relation::check_D25,
        &Relation::check_D26,
        &Relation::check_D27,
        &Relation::check_D28,
        /* D-29 form is not handled by checker */ // TODO
        &Relation::check_D30,
        /* D-31 we won't find */

        /* ... */

        /* D-35 and D-36 we won't find */
        &Relation::check_D37,
        &Relation::check_D38,
        &Relation::check_D39,
        &Relation::check_D40,
        &Relation::check_D41,
        &Relation::check_D42,

        /* ... */

        /* D-44 and D-45 we won't find */
        &Relation::check_D46,
        &Relation::check_D47,
        /* D-48 we won't find */

        /* ... */

        &Relation::check_D50,

        /* ... */

        &Relation::check_D52,
        &Relation::check_D53,
        /* D-54 we won't find */

        /* ... */

        /* D-56 & D-57 we won't find */
        &Relation::check_D58,
        /* D-59 and later we won't find */

        /*
         * C formulae, found with CrazySums
         */
        &Relation::check_C1,
        &Relation::check_C11,
        &Relation::check_C13,
        &Relation::check_C14,
        &Relation::check_C15,
        &Relation::check_C17,
        &Relation::check_C18,
        &Relation::check_C19,
    };

public:
    void classify() {
        const RelationSummary summary(elements);
        #if DEBUG_TIME_CLASSIFY
        size_t classify_debug_idx = 0;
        #endif /* DEBUG_TIME_CLASSIFY */
        for (auto classifier: classifiers) {
            #if DEBUG_TIME_CLASSIFY
            auto t1 = std::chrono::high_resolution_clock::now();
            #endif /* DEBUG_TIME_CLASSIFY */
            std::string found_formula;
            if ((this->*classifier)(summary, found_formula)) {
                known_formula = found_formula;
                known = true;
            }
            #if DEBUG_TIME_CLASSIFY
            auto t2 = std::chrono::high_resolution_clock::now();
            std::chrono::duration<float> e21 = t2 - t1;
            if (classify_debug_idx >= relation_time_classify_debug.size()) {
                assert(classify_debug_idx == relation_time_classify_debug.size());
                relation_time_classify_debug.push_back(e21);
            } else {
                relation_time_classify_debug[classify_debug_idx] += e21;
            }
            classify_debug_idx++;
            #endif /* DEBUG_TIME_CLASSIFY */
            if (known) {
                return;
            }
        }
    }
};
