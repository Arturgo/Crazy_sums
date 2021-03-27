#pragma once
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include "print.h"

class Relation {
private:
    typedef std::pair<const FormulaName*, Rational> Element;
    std::vector<Element> elements;
    bool known = false;
    std::string known_formula;

    void normalize() {
        /* To be called once at the end of the constructor. Improves printing and formula detection */
        /* Try to normalize the L-side. We want the shortest product on bottom, then the bigger exponent on bottom. */
        size_t l_top = 0, l_bottom = 0;
        const FormulaName* top_n = NULL;
        const FormulaName* bottom_n = NULL;
        for (const Element& e: elements) {
            if (e.first->isLFunc() && !e.first->isZeta()) {
                if (is_positive(e.second)) {
                    l_top++;
                    top_n = e.first;
                } else {
                    l_bottom++;
                    bottom_n = e.first;
                }
            }
        }
        if ((l_top < l_bottom) ||
            ((l_top == 1) && (l_bottom == 1) && (top_n->getLFuncExponent() > bottom_n->getLFuncExponent()))) {
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
        element.first->print_inside(out, latex);
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
        typedef std::function<bool (const RelationSummary&)> instance_early_bailout;
        static bool no_early_bailout([[maybe_unused]] const RelationSummary&) {
            return false;
        }
    };

    class SymbolicInstantiation {
    public:
        typedef std::unordered_map<std::string, SomeInt> assignment;
        assignment variables;
        std::unordered_set<const Element*> relation_elements;
        std::unordered_set<const Element*> formula_elements;

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
    Relation(const vector<pair<size_t, Rational>>& relation_row, vector<const FormulaName*>& names) {
        for(const pair<size_t, Rational>& coeff: relation_row) {
            if(!(coeff.second == Rational(0))) {
                elements.push_back(make_pair(names[coeff.first], coeff.second));
            }
        }
        normalize();
    }

    Relation(const vector<pair<const FormulaName*, Rational>>& relation_elements) {
        for(auto relation_element: relation_elements) {
            elements.push_back(relation_element);
        }
        normalize();
    }

    Relation(const vector<Rational>& relation_row, vector<const FormulaName*>& names,
             const vector<size_t>& iCol_in_rows) {
        for(size_t iCol = 0;iCol < relation_row.size();iCol++) {
            if(!(relation_row[iCol] == Rational(0))) {
                elements.push_back(make_pair(names[iCol_in_rows[iCol]], relation_row[iCol]));
            }
        }
        normalize();
    }

    std::ostream& print(std::ostream& out, bool latex) const {
        std::ostringstream left;
        bool is_simple;
        print_frac(left, latex, false, [](const Element& e){ return !e.first->isZeta(); }, is_simple);

        if (known) {
            if (latex) {
                out << "\\textcolor{DarkGreen}{\\texttt{[" << known_formula << "]}}";
            } else {
                std::ostringstream ss;
                out << KGRN << "[" << known_formula << "]" KRST;
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
        print_frac(out, latex, true , [](const Element& e){ return  e.first->isZeta(); }, is_simple);

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

    int iincr(int debug) const {
        if (debug >= 0) {
            return debug+1;
        }
        return debug;
    }

    bool try_instantiate(SomeInt value, FormulaName::MaybeSymbolic maybe_symbolic, SymbolicInstantiation& instantiation, int debug) const {
        if (debug>=0) { cerr << string(debug, ' ') << __func__ << " I " << value << " wrt " << maybe_symbolic << endl; }
        if (!maybe_symbolic.is_symbolic()) {
            return value == maybe_symbolic.extract_value();
        }
        FormulaName::Symbolic symbolic = maybe_symbolic.extract_symbol();
        vector<string> n_symbols = instantiate_split_helper(symbolic.str, "+");
        int uninstanticated_var_cnt = 0;
        string uninstanticated_var_name;
        SomeInt uninstanticated_var_times = 0;
        SomeInt sum = 0;
        for (auto n_symbol: n_symbols) {
            vector<string> tmp = instantiate_split_helper(n_symbol, "*");
            SomeInt times;
            string name;
            if (tmp.size() == 1) {
                times = 1;
                name = tmp[0];
            } else {
                assert(tmp.size() == 2);
                times = std::stoi(tmp[0]);
                name = tmp[1];
            }

            auto it = instantiation.variables.find(name);
            if (it != instantiation.variables.end()) {
                sum += times * (it->second);
            } else {
                if (uninstanticated_var_cnt > 0) {
                    /* TODO: we cannot instantiate more than one, bail out for now */
                    if (debug>=0) { cerr << string(debug, ' ') << __func__ << " I " KGRY "TMV" KRST << endl; }
                    return false;
                }
                uninstanticated_var_name = name;
                uninstanticated_var_times = times;
                uninstanticated_var_cnt++;
            }
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

    bool try_instantiate(const FormulaName* rel, const FormulaName* form, SymbolicInstantiation& instantiation, int debug) const {
        if (debug>=0) { cerr << string(debug, ' ') << __func__ << " TX " << *rel << " wrt " << *form << endl; }

        if (rel->isLeaf() && form->isLeaf()) {
            if (debug>=0) { cerr << string(debug, ' ') << __func__ << " TX L<->L" << endl; }
            if (!rel->isLeafSameAs(form)) {
                return false;
            }
            SymbolicInstantiation new_instantiation = instantiation;
            bool good = try_instantiate(rel->getLeafK_dangerous(), form->getLeafKS_dangerous(), new_instantiation, iincr(debug));
            if (good) {
                good = try_instantiate(rel->getLeafL_dangerous(), form->getLeafLS_dangerous(), new_instantiation, iincr(debug));
                if (good) {
                    instantiation = new_instantiation;
                }
            }
            return good;
        } else if (rel->isLeaf() && form->isPower()) {
            if (debug>=0) { cerr << string(debug, ' ') << __func__ << " TX L<->P" << endl; }
            if ((!form->getPowerS().is_symbolic()) && (form->getPowerS().extract_value() == 1)) {
                return try_instantiate(rel, form->getPowerInner(), instantiation, iincr(debug));
            }
        } else if (rel->isPower() && form->isLeaf()) {
            if (debug>=0) { cerr << string(debug, ' ') << __func__ << " TX P<->L" << endl; }
            if (rel->getPower() == 1) {
                return try_instantiate(rel->getPowerInner(), form, instantiation, iincr(debug));
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
            return false; //TODO: try all combinations...
        } else if (rel->isZeta() && form->isZeta()) {
            return try_instantiate(rel->getZetaExponent(), form->getZetaExponentS(), instantiation, iincr(debug));
        } else if (rel->isLFuncNonZeta() && form->isLFuncNonZeta()) {
            SymbolicInstantiation new_instantiation = instantiation;
            bool good = try_instantiate(rel->getLFuncExponent(), form->getLFuncExponentS(), new_instantiation, iincr(debug));
            if (good) {
                good = try_instantiate(rel->getLFuncProduct(), form->getLFuncProduct(), new_instantiation, iincr(debug));
                if (good) {
                    instantiation = new_instantiation;
                }
            }
            return good;
        }

        if (debug>=0) { cerr << string(debug, ' ') << __func__ << " TX " KGRY "fail" KRST << endl; }
        return false; /* We could not do it. This does not mean this is not doable */
    }

    bool try_instantiate(const Element& rel, const Element& form, SymbolicInstantiation& instantiation, int debug) const {
        if (debug>=0) {
            cerr << string(debug, ' ') << __func__ << " TI ";
            print_element(cerr, false, rel, false);
            cerr << " wrt ";
            print_element(cerr, false, form, false);
            cerr << endl;
        }
        if (rel.second != form.second) {
            if (debug>=0) { cerr << string(debug, ' ') << __func__ << KBLD " TI DIFF-PWR" KRST << endl; }
            return false; /* Not the same relation power... maybe we will need to do better here someday */
        }
        return try_instantiate(rel.first, form.first, instantiation, iincr(debug));
    }

    bool is_instance_of(SymbolicInstantiation& instantiation, int debug) const {
        if (debug>=0) cerr << string(debug, ' ') << __func__ << " IIO*" << endl;
        if (instantiation.relation_elements.size() == 0) {
            if (instantiation.formula_elements.size() == 0) {
                return true; /* YA RLY! */
            } else {
                /* It matches if every variable has been instantiated and the remaining product equals 1 */
                return false; //TODO: handle me!
            }
        }
        assert(instantiation.formula_elements.size() > 0);

        for (auto rel_elem: instantiation.relation_elements) {
            for (auto form_elem: instantiation.formula_elements) {
                SymbolicInstantiation new_instantiation = instantiation;
                if (try_instantiate(*rel_elem, *form_elem, new_instantiation, iincr(debug))) {
                    new_instantiation.relation_elements.erase(rel_elem);
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
                        const Relation &formula, SymbolicInstantiation::assignment* assignment, int debug) {
        if (debug>=0) cerr << string(debug, ' ') << __func__ << " IIO " << *this << " wrt "<< formula << endl;
        if (formula.elements.size() < elements.size()) {
            return false;
        }
        if (early_bailout(summary)) {
            return false;
        }
        SymbolicInstantiation instantiation = SymbolicInstantiation(elements, formula.elements);
        bool good = is_instance_of(instantiation, iincr(debug));
        if (good && (assignment != NULL)) {
#if 0
            cerr << KRED << "Assignement:";
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

    bool check_D2(const RelationSummary& summary, string& out_name) {
        std::string name = "D-2 ";
        vector<pair<const FormulaName*, Rational>> vect{
            {new FormulaNameLFunction(new FormulaNameProduct(new FormulaNameLeaf(FormulaName::LEAF_MU)), FormulaName::Symbolic("s")), Rational(1)},
            {new FormulaNameLFunction(new FormulaName(), FormulaName::Symbolic("s")), Rational(1)},
        };
        Relation formula = Relation(vect);
        formula.classify_raw(name);
        bool good = is_instance_of(RelationSummary::no_early_bailout, summary, formula, NULL, -1);
        out_name = name;
        return good;
    }

    bool check_D4_D6(const RelationSummary& summary, string& out_name) {
        std::string name = "D-6 ";
        vector<pair<const FormulaName*, Rational>> vect{
            {new FormulaNameLFunction(new FormulaNameProduct(new FormulaNameLeaf(
                FormulaName::LEAF_SIGMA, (FormulaName::LeafExtraArg){.k = FormulaName::Symbolic("1*k"), .l = 0})), FormulaName::Symbolic("s")), Rational(1)},
            {new FormulaNameLFunction(new FormulaName(), FormulaName::Symbolic("1*s")), Rational(-1)},
            {new FormulaNameLFunction(new FormulaName(), FormulaName::Symbolic("1*s+-1*k")), Rational(-1)},
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
            out_name = "D-4 ";
        }
        return good;
    }

    bool check_D10(const RelationSummary& summary, string& out_name) {
        std::string name = "D-10";
        vector<pair<const FormulaName*, Rational>> vect{
            {new FormulaNameLFunction(new FormulaNameProduct(new FormulaNamePower(new FormulaNameLeaf(FormulaName::LEAF_MU), 2)), FormulaName::Symbolic("s")), Rational(1)},
            {new FormulaNameLFunction(new FormulaName(), FormulaName::Symbolic("2*s")), Rational(1)},
            {new FormulaNameLFunction(new FormulaName(), FormulaName::Symbolic("1*s")), Rational(-1)},
        };
        Relation formula = Relation(vect);
        formula.classify_raw(name);
        bool good = is_instance_of(RelationSummary::no_early_bailout, summary, formula, NULL, -1);
        out_name = name;
        return good;
    }

    bool check_D18(const RelationSummary& summary, string& out_name) {
        std::string name = "D-18";
        vector<pair<const FormulaName*, Rational>> vect{
            {new FormulaNameLFunction(new FormulaNameProduct(new FormulaNamePower(new FormulaNameLeaf(FormulaName::LEAF_THETA), 1)), FormulaName::Symbolic("s")), Rational(1)},
            {new FormulaNameLFunction(new FormulaName(), FormulaName::Symbolic("2*s")), Rational(1)},
            {new FormulaNameLFunction(new FormulaName(), FormulaName::Symbolic("1*s")), Rational(-2)},
        };
        Relation formula = Relation(vect);
        formula.classify_raw(name);
        bool good = is_instance_of(RelationSummary::no_early_bailout, summary, formula, NULL, -1);
        out_name = name;
        return good;
    }

    bool check_D52(const RelationSummary& summary, string& out_name) {
        std::string name = "D-52";
        vector<pair<const FormulaName*, Rational>> vect{
            {new FormulaNameLFunction(new FormulaNameProduct(new FormulaNameLeaf(
                FormulaName::LEAF_JORDAN_T, (FormulaName::LeafExtraArg){.k = FormulaName::Symbolic("1*k"), .l = 0})), FormulaName::Symbolic("s")), Rational(1)},
            {new FormulaNameLFunction(new FormulaName(), FormulaName::Symbolic("1*s+-1*k")), Rational(-1)},
            {new FormulaNameLFunction(new FormulaName(), FormulaName::Symbolic("1*s")), Rational(1)},
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
            out_name = "D-5 ";
        }
        return good;
    }

    bool check_D58(string& out_name) {
        bool debug = false;
        if (debug) cerr << __func__ << " A" << endl;
        int Lsigma_exponent = 0;
        int sigmaA_k = 0;
        int sigmaB_k = 0;
        vector<int> zeta_numbers_down;
        int zeta_number_up = 0;

        for(auto element: elements) {
            if (debug) cerr << __func__ << " L" << endl;
            const FormulaName* name = element.first;
            Rational power = element.second;
            if (name->isZeta()) {
                int number = name->getZetaExponent();
                if (!is_integer(power)) {
                    if (debug) cerr << __func__ << " B" << endl;
                    return false;
                }
                if (is_positive(power)) {
                    if (zeta_number_up != 0) {
                        if (debug) cerr << __func__ << " C" << endl;
                        return false;
                    }
                    zeta_number_up = number;
                    if (debug) cerr << __func__ << " zeta_number_up " << zeta_number_up << endl;
                }
                for (SomeInt i=0; i<-power.getNumerator(); i++) {
                    if (zeta_numbers_down.size() < 4) {
                        zeta_numbers_down.push_back(number);
                        if (debug) cerr << __func__ << " zeta_number_down[] " << number << endl;
                    } else {
                        if (debug) cerr << __func__ << " D" << endl;
                        return false;
                    }
                }
            } else {
                if (sigmaA_k != 0) {
                    if (debug) cerr << __func__ << " E" << endl;
                    return false;
                }
                if (power != Rational(1)) {
                    if (debug) cerr << __func__ << " F" << endl;
                    return false;
                }
                const FormulaName* infunc = name->getLFuncProduct();
                if (infunc->getProductSize() == 1) {
                    const FormulaName* inner = infunc->getProductElem(0);
                    if (!inner->isSigma()) {
                        if (debug) cerr << __func__ << " G" << endl;
                        return false;
                    } else if (inner->getPower() != 2) {
                        if (debug) cerr << __func__ << " H" << endl;
                        return false;
                    } else {
                        assert(sigmaA_k == 0);
                        assert(Lsigma_exponent == 0);
                        Lsigma_exponent = name->getLFuncExponent();
                        sigmaA_k = inner->getSigmaK();
                        sigmaB_k = inner->getSigmaK();
                    }
                } else if (infunc->getProductSize() == 2) {
                    for (size_t idx=0; idx<2; idx++) {
                        const FormulaName* inner = infunc->getProductElem(idx);
                        if (!inner->isSigma()) {
                            if (debug) cerr << __func__ << " I" << endl;
                            return false;
                        } else if (inner->getPower() != 1) {
                            if (debug) cerr << __func__ << " J" << endl;
                            return false;
                        } else {
                            if (sigmaA_k == 0) {
                                Lsigma_exponent = name->getLFuncExponent();
                                sigmaA_k = inner->getSigmaK();
                            } else {
                                assert(sigmaB_k == 0);
                                if(Lsigma_exponent != name->getLFuncExponent()) {
                                    if (debug) cerr << __func__ << " K" << endl;
                                    return false;
                                }
                                sigmaB_k = inner->getSigmaK();
                            }
                        }
                    }
                } else {
                    if (debug) cerr << __func__ << " L" << endl;
                    return false;
                }
            }
        }

        if ((Lsigma_exponent == 0) || (sigmaA_k == 0) || (sigmaB_k == 0)) {
            if (debug) cerr << __func__ << " M" << endl;
            return false;
        }
        if (!(
              ((zeta_numbers_down.size() == 3) && (zeta_number_up == 0))
              ||
              ((zeta_numbers_down.size() == 4) && (zeta_number_up != 0))
              )) {
            if (debug) cerr << __func__ << " N" << endl;
            return false;
        }
        std::sort(zeta_numbers_down.begin(), zeta_numbers_down.end());

        if (debug) cerr << __func__ << " Z " << Lsigma_exponent << " " << sigmaA_k << " " << sigmaB_k << endl;

        /* Tandem, the final mathemagical check */
        int a = sigmaA_k;
        int b = sigmaB_k;
        int s = Lsigma_exponent;
        vector<int> to_find_down { s, s-a, s-b, s-a-b };
        std::sort(to_find_down.begin(), to_find_down.end());
        int to_find_up = 2*s-a-b;

        if (to_find_up < 0) {
            if (debug) cerr << __func__ << " P" << endl;
            return false;
        }
        for(auto n: to_find_down) {
            if (n < 0) {
                if (debug) cerr << __func__ << " Q" << endl;
                return false;
            }
        }
        for (size_t idx=0; idx<to_find_down.size(); idx++) {
            if (to_find_down[idx] == to_find_up) {
                to_find_down.erase(to_find_down.begin() + idx);
                to_find_up = 0;
                break;
            }
        }

        if (to_find_up != zeta_number_up) {
            if (debug) cerr << __func__ << " R" << endl;
            return false;
        }

        if (to_find_down.size() != zeta_numbers_down.size()) {
            if (debug) cerr << __func__ << " S" << endl;
            return false;
        }
        for (size_t idx=0; idx<to_find_down.size(); idx++) {
            if (to_find_down[idx] != zeta_numbers_down[idx]) {
                if (debug) cerr << __func__ << " T" << endl;
                return false;
            }
        }

        /* Incredible */
        out_name = "D-58";
        return true;
    }


public:
    void classify() {
        size_t zeta = 0;
        size_t nb = elements.size();

        for(auto element: elements) {
            const FormulaName* name = element.first;
            if (!name->isLFunc()) {
                return; /* What is this??? */
            }
            if (name->isZeta()) {
                zeta++;
            }
        }

        vector<relation_classifier> classifiers{
            /* D-1 we won't find */
            &Relation::check_D2,
            /* Check for D-3 */ ///TODO: implement
            &Relation::check_D4_D6,
            /* Check for D-5: handled in D-52 */
            /* Check for D-6: handled in D-4 */
            /* D-7 & D-8 we won't find */
            /* Check for D-9: this is D-18 */
            &Relation::check_D10,

            /* ... */

            &Relation::check_D18,

            /* ... */

            &Relation::check_D52,

            /* ... */

            /* Check for D-58: TODO: handled separately below */
            /* D-59 and later we won't find */
        };

        const RelationSummary dummy; //TODO
        for (auto classifier: classifiers) {
            std::string found_formula;
            if ((this->*classifier)(dummy, found_formula)) {
                known_formula = found_formula;
                known = true;
                return;
            }
        }

        if (zeta + 1 == nb) {
            /* Check for D-58 */
            if (check_D58(known_formula)) {
                known = true;
                return;
            }
        }
    }
};
