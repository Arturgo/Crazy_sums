#pragma once
#include <iostream>
#include "print.h"

class Relation {
private:
    typedef std::pair<const FormulaName*, Rational> Element;
    std::vector<Element> elements;
    bool known = false;
    std::string known_formula;

    std::ostream& print_element(std::ostream& out, bool latex, Element& element, bool neg_power) const {
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

    void print_frac(std::ostream& out, bool latex, bool neg_power, print_condition condition) const {
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
    }

public:
    Relation(const vector<Rational>& relation_row, vector<const FormulaName*>& names,
             const vector<size_t>& iCol_in_rows) {
        for(size_t iCol = 0;iCol < relation_row.size();iCol++) {
            if(!(relation_row[iCol] == Rational(0))) {
                elements.push_back(make_pair(names[iCol_in_rows[iCol]], relation_row[iCol]));
            }
        }
    }

    std::ostream& print(std::ostream& out, bool latex) const {
        if (known) {
            if (latex) {
                out << "\\texttt{[" << known_formula << "]}";
            } else {
                std::ostringstream ss;
                out << KGRN << "[" << known_formula << "]" KRST;
            }
        } else {
            if (latex) {
                out << "\\texttt{[****]}";
            } else {
                out << KMAG << "[****]" KRST;
            }
        }
        out << "  ";

        if (latex) {
            out << "$";
        }

        /* Print all non-zeta */
        print_frac(out, latex, false, [](const Element& e){ return !e.first->isZeta(); });

        out << (latex ? "=" : KRED " = " KRST);

        /* Then print all zeta */
        print_frac(out, latex, true , [](const Element& e){ return  e.first->isZeta(); });

        if (latex) {
            out << "$" << endl << endl;
        }
        out.flush(); /* Better for debug purposes */
        return out;
    }

    friend std::ostream& operator << (std::ostream& out, const Relation &r) {
        return r.print(out, false);
    }

    bool check_D5(string& out_name) {
        bool debug = false;
        if (debug) cerr << __func__ << " A" << endl;
        int Lphi_exponent = 0;
        int zetaA_number = 0;
        int zetaB_number = 0;

        for(auto element: elements) {
            if (debug) cerr << __func__ << " L" << endl;
            const FormulaName* name = element.first;
            Rational power = element.second;
            if (name->isZeta()) {
                int number = name->getZetaExponent();
                if ((zetaA_number == 0) && (power == Rational(1))) {
                    zetaA_number = number;
                } else if ((zetaB_number == 0) && (power == Rational(-1))) {
                    zetaB_number = number;
                } else {
                    if (debug) cerr << __func__ << " B" << endl;
                    return false;
                }
            } else {
                if (power != Rational(1)) {
                    if (debug) cerr << __func__ << " C" << endl;
                    return false;
                }
                const FormulaName* infunc = name->getLFuncProduct();
                if (infunc->getProductSize() != 1) {
                    if (debug) cerr << __func__ << " D" << endl;
                    return false;
                }
                const FormulaName* inner = infunc->getProductElem(0);
                if (!inner->isPhi()) {
                    if (debug) cerr << __func__ << " E" << endl;
                    return false;
                } else if (inner->getPower() != 1) {
                    if (debug) cerr << __func__ << " F" << endl;
                    return false;
                } else if (Lphi_exponent != 0) {
                    if (debug) cerr << __func__ << " G" << endl;
                    return false;
                } else {
                    Lphi_exponent = name->getLFuncExponent();
                }
            }
        }

        if ((Lphi_exponent == 0) || (zetaA_number == 0) || (zetaB_number == 0)) {
            if (debug) cerr << __func__ << " H" << endl;
            return false;
        }
        if (debug) cerr << __func__ << " Z" << Lphi_exponent << " " << zetaA_number << " " << zetaB_number << endl;
        bool good = (Lphi_exponent == zetaA_number) && (Lphi_exponent == zetaB_number+1);
        if (!good) {
            return false;
        }
        out_name = "D-5 ";
        return true;
    }

    bool check_D6(string& out_name) {
        bool debug = false;
        if (debug) cerr << __func__ << " A" << endl;
        int Lsigma_exponent = 0;
        int sigma_k = 0;
        int zetaA_number = 0;
        int zetaB_number = 0;

        for(auto element: elements) {
            if (debug) cerr << __func__ << " L" << endl;
            const FormulaName* name = element.first;
            Rational power = element.second;
            if (name->isZeta()) {
                int number = name->getZetaExponent();
                if ((zetaA_number == 0) && (power == Rational(-1))) {
                    zetaA_number = number;
                } else if ((zetaB_number == 0) && (power == Rational(-1))) {
                    zetaB_number = number;
                } else {
                    if (debug) cerr << __func__ << " B" << endl;
                    return false;
                }
            } else {
                if (power != Rational(1)) {
                    if (debug) cerr << __func__ << " C" << endl;
                    return false;
                }
                const FormulaName* infunc = name->getLFuncProduct();
                if (infunc->getProductSize() != 1) {
                    if (debug) cerr << __func__ << " D" << endl;
                    return false;
                }
                const FormulaName* inner = infunc->getProductElem(0);
                if (!inner->isSigma()) {
                    if (debug) cerr << __func__ << " E" << endl;
                    return false;
                } else if (inner->getPower() != 1) {
                    if (debug) cerr << __func__ << " F" << endl;
                    return false;
                } else if (Lsigma_exponent != 0) {
                    if (debug) cerr << __func__ << " G" << endl;
                    return false;
                } else {
                    Lsigma_exponent = name->getLFuncExponent();
                    assert(sigma_k == 0);
                    sigma_k = inner->getSigmaK();
                }
            }
        }

        if ((Lsigma_exponent == 0) || (zetaA_number == 0) || (zetaB_number == 0)) {
            if (debug) cerr << __func__ << " H" << endl;
            return false;
        }
        if (zetaA_number > zetaB_number) {
            auto tmp = zetaA_number;
            zetaA_number = zetaB_number;
            zetaB_number = tmp;
        }

        if (debug) cerr << __func__ << " Z" << Lsigma_exponent << " " << zetaA_number << " " << zetaB_number << endl;
        bool good = (Lsigma_exponent == zetaA_number+sigma_k) && (Lsigma_exponent == zetaB_number);
        if (!good) {
            return false;
        }
        out_name = "D-6 ";
        if (sigma_k == 1) {
            out_name = "D-4 ";
        }
        return true;
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
                for (int i=0; i<-power.getNumerator(); i++) {
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

        if (zeta + 1 == nb) {
            /* Check for D-5 */
            if (check_D5(known_formula)) {
                known = true;
                return;
            }
            /* Check for D-6 */
            if (check_D6(known_formula)) {
                known = true;
                return;
            }
            /* Check for D-58 */
            if (check_D58(known_formula)) {
                known = true;
                return;
            }
        }
    }
};
