#pragma once
#include <cassert>
#include <memory>
#include <regex>

class HFormula {
public:
    class Node;

protected:
    std::shared_ptr<const Node> formula;
    HFormula() {}

public:
    const Node* get() const {
        return formula.get();
    }
};

class HFormula::Node {
public:
    class Symbolic {
    public:
        std::string str;

        Symbolic() {
            /* Nothing to do */
        }

        Symbolic(std::string value) {
            str = value;
        }

        friend std::ostream& operator << (std::ostream& out, const Symbolic &s) {
            out << s.str;
            return out;
        }
    };

    class MaybeSymbolic {
    private:
        bool _is_symbolic;
        int value;
        Symbolic symbolic;

    public:
        MaybeSymbolic() {
            _is_symbolic = false;
            value = 0;
        }

        MaybeSymbolic(int init_value) {
            _is_symbolic = false;
            value = init_value;
        }

        MaybeSymbolic(Symbolic init_value) {
            _is_symbolic = true;
            symbolic = init_value;
        }

        bool is_symbolic(void) const {
            return _is_symbolic;
        }

        int extract_value(void) const {
            assert(!is_symbolic());
            return value;
        }

        Symbolic extract_symbol(void) const {
            assert(is_symbolic());
            return symbolic;
        }

        friend std::ostream& operator << (std::ostream& out, const MaybeSymbolic &s) {
            if (s.is_symbolic()) {
                out << s.extract_symbol();
            } else {
                out << s.extract_value();
            }
            return out;
        }
    };

    enum LeafType { LEAF_THETA, LEAF_ZETAK, LEAF_MU, LEAF_NU, LEAF_RHO,
                    LEAF_SIGMA, LEAF_XI, LEAF_SIGMA_PRIME, LEAF_TAUK, LEAF_PSI,
                    LEAF_JORDAN_T, LEAF_LIOUVILLE,
                    LEAF_UNKNOWN
                };
    typedef struct LeafExtraArg {
        MaybeSymbolic k;
        MaybeSymbolic l;
    } LeafExtraArg;

protected:
    enum FormulaType { FORM_ONE, FORM_LEAF, FORM_LFUNC, FORM_POWER, FORM_PRODUCT };
    FormulaType formula_type = FORM_ONE;
    LeafType leaf_type;
    LeafExtraArg leaf_extra;
    MaybeSymbolic power = MaybeSymbolic(0);
    MaybeSymbolic exponent = MaybeSymbolic(0);
    std::vector<HFormula> sub_formula;

    std::string textify(LeafType in, bool latex) const {
        std::string ret;
        switch(in) {
            case LEAF_THETA:
                ret = (latex ? "\\theta{}" : "θ");
                break;
            case LEAF_LIOUVILLE:
                ret = (latex ? "\\lambda{}" : "λ");
                break;
            case LEAF_MU:
                ret = (latex ? "\\mu{}" : "µ");
                if (leaf_extra.k.is_symbolic()) {
                    ret += (latex ? "_{" : "_{") + leaf_extra.k.extract_symbol().str + (latex ? "}" : "}");
                } else {
                    auto value = leaf_extra.k.extract_value();
                    if (value != 1) {
                        ret += "_" + std::to_string(value);
                    }
                }
                break;
            case LEAF_SIGMA:
                ret = (latex ? "\\sigma{}" : "σ");
                if (leaf_extra.k.is_symbolic()) {
                    ret += (latex ? "_{" : "_{") + leaf_extra.k.extract_symbol().str + (latex ? "}" : "}");
                } else {
                    auto value = leaf_extra.k.extract_value();
                    if (value != 1) {
                        ret += "_" + std::to_string(value);
                    }
                }
                break;
            case LEAF_SIGMA_PRIME:
                ret = (latex ? "\\sigma{}'" : "σ'");
                if (leaf_extra.k.is_symbolic()) {
                    ret += (latex ? "_{" : "_{") + leaf_extra.k.extract_symbol().str + (latex ? "}" : "}");
                } else {
                    auto value = leaf_extra.k.extract_value();
                    if (value != 1) {
                        ret += "_" + std::to_string(value);
                    }
                }
                break;
            case LEAF_ZETAK:
                ret = (latex ? "\\zeta{}" : "ζ");
                if (leaf_extra.k.is_symbolic()) {
                    ret += (latex ? "_{" : "_{") + leaf_extra.k.extract_symbol().str + (latex ? "}" : "}");
                } else {
                    auto value = leaf_extra.k.extract_value();
                    ret += "_" + std::to_string(value);
                }
                break;
            case LEAF_JORDAN_T:
                if (leaf_extra.k.is_symbolic()) {
                    ret = (latex ? "J" : "J");
                    ret += (latex ? "_{" : "_{") + leaf_extra.k.extract_symbol().str + (latex ? "}" : "}");
                } else {
                    auto value = leaf_extra.k.extract_value();
                    if (value == 1) {
                        ret = (latex ? "\\phi{}" : "φ");
                    } else {
                        ret = (latex ? "J" : "J");
                        ret += "_" + std::to_string(value);
                    }
                }
                break;
            case LEAF_TAUK:
                ret = (latex ? "\\tau" : "τ");
                if (leaf_extra.k.is_symbolic()) {
                    ret += (latex ? "_{" : "_{") + leaf_extra.k.extract_symbol().str + (latex ? "}" : "}");
                } else {
                    auto value = leaf_extra.k.extract_value();
                    if (value != 2) {
                        ret += "_" + std::to_string(value);
                    }
                }
                break;
            case LEAF_XI:
                ret = (latex ? "\\xi" : "ξ");
                if (leaf_extra.k.is_symbolic()) {
                    ret += (latex ? "_{" : "_{") + leaf_extra.k.extract_symbol().str + (latex ? "}" : "}");
                } else {
                    auto value = leaf_extra.k.extract_value();
                    if (value != 1) {
                        ret += "_" + std::to_string(value);
                    }
                }
                break;
            case LEAF_PSI:
                ret = (latex ? "\\psi" : "ψ");
                if (leaf_extra.k.is_symbolic()) {
                    ret += (latex ? "_{" : "_{") + leaf_extra.k.extract_symbol().str + (latex ? "}" : "}");
                } else {
                    auto value = leaf_extra.k.extract_value();
                    if (value != 1) {
                        ret += "_" + std::to_string(value);
                    }
                }
                break;
            case LEAF_NU:
                ret = (latex ? "\\nu" : "ν");
                if (leaf_extra.k.is_symbolic()) {
                    ret += (latex ? "_{" : "_{") + leaf_extra.k.extract_symbol().str + (latex ? "}" : "}");
                } else {
                    auto value = leaf_extra.k.extract_value();
                    if (value != 1) {
                        ret += "_" + std::to_string(value);
                    }
                }
                break;
            case LEAF_RHO:
                ret = (latex ? "\\rho_{" : "ρ_{");
                if (leaf_extra.k.is_symbolic()) {
                    ret += leaf_extra.k.extract_symbol().str;
                } else {
                    auto value = leaf_extra.k.extract_value();
                    ret += std::to_string(value);
                }
                ret += ",";
                if (leaf_extra.l.is_symbolic()) {
                    ret += leaf_extra.l.extract_symbol().str;
                } else {
                    auto value = leaf_extra.l.extract_value();
                    ret += std::to_string(value);
                }
                ret += (latex ? "}" : "}");
                break;
            default:
                ret = (latex ? "?" : "?");
        }
        return ret;
    }

    std::ostream& print_inner(std::ostream& out, bool latex) const {
        switch (formula_type) {
            case FORM_ONE:
                out << "1";
                break;
            case FORM_LEAF:
                out << textify(leaf_type, latex);
                break;
            case FORM_LFUNC: {
                const Node* sub_func = sub_formula[0].get();
                if (sub_func->isOne()) {
                    out << (latex ? "\\zeta" : "ζ");
                } else {
                    out << (latex ? "L" : "L");
                }
                out << (latex ? "\\left(" : "(");
                if (!sub_func->isOne()) {
                    sub_func->print_inner(out, latex);
                    out << ", ";
                }
                if (exponent.is_symbolic()) {
                    out << exponent.extract_symbol();
                } else {
                    out << std::to_string(exponent.extract_value());
                }
                out << (latex ? "\\right)" : ")");
                break;
            }
            case FORM_POWER: {
                bool inner_is_mu = sub_formula[0].get()->isMu();
                int power_value = power.extract_value();
                if (power_value != 0) {
                   if (inner_is_mu && (power_value == 2)) {
                      out << (latex ? "\\abs{" : "|");
                   }
                   sub_formula[0].get()->print_inner(out, latex);
                   if (inner_is_mu && (power_value == 2)) {
                      out << (latex ? "}" : "|");
                   } else if (power_value != 1) {
                      out << (latex ? "^{" : "^") << std::to_string(power_value) << (latex ? "}" : "");
                   }
                }
                break;
            }
            case FORM_PRODUCT: {
                bool first = true;
                for (auto sub : sub_formula) {
                    if (!first) {
                        out << " ";
                    }
                    sub.get()->print_inner(out, latex);
                    first = false;
                }
                break;
            }
            default:
                out << "TODO";
        }
        return out;
    }

    bool isLeafOfType(LeafType requested_type) const {
        assert(isPower() || isLeaf());
        if (isPower()) {
            assert(sub_formula[0].get()->isLeaf());
            return sub_formula[0].get()->isLeafOfType(requested_type);
        }
        return leaf_type == requested_type;
    }

    int getLeafK(void) const {
        if (isPower()) {
            assert(sub_formula[0].get()->isLeaf());
            return sub_formula[0].get()->getLeafK();
        }
        return leaf_extra.k.extract_value();
    }

    int getLeafK(LeafType requested_type) const {
        assert(isLeafOfType(requested_type));
        return getLeafK();
    }

    Node() { }

public:
    std::ostream& print_inside(std::ostream& out, bool latex) const {
        print_inner(out, latex);
        return out;
    }

    std::ostream& print_full(std::ostream& out, bool latex) const {
        if (latex) {
            out << "$";
        }
        print_inside(out, latex);
        if (latex) {
            out << "$" << std::endl << std::endl;
        }
        out.flush(); /* Better for debug purposes */
        return out;
    }

    void debug_print_type(std::ostream& out) const {
        switch (formula_type) {
            case FORM_ONE:
                out << "FORM_ONE";
                break;
            case FORM_LEAF:
                out << "FORM_LEAF";
                break;
            case FORM_LFUNC:
                out << "FORM_LFUNC";
                break;
            case FORM_POWER:
                out << "FORM_POWER";
                break;
            case FORM_PRODUCT:
                out << "FORM_PRODUCT";
                break;
            default:
                out << "FORM_???";
        }
    }

    friend std::ostream& operator << (std::ostream& out, const Node &r) {
        return r.print_full(out, false);
    }

    bool isOne() const {
        return formula_type == FORM_ONE;
    }

    bool isLeaf() const {
        return formula_type == FORM_LEAF;
    }

    bool isProduct() const {
        return formula_type == FORM_PRODUCT;
    }

    bool isPower() const {
        return formula_type == FORM_POWER;
    }

    int getPower() const {
        if (isOne() || isLeaf()) {
            return 1;
        }
        assert(isPower());
        return power.extract_value();
    }

    MaybeSymbolic getPowerS() const {
        assert(isPower());
        return power;
    }

    HFormula getPowerInner() const {
        assert(isPower());
        return sub_formula[0];
    }

    size_t getProductSize() const {
        assert(isProduct());
        return sub_formula.size();
    }

    HFormula getProductElem(size_t idx) const {
        assert(isProduct() && (idx < sub_formula.size()));
        return sub_formula[idx];
    }

    bool isLFunc() const {
        return formula_type == FORM_LFUNC;
    }

    HFormula getLFuncProduct() const {
        assert(isLFunc() && sub_formula[0].get()->isProduct());
        return sub_formula[0];
    }

    MaybeSymbolic getLFuncExponentS() const {
        assert(isLFunc());
        return exponent;
    }

    int getLFuncExponent() const {
        return getLFuncExponentS().extract_value();
    }

    bool isZeta() const {
        return isLFunc() && sub_formula[0].get()->isOne();
    }

    bool isLFuncNonZeta() const {
        return isLFunc() && !isZeta();
    }

    MaybeSymbolic getZetaExponentS() const {
        assert(isZeta());
        return getLFuncExponentS();
    }

    int getZetaExponent() const {
        return getZetaExponentS().extract_value();
    }

    bool isLeafSameAs(const Node* other) const {
        assert(isLeaf());
        assert(other->isLeaf());
        return leaf_type == other->leaf_type;
    }

    bool isMu() const {
        return isLeafOfType(LEAF_MU);
    }

    bool isNu() const {
        return isLeafOfType(LEAF_NU);
    }

    bool isSigma() const {
        return isLeafOfType(LEAF_SIGMA);
    }

    bool isSigmaPrime() const {
        return isLeafOfType(LEAF_SIGMA_PRIME);
    }

    bool isTheta() const {
        return isLeafOfType(LEAF_THETA);
    }

    bool isJordanT() const {
        return isLeafOfType(LEAF_JORDAN_T);
    }

    bool isLiouville() const {
        return isLeafOfType(LEAF_LIOUVILLE);
    }

    bool isTauK() const {
        return isLeafOfType(LEAF_TAUK);
    }

    int getLeafK_dangerous() const {
        assert(isLeaf());
        return leaf_extra.k.extract_value();
    }

    int getLeafL_dangerous() const {
        assert(isLeaf());
        return leaf_extra.l.extract_value();
    }

    MaybeSymbolic getLeafKS_dangerous() const {
        assert(isLeaf());
        return leaf_extra.k;
    }

    MaybeSymbolic getLeafLS_dangerous() const {
        assert(isLeaf());
        return leaf_extra.l;
    }

    std::vector<HFormula> getSubFormula() const {
        return sub_formula;
    }
};

class NodeOne : public HFormula::Node
{
public:
    NodeOne() {
        formula_type = FORM_ONE;
    }
};

class HFormulaOne : public HFormula
{
public:
    HFormulaOne() {
        formula = std::make_shared<const Node>(NodeOne());
    }
};

class NodeLeaf : public HFormula::Node
{
public:
    NodeLeaf(LeafType type) {
        assert(type==LEAF_THETA || type==LEAF_LIOUVILLE);
        formula_type = FORM_LEAF;
        leaf_type = type;
        leaf_extra = (Node::LeafExtraArg){.k = 0, .l = 0}; /* Force init to zero, helps comparison */
    }

    NodeLeaf(LeafType type, LeafExtraArg extra) {
        assert( (   (!extra.k.is_symbolic() && (extra.k.extract_value()== 0))
                 && (!extra.l.is_symbolic() && (extra.l.extract_value()== 0)))
                || type==LEAF_MU || type==LEAF_NU || type==LEAF_SIGMA || type==LEAF_SIGMA_PRIME
                || type==LEAF_ZETAK|| type==LEAF_JORDAN_T || type==LEAF_PSI || type==LEAF_TAUK || type==LEAF_XI
                || type==LEAF_RHO
              );
        formula_type = FORM_LEAF;
        leaf_type = type;
        leaf_extra = extra;
    }
};

class HFormulaLeaf : public HFormula
{
public:
    HFormulaLeaf(Node::LeafType type) {
        formula = std::make_shared<const Node>(NodeLeaf(type));
    }

    HFormulaLeaf(Node::LeafType type, Node::LeafExtraArg extra) {
        formula = std::make_shared<const Node>(NodeLeaf(type, extra));
    }
};

class NodeLFunction : public HFormula::Node
{
public:
    NodeLFunction(const HFormula& sub_func, int sub_exponent) {
        formula_type = FORM_LFUNC;
        sub_formula.push_back(sub_func);
        exponent = sub_exponent;
    }
    NodeLFunction(const HFormula& sub_func, Symbolic sub_exponent) {
        formula_type = FORM_LFUNC;
        sub_formula.push_back(sub_func);
        exponent = sub_exponent;
    }
};

class HFormulaLFunction : public HFormula
{
public:
    HFormulaLFunction(const HFormula& sub_func, int sub_exponent) {
        formula = std::make_shared<const Node>(NodeLFunction(sub_func, sub_exponent));
    }

    HFormulaLFunction(const HFormula& sub_func, Node::Symbolic sub_exponent) {
        formula = std::make_shared<const Node>(NodeLFunction(sub_func, sub_exponent));
    }
};

class NodePower : public HFormula::Node
{
public:
    NodePower(const HFormula& part, int new_power) {
        if (part.get()->isOne() || (new_power == 0)) {
            formula_type = FORM_ONE;
        } else {
            formula_type = FORM_POWER;
            power = new_power;
            sub_formula.push_back(part);
        }
    }
};

class HFormulaPower : public HFormula
{
public:
    HFormulaPower(const HFormula& part, int power) {
        formula = std::make_shared<const Node>(NodePower(part, power));
    }
};

class NodeProduct : public HFormula::Node
{
private:
    void add_subformula(const HFormula& form) {
        if (form.get()->isOne()) {
            return;
        }
        if (form.get()->isProduct()) {
            for (auto sub : form.get()->getSubFormula()) {
                assert(!sub.get()->isProduct());
                sub_formula.push_back(sub);
            }
        } else {
            sub_formula.push_back(form);
        }
    }

public:
    NodeProduct(const HFormula& unique) {
        if (unique.get()->isOne()) {
            formula_type = FORM_ONE;
        } else {
            formula_type = FORM_PRODUCT;
            add_subformula(unique);
        }
    }

    NodeProduct(const HFormula& left, const HFormula& right) {
        if (left.get()->isOne() && right.get()->isOne()) {
            formula_type = FORM_ONE;
        } else {
            formula_type = FORM_PRODUCT;
            add_subformula(left);
            add_subformula(right);
        }
    }

    NodeProduct(const HFormula& a, const HFormula& b, const HFormula& c) {
        if (a.get()->isOne() && b.get()->isOne() && c.get()->isOne()) {
            formula_type = FORM_ONE;
        } else {
            formula_type = FORM_PRODUCT;
            add_subformula(a);
            add_subformula(b);
            add_subformula(c);
        }
    }
};

class HFormulaProduct : public HFormula
{
public:
    HFormulaProduct(const HFormula& unique) {
        formula = std::make_shared<const Node>(NodeProduct(unique));
    }

    HFormulaProduct(const HFormula& left, const HFormula& right) {
        formula = std::make_shared<const Node>(NodeProduct(left, right));
    }

    HFormulaProduct(const HFormula& a, const HFormula& b, const HFormula& c) {
        formula = std::make_shared<const Node>(NodeProduct(a, b, c));
    }
};

std::ostream& operator << (std::ostream& out, const HFormula &h) {
    return h.get()->print_full(out, false);
}

using FormulaNode = HFormula::Node;