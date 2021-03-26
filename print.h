#pragma once
#include <regex>

#if __GNUC__ > 8

#include <filesystem>
namespace fs = std::filesystem;

#else

#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;

#endif

#ifdef HAS_COLOR
#define KBLD "\x1b[1m"
#define KDIM "\x1b[2m"
#define KRED "\x1b[31m"
#define KGRN "\x1b[32m"
#define KYLW "\x1b[33m"
#define KBLU "\x1b[34m"
#define KMAG "\x1b[35m"
#define KCYN "\x1b[36m"
#define KGRY "\x1b[90m"
#define KRST "\x1b[0m"
#else
#define KBLD ""
#define KDIM ""
#define KRED ""
#define KGRN ""
#define KYLW ""
#define KBLU ""
#define KMAG ""
#define KCYN ""
#define KGRY ""
#define KRST ""
#endif

class FormulaName {
public:
    enum LeafType { LEAF_MU, LEAF_SIGMA, LEAF_THETA, LEAF_JORDAN_T, LEAF_UNKNOWN };
    typedef struct LeafExtraArg {
        int k;
        int l;
    } LeafExtraArg;

protected:
    enum FormulaType { FORM_ONE, FORM_LEAF, FORM_LFUNC, FORM_POWER, FORM_PRODUCT };
    FormulaType formula_type = FORM_ONE;
    LeafType leaf_type;
    LeafExtraArg leaf_extra;
    int power = 666;
    int exponent = 777;
    vector<const FormulaName*> sub_formula;

    string textify(LeafType in, bool latex) const {
        string ret;
        switch(in) {
            case LEAF_THETA:
                ret = (latex ? "\\theta{}" : "θ");
                break;
            case LEAF_MU:
                ret = (latex ? "\\mu{}" : "µ");
                break;
            case LEAF_SIGMA:
                ret = (latex ? "\\sigma{}" : "σ");
                if (leaf_extra.k != 1) {
                    ret += "_" + to_string(leaf_extra.k);
                }
                break;
            case LEAF_JORDAN_T:
                if (leaf_extra.k == 1) {
                    ret = (latex ? "\\phi{}" : "φ");
                } else {
                    ret = (latex ? "J" : "J");
                    ret += "_" + to_string(leaf_extra.k);
                }
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
                const FormulaName* sub_func = sub_formula[0];
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
                out << std::to_string(exponent);
                out << (latex ? "\\right)" : ")");
                break;
            }
            case FORM_POWER: {
                bool inner_is_mu = sub_formula[0]->isMu();
                if (power != 0) {
                   if (inner_is_mu && (power == 2)) {
                      out << (latex ? "\\abs{" : "|");
                   }
                   sub_formula[0]->print_inner(out, latex);
                   if (inner_is_mu && (power == 2)) {
                      out << (latex ? "}" : "|");
                   } else if (power != 1) {
                      out << (latex ? "^{" : "^") << std::to_string(power) << (latex ? "}" : "");
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
                    sub->print_inner(out, latex);
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
            assert(sub_formula[0]->isLeaf());
            return sub_formula[0]->isLeafOfType(requested_type);
        }
        return leaf_type == requested_type;
    }

    int getLeafK(LeafType requested_type) const {
        assert(isLeafOfType(requested_type));
        if (isPower()) {
            assert(sub_formula[0]->isLeaf());
            return sub_formula[0]->getLeafK(requested_type);
        }
        return leaf_extra.k;
    }

public:
    FormulaName() {
        formula_type = FORM_ONE;
    }

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
            out << "$" << endl << endl;
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

    friend std::ostream& operator << (std::ostream& out, const FormulaName &r) {
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
        return power;
    }

    size_t getProductSize() const {
        assert(isProduct());
        return sub_formula.size();
    }

    const FormulaName* getProductElem(size_t idx) const {
        assert(isProduct() && (idx < sub_formula.size()));
        return sub_formula[idx];
    }

    bool isLFunc() const {
        return formula_type == FORM_LFUNC;
    }

    const FormulaName* getLFuncProduct() const {
        assert(isLFunc() && sub_formula[0]->isProduct());
        return sub_formula[0];
    }

    int getLFuncExponent() const {
        assert(isLFunc());
        return exponent;
    }

    bool isZeta() const {
        return isLFunc() && sub_formula[0]->isOne();
    }

    int getZetaExponent() const {
        assert(isZeta());
        return getLFuncExponent();
    }

    bool isTheta() const {
        return isLeafOfType(LEAF_THETA);
    }

    bool isMu() const {
        return isLeafOfType(LEAF_MU);
    }

    bool isJordanT() const {
        return isLeafOfType(LEAF_JORDAN_T);
    }

    int getJordanTK() const {
        return getLeafK(LEAF_JORDAN_T);
    }

    bool isSigma() const {
        return isLeafOfType(LEAF_SIGMA);
    }

    int getSigmaK() const {
        return getLeafK(LEAF_SIGMA);
    }

    vector<const FormulaName*> getSubFormula() const {
        return sub_formula;
    }
};

class FormulaNameLeaf : public FormulaName
{
public:
    FormulaNameLeaf(LeafType type) {
        formula_type = FORM_LEAF;
        leaf_type = type;
    }

    FormulaNameLeaf(LeafType type, LeafExtraArg extra) {
        formula_type = FORM_LEAF;
        leaf_type = type;
        leaf_extra = extra;
    }
};

class FormulaNameLFunction : public FormulaName
{
public:
    FormulaNameLFunction(const FormulaName* sub_func, int sub_exponent) {
        formula_type = FORM_LFUNC;
        sub_formula.push_back(sub_func);
        exponent = sub_exponent;
    }
};

class FormulaNamePower : public FormulaName
{
public:
    FormulaNamePower(const FormulaName* part, int new_power) {
        if (part->isOne() || (new_power == 0)) {
            formula_type = FORM_ONE;
        } else {
            formula_type = FORM_POWER;
            power = new_power;
            sub_formula.push_back(part);
        }
    }
};

class FormulaNameProduct : public FormulaName
{
private:
    void add_subformula(const FormulaName* form) {
        if (form->isOne()) {
            return;
        }
        if (form->isProduct()) {
            for (auto sub : form->getSubFormula()) {
                assert(!sub->isProduct());
                sub_formula.push_back(sub);
            }
        } else {
            sub_formula.push_back(form);
        }
    }

public:
    FormulaNameProduct(const FormulaName* left, const FormulaName* right) {
        if (left->isOne() && right->isOne()) {
            formula_type = FORM_ONE;
        } else {
            formula_type = FORM_PRODUCT;
            add_subformula(left);
            add_subformula(right);
        }
    }
};


std::ofstream latex;

void latex_init(void) {
    fs::create_directory("tex");
    latex.open ("tex/logs.tex", std::ofstream::out);
    latex << "\\documentclass{minimal}" << endl;
    latex << "\\usepackage[utf8]{inputenc}" << endl;
    latex << "\\usepackage{amssymb}" << endl;
    latex << "\\usepackage{mathtools}" << endl;
    latex << "\\DeclarePairedDelimiter\\abs{\\lvert}{\\rvert}%" << endl;
    latex << "\\makeatletter" << endl;
    latex << "\\let\\oldabs\\abs" << endl;
    latex << "\\def\\abs{\\@ifstar{\\oldabs}{\\oldabs*}}" << endl;
    latex << "%" << endl;
    latex << "\\usepackage[svgnames]{xcolor}" << endl;
    latex << "\\begin{document}" << endl;
}

void latex_end(void) {
    latex << "\\end{document}" << endl;
    latex.close();
}
