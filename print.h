#pragma once
#include <filesystem>
#include <regex>
namespace fs = std::filesystem;

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
protected:
    enum FormulaType { FORM_ONE, FORM_LEAF, FORM_LFUNC, FORM_POWER, FORM_PRODUCT };
    FormulaType formula_type = FORM_ONE;
    std::string part_name;
    int power;
    vector<const FormulaName*> sub_formula;

    string textify(string const &in, bool latex) const {
        if (latex) {
            return in;
        }
        string ret = in;
        ret = std::regex_replace(ret, std::regex("\\\\sigma\\{\\}"), "σ");
        ret = std::regex_replace(ret, std::regex("\\\\phi\\{\\}"), "φ");
        return ret;
    }

    std::ostream& print_inner(std::ostream& out, bool latex) const {
        switch (formula_type) {
            case FORM_ONE:
                out << "1";
                break;
            case FORM_LEAF:
                if (part_name == "") {
                    out << "BADID";
                } else {
                    out << textify(part_name, latex);
                }
                break;
            case FORM_LFUNC: {
                const FormulaName* sub_func = sub_formula[0];
                const FormulaName* exponent = sub_formula[1];
                if (sub_func->isOne()) {
                    out << (latex ? "\\zeta" : "ζ");
                } else {
                    out << (latex ? "L" : "L");
                }
                out << textify(part_name, latex);
                out << (latex ? "\\left(" : "(");
                if (!sub_func->isOne()) {
                    sub_func->print_inner(out, latex);
                    out << ", ";
                }
                exponent->print_inner(out, latex);
                out << (latex ? "\\right)" : ")");
                break;
            }
            case FORM_POWER:
                if (power != 0) {
                   sub_formula[0]->print_inner(out, latex);
                   if (power != 1) {
                      out << (latex ? "^{" : "^") << std::to_string(power) << (latex ? "}" : "");
                   }
                }
                break;
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

    friend std::ostream& operator << (std::ostream& out, const FormulaName &r) {
        return r.print_full(out, false);
    }

    bool isOne() const {
        return formula_type == FORM_ONE;
    }

    bool isProduct() const {
        return formula_type == FORM_PRODUCT;
    }

    bool isLFunc() const {
        return formula_type == FORM_LFUNC;
    }

    bool isZeta() const {
        return isLFunc() && sub_formula[0]->isOne();
    }

    vector<const FormulaName*> getSubFormula() const {
        return sub_formula;
    }
};

class FormulaNameLeaf : public FormulaName
{
public:
    FormulaNameLeaf(std::string leaf_name) {
        formula_type = FORM_LEAF;
        part_name = leaf_name;
    }

    FormulaNameLeaf(int leaf_name) {
        formula_type = FORM_LEAF;
        part_name = std::to_string(leaf_name);
    }
};

class FormulaNameLFunction : public FormulaName
{
public:
    FormulaNameLFunction(const FormulaName* sub_func, const FormulaName* exponent) {
        formula_type = FORM_LFUNC;
        sub_formula.push_back(sub_func);
        sub_formula.push_back(exponent);
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
            if (!left->isOne()) {
                sub_formula.push_back(left);
            }
            if (!right->isOne()) {
                sub_formula.push_back(right);
            }
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
    latex << "\\begin{document}" << endl;
}

void latex_end(void) {
    latex << "\\end{document}" << endl;
    latex.close();
}
