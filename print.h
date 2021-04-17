#pragma once
#include <fstream>
#if defined(__GNUC__) && !defined(__llvm__) && (__GNUC__ < 8)
    #include <experimental/filesystem>
    namespace fs = std::experimental::filesystem;
#else
    #include <filesystem>
    namespace fs = std::filesystem;
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

class Latex {
public:
    std::ofstream stream;

    Latex() {
        fs::create_directory("tex");
        stream.open ("tex/logs.tex", std::ofstream::out);
        stream << "\\documentclass{minimal}" << std::endl;
        stream << "\\usepackage[utf8]{inputenc}" << std::endl;
        stream << "\\usepackage{amssymb}" << std::endl;
        stream << "\\usepackage{mathtools}" << std::endl;
        stream << "\\DeclarePairedDelimiter\\abs{\\lvert}{\\rvert}%" << std::endl;
        stream << "\\makeatletter" << std::endl;
        stream << "\\let\\oldabs\\abs" << std::endl;
        stream << "\\def\\abs{\\@ifstar{\\oldabs}{\\oldabs*}}" << std::endl;
        stream << "%" << std::endl;
        stream << "\\usepackage[svgnames]{xcolor}" << std::endl;
        stream << "\\begin{document}" << std::endl;
    }

    ~Latex() {
        stream << "\\end{document}" << std::endl;
        stream.close();
    }
};
