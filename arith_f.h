#pragma once
#include <cassert>
#include "matrix.h"
#include "polynomial.h"
using namespace std;

Univariate X, U, Z;
Fraction<Univariate> x, u, z;

std::chrono::duration<float> v1, v2;

typedef Matrix<Fraction<Univariate>> FArithMatrix;

struct FArith {
    FArithMatrix A;
    FArithMatrix u;

    Fraction<Univariate> get_fraction() const {
        auto id = identity<Fraction<Univariate>>(A.nbRows());
        auto t0 = std::chrono::high_resolution_clock::now();
        auto mat = inverse(id - A);
        auto t1 = std::chrono::high_resolution_clock::now();
        auto res = single_product_element(mat, u, 0, 0);
        auto t2 = std::chrono::high_resolution_clock::now();

        v1 += t1 - t0;
        v2 += t2 - t1;

        return res;
    }

    void remove_row(size_t row) {
        FArithMatrix nA(A.coeffs.size() - 1, A.coeffs.size() - 1);
        FArithMatrix nu(A.coeffs.size() - 1, 1);

        size_t nRow = 0;
        for(size_t iRow = 0;iRow < A.coeffs.size();iRow++) {
            if(iRow == row)
                continue;

            nu.coeffs[nRow].setCoeff(0, u.coeffs[iRow].getCoeff(0));

            size_t nCol = 0;
            for(size_t iCol = 0;iCol < A.coeffs.size();iCol++) {
                if(iCol == row)
                    continue;

                nA.coeffs[nRow].setCoeff(nCol, A.coeffs[iRow].getCoeff(iCol));

                nCol++;
            }
            nRow++;
        }

        A = nA;
        u = nu;
    }

    void simplify() {
        for(size_t iCol = 1;iCol < A.nbCols();iCol++) {
            bool is_col_zero = true;

            for(size_t iRow = 0;is_col_zero && iRow < A.nbRows();iRow++) {
                if(iRow != iCol && !is_zero(A.coeffs[iRow].getCoeff(iCol))) {
                    is_col_zero = false;
                }
            }

            if(is_col_zero) {
                remove_row(iCol);
                simplify();
                return;
            }
        }

        FArithMatrix new_mat = transpose(A);
        new_mat.coeffs.push_back(transpose(u).coeffs[0]);
        new_mat = transpose(new_mat);

        FArithMatrix basis = kernel_basis(new_mat);

        if(basis.nbRows() == 0) {
            return;
        }

        size_t row = basis.coeffs[0].coeffs[0].first;
        if(row == 0)
            row = basis.coeffs[0].coeffs[1].first;

        Fraction<Univariate> bCoeff = basis.coeffs[0].getCoeff(A.nbCols());
        basis.coeffs[0].setCoeff(A.nbCols(), 0);

        for(size_t iRow = 0;iRow < A.nbRows();iRow++) {
            if(iRow != row) {
                u.coeffs[iRow].setCoeff(0, u.coeffs[iRow].getCoeff(0) - A.coeffs[iRow].getCoeff(row) / basis.coeffs[0].getCoeff(row) * bCoeff);

                A.coeffs[iRow] = A.coeffs[iRow] - A.coeffs[iRow].getCoeff(row) / basis.coeffs[0].getCoeff(row) * basis.coeffs[0];
            }
        }

        simplify();
    }

    friend std::ostream& operator << (std::ostream& out, const FArith &a) {
        out << "FArith:" << endl;
        out << "A =" << endl;
        out << a.A << endl;
        out << "u =" << endl;
        out << a.u << endl;
        return out;
    }
};

FArith operator * (const FArith &a, const FArith &b) {
    FArith res = {.A = tensor(a.A, b.A), .u = tensor(a.u, b.u)};
    res.simplify();
    return res;
}

FArith operator ^ (const FArith &a, const FArith &b) {
    FArithMatrix tensAId = tensor(a.A, identity<Fraction<Univariate>>(b.A.nbRows()));
    FArithMatrix tensIdB = tensor(identity<Fraction<Univariate>>(a.A.nbRows()), b.A);
    FArithMatrix cross_mat = tensAId - tensIdB;

    FArithMatrix v = tensor(a.u, b.u);
    v = inverse(cross_mat) * v;

    FArithMatrix va = tensAId * v;
    FArithMatrix vb = -u * tensIdB * v;

    FArithMatrix matrix = magic_op(tensAId, tensIdB);

    FArithMatrix vfinal = va;

    for(auto& row : vb.coeffs) {
        vfinal.coeffs.push_back(row);
    }

    vfinal.coeffs[0].setCoeff(0, vfinal.coeffs[0].getCoeff(0) + vb.coeffs[0].getCoeff(0));

    FArith res = {.A = matrix, .u = vfinal};
    res.simplify();
    return res;
}

FArith one() {
    return {
        .A = FArithMatrix({
            {u}
        }),
        .u = FArithMatrix({
            {u}
        })
    };
}

FArith pow(const FArith &a, size_t exp) {
    /* Exponentiation by squaring */
    if (exp <= 0) {
        return one();
    }
    FArith res = pow(a, exp/2);
    res = res * res;
    if (exp%2 != 0) {
        res = res * a;
    }
    return res;
}

/* Generate some usual function */

FArith id() {
    return {
        .A = FArithMatrix({
            {x}
        }),
        .u = FArithMatrix({
            {u}
        })
    };
}

FArith inv_id() {
    return {
        .A = FArithMatrix({
            {u / x}
        }),
        .u = FArithMatrix({
            {u}
        })
    };
}

FArith mobius_k(size_t k) {
    FArith res = {
        .A = FArithMatrix(k + 1, k + 1),
        .u = FArithMatrix(k + 1, 1)
    };

    for(size_t i = 0;i < k;i++) {
        res.A.coeffs[i].setCoeff(i + 1, u);
    }
    res.u.coeffs[k].setCoeff(0, -u);
    res.u.coeffs[0].setCoeff(0, u);
    return res;
}

FArith mobius() {
    return mobius_k(1);
}

FArith nu_k(size_t k) {
    FArith res = {
        .A = FArithMatrix(k, k),
        .u = FArithMatrix(k, 1)
    };

    for(size_t i = 0;i < k - 1;i++) {
        res.A.coeffs[i].setCoeff(i + 1, u);
    }
    res.A.coeffs[k - 1].setCoeff(0, u);

    res.u.coeffs[0].setCoeff(0, u);
    return res;
}

FArith jordan_totient(size_t k) {
    assert(k > 0);
    return pow(id(), k) ^ mobius();
}

FArith phi() {
    return jordan_totient(1);
}

FArith sigma_k(size_t k) {
    assert(k > 0);
    return one() ^ pow(id(), k);
}

FArith sigma_prime_k(size_t k) {
    return {
        .A = FArithMatrix({
            {z, Fraction<Univariate>(U, U + (U << k)), Fraction<Univariate>(U, U + (U << k)) },
            {z, -(U << k), z},
            {z, z, u}
        }),
        .u = FArithMatrix({
            {u},
            {- (U << (2 * k))},
            {u}
        })
    };
}

FArith theta() {
    return pow(mobius(), 2) ^ one();
}

FArith liouville() {
    return {
        .A = FArithMatrix({
            {-u}
        }),
        .u = FArithMatrix({
            {u}
        })
    };
}

FArith zeta_1() {
    return {
        .A = FArithMatrix({
            {x}
        }),
        .u = FArithMatrix({
            {u}
        })
    };
}

FArith simple_factor(size_t k) {
    return {
        .A = FArithMatrix({
            {u, u},
            {z, u},
        }),
        .u = FArithMatrix({
            {u},
            {Fraction<Univariate>(1, k)},
        })
    };
}

FArith tau(size_t k) {
    FArith res = one();
    for(size_t i = 1;i < k;i++) {
        res = res * simple_factor(i);
    }
    return res;
}

FArith nb_divisors() {
    return tau(2);
}

FArith xi_k(size_t k) {
    FArith res = {
        .A = FArithMatrix(k, k),
        .u = FArithMatrix(k, 1)
    };

    for(size_t i = 0;i < k - 1;i++) {
        res.A.coeffs[i].setCoeff(i + 1, u);
        res.u.coeffs[i].setCoeff(0, u);
    }
    res.u.coeffs[k - 1].setCoeff(0, u);

    return res;
}

FArith rho_k_s(size_t k, size_t s) {
    return pow(id(), k) ^ nu_k(s);
}

FArith psi_k(size_t k) {
    return pow(id(), k) ^ (mobius() * mobius());
}

/********************************************************
 *                Useless functions                     *
 ********************************************************/

/* We use directly liouville && sigma_prime_k */
FArith beta_k(size_t k) {
    return sigma_prime_k(k) * liouville();
}

/* We use directly nu_k */
FArith precompose_with_kth_power(FArith f, size_t k) {
    return f * nu_k(k);
}
