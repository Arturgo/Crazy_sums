#pragma once
#include <cassert>
#include "matrix.h"
#include "polynomial.h"
using namespace std;

Univariate X, U, Z;
Fraction<Univariate> x, u, z;

typedef Matrix<Fraction<Univariate>> FArithMatrix;

struct FArith {
    FArithMatrix A;
    FArithMatrix u;

    Fraction<Univariate> get_fraction() {
        return (inverse(identity<Fraction<Univariate>>(A.nbRows()) - A) * u).coeffs[0].getCoeff(0);
    }
};

FArith operator * (const FArith &a, const FArith &b) {
    return {.A = tensor(a.A, b.A), .u = tensor(a.u, b.u)};
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

    return {.A = matrix, .u = vfinal};
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

FArith mobius() {
    return {
        .A = FArithMatrix({
            {z, u},
            {z, z}
        }),
        .u = FArithMatrix({
            {u},
            {-u}
        })
    };
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

FArith theta() {
    return pow(mobius(), 2) ^ one();
}

