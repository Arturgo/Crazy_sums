#pragma once
#include <cassert>
#include <iostream>
#include "polynomial.h"
#include "matrix.h"
using namespace std;

Univariate X, U, Z;
Fraction<Univariate> x, u, z;

struct FArith {
	Matrix<Fraction<Univariate>> A;
	Matrix<Fraction<Univariate>> u;
	
	Fraction<Univariate> get_fraction() {
		return (inverse(identity<Fraction<Univariate>>(A.nbRows()) - A) * u).coeffs[0].getCoeff(0);
	}
};

FArith operator * (const FArith &a, const FArith &b) {
    return {tensor(a.A, b.A), tensor(a.u, b.u)};
}

FArith operator ^ (const FArith &a, const FArith &b) {
	auto tensAId = tensor(a.A, identity<Fraction<Univariate>>(b.A.nbRows()));
	auto tensIdB = tensor(identity<Fraction<Univariate>>(a.A.nbRows()), b.A);
	Matrix<Fraction<Univariate>> cross_mat = tensAId - tensIdB;
	
	Matrix<Fraction<Univariate>> v = tensor(a.u, b.u);
	v = inverse(cross_mat) * v;

	auto va = tensAId * v;
	auto vb = -u * tensIdB * v;
	
	auto matrix = magic_op(
		tensAId, 
		tensIdB
	);
	
	auto vfinal = va;
	
	for(auto& row : vb.coeffs) {
		vfinal.coeffs.push_back(row);
	}
	
	vfinal.coeffs[0].setCoeff(0, vfinal.coeffs[0].getCoeff(0) + vb.coeffs[0].getCoeff(0));
	
   return {matrix, vfinal};
}

/* Génération de quelques fonctions classiques */

FArith id() {
	return {Matrix<Fraction<Univariate>>({
			{x}
		}),
		Matrix<Fraction<Univariate>>({
			{u}
		})
	};
}

FArith inv_id() {
	return {Matrix<Fraction<Univariate>>({
			{u / x}
		}),
		Matrix<Fraction<Univariate>>({
			{u}
		})
	};
}

FArith mobius() {
    return {Matrix<Fraction<Univariate>>({
			{z, u},
			{z, z}
		}),
		Matrix<Fraction<Univariate>>({
			{u},
			{-u}
		})
	};
}

FArith one() {
    return {Matrix<Fraction<Univariate>>({
			{u}
		}),
		Matrix<Fraction<Univariate>>({
			{u}
		})
	};
}

FArith pow(const FArith &a, size_t exp) {
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

