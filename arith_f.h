#include <iostream>
#include "polynomial.h"
using namespace std;

Univariate X, U;

class RecO1 {
	public:
	size_t first_id;
	Fraction<Univariate> first_term, ratio;
	
	RecO1() {}
	
	RecO1(size_t _first_id, Fraction<Univariate> _first_term, Fraction<Univariate> _ratio) {
		first_id = _first_id;
		first_term = _first_term;
		ratio = _ratio;
	}
	
	Fraction<Univariate> get_fraction() {
		return Fraction<Univariate>(first_term) / 
			(Fraction<Univariate>(Univariate(1)) - ratio);
	}
};

RecO1 mult(const RecO1 &a, const RecO1 &b) {
	if(a.first_id > b.first_id)
		return mult(b, a);
	
	RecO1 res;
	
	res.first_id = b.first_id;
	res.ratio = a.ratio * b.ratio;
	
	Fraction<Univariate> term = a.first_term;
	
	for(size_t id = a.first_id;id < b.first_id;id++) {
		term = term * a.ratio;
	}
	
	res.first_term = term * b.first_term;
	
	return res;
}

class FArith {
public:	
	Fraction<Univariate> cst;
	vector<RecO1> recurrences;
	
	FArith() {
		cst = Fraction<Univariate>(1);
	}
	
	Fraction<Univariate> get_fraction() {
		Fraction<Univariate> sum = cst;
		
		for(RecO1 rec : recurrences) {
			sum = sum + rec.get_fraction();
		}
		
		return sum;
	}
};

FArith operator * (const FArith &a, const FArith &b) {
	FArith res;
	
	for(RecO1 recA : a.recurrences) {
		for(RecO1 recB : b.recurrences) {
			RecO1 p = mult(recA, recB);
			res.recurrences.push_back(p);
		}
	}
	
	return res;
}

/* Génération de quelques fonctions classiques */

FArith id() {
	FArith res;
	res.recurrences.push_back(RecO1(1, 
   	X, 
   	X
   ));
   return res;
}

FArith inv_id() {
	FArith res;
	res.recurrences.push_back(RecO1(1, 
   	Fraction<Univariate>(U, X),
   	Fraction<Univariate>(U, X)
   ));
   return res;
}

FArith phi() {
	FArith res;
	res.recurrences.push_back(RecO1(1,
   	X - U,
   	X
   ));
   return res;
}

FArith sigma_k(size_t k) {
	FArith res;
	res.recurrences.push_back(RecO1(1,
   	Fraction<Univariate>(-U, (U << k) - U),
   	U
   ));
   res.recurrences.push_back(RecO1(1,
   	Fraction<Univariate>((U << (2 * k))	, (U << k) - U),
   	U << k
   ));
   return res;
}

FArith one() {
	FArith res;
	res.recurrences.push_back(RecO1(1,
		U,
		U
	));
	return res;
}

FArith pow(const FArith &a, size_t exp) {
	FArith res = one();
	for(size_t i = 0;i < exp;i++) {
		res = res * a;
	}
	return res;
}

