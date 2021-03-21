#pragma once
#include <cassert>
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
            (Fraction<Univariate>(U) - ratio);
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

vector<RecO1> conv(const RecO1 &a, const RecO1 &b) {
	vector<RecO1> recs;

	RecO1 recA, recB;
	
	recA.first_id = a.first_id + b.first_id;
	recA.first_term = a.first_term * b.first_term * a.ratio / (a.ratio - b.ratio);
	recA.ratio = a.ratio;
	
	recB.first_id = a.first_id + b.first_id;
	recB.first_term = - a.first_term * b.first_term * b.ratio / (a.ratio - b.ratio);
	recB.ratio = b.ratio;
	
	recs.push_back(recA);
	recs.push_back(recB);
	
	return recs;
}

class FArith {
public:
	vector<Fraction<Univariate>> csts;
    vector<RecO1> recurrences;
    
    FArith(Fraction<Univariate> cst = Fraction<Univariate>(0)) {
    	csts = {cst};
    }
    
    void reduce() {
    	while(!csts.empty() && csts.back() == Fraction<Univariate>(0))
			csts.pop_back();
    }
    
    void ajoute(RecO1 rec) {
    	if(rec.ratio == Fraction<Univariate>(0)) {
    		csts.resize(max(rec.first_id + 1, csts.size()), 0);
    		csts[rec.first_id] = csts[rec.first_id] + rec.first_term;
    		return;
    	}
    
    	for(size_t iRec = 0;iRec < recurrences.size();iRec++) {
    		RecO1 r = recurrences[iRec];
    		if(r.ratio == rec.ratio && r.first_id == rec.first_id) {
    			csts.resize(max(max(r.first_id, rec.first_id) + 1, csts.size()));
    			
    			while(r.first_id < rec.first_id) {
					csts[r.first_id] = csts[r.first_id] + r.first_term;
    				r.first_term = r.first_term * r.ratio;
    				r.first_id++;
    			}
    			
    			while(rec.first_id < r.first_id) {
    				csts[rec.first_id] = csts[rec.first_id] + rec.first_term;
    				rec.first_term = rec.first_term * rec.ratio;
    				rec.first_id++;
    			}
    			
    			r.first_term = r.first_term + rec.first_term;
    			
    			recurrences[iRec] = r;
    			
    			return;
			}
		}
		
		recurrences.push_back(rec);
    }

    Fraction<Univariate> get_fraction() {
        Fraction<Univariate> sum;
        
        for(auto cst : csts) {
        	sum = sum + cst;
        }

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
            res.ajoute(p);
        }
    }
    
    res.csts.resize(max(res.csts.size(), max(a.csts.size(), b.csts.size())), 0);
    
    for(RecO1 recB : b.recurrences) {
    	Fraction<Univariate> value = recB.first_term;
    	for(size_t iCst = recB.first_id;iCst < a.csts.size();iCst++) {
			res.csts[iCst] = res.csts[iCst] + a.csts[iCst] * value;
			value = value * recB.ratio;
		}
	}
	
	for(RecO1 recA : a.recurrences) {
    	Fraction<Univariate> value = recA.first_term;
    	for(size_t iCst = recA.first_id;iCst < b.csts.size();iCst++) {
			res.csts[iCst] = res.csts[iCst] + b.csts[iCst] * value;
			value = value * recA.ratio;
		}
	}
	
	for(size_t iCst = 0;iCst < min(a.csts.size(), b.csts.size());iCst++) {
		res.csts[iCst] = res.csts[iCst] + a.csts[iCst] * b.csts[iCst];
	}

	res.reduce();
    return res;
}

FArith operator ^ (const FArith &a, const FArith &b) {	
    FArith res;

    for(RecO1 recA : a.recurrences) {
        for(RecO1 recB : b.recurrences) {
            vector<RecO1> c = conv(recA, recB);
            for(RecO1 p : c) {
           		res.ajoute(p);
           	}
        }
    }
    
    for(size_t iCst = 0;iCst < a.csts.size();iCst++) {
		for(RecO1 recB : b.recurrences) {
			recB.first_id = recB.first_id + iCst;
			recB.first_term = recB.first_term * a.csts[iCst];
			res.ajoute(recB);
		}
	}
	
	for(size_t iCst = 0;iCst < b.csts.size();iCst++) {
		for(RecO1 recA : a.recurrences) {
			recA.first_id = recA.first_id + iCst;
			recA.first_term = recA.first_term * b.csts[iCst];
			res.ajoute(recA);
		}
	}
	
	res.csts.resize(max(res.csts.size(), a.csts.size() + b.csts.size()), 0);
	
	for(size_t iCstA = 0;iCstA < a.csts.size();iCstA++) {
		for(size_t iCstB = 0;iCstB < b.csts.size();iCstB++) {
			res.csts[iCstA + iCstB] = res.csts[iCstA + iCstB] + a.csts[iCstA] * b.csts[iCstB];
		}
	}

	res.reduce();
    return res;
}

/* Génération de quelques fonctions classiques */

FArith id() {
    FArith res(U);
    res.ajoute(RecO1(1,
        X,
        X
    ));
    return res;
}

FArith inv_id() {
    FArith res(U);
    res.ajoute(RecO1(1,
        Fraction<Univariate>(U, X),
        Fraction<Univariate>(U, X)
    ));
    return res;
}

FArith mobius() {
    FArith res(U);
    res.ajoute(RecO1(1,
        -U,
        0
    ));
    return res;
}

FArith one() {
    FArith res(U);
    res.ajoute(RecO1(1,
        U,
        U
    ));
    return res;
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
