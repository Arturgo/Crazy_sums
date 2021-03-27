#include "polynomial.h"
#include "matrix.h"

#if 0
Univariate berlekamp(Univariate poly) {
	Univariate derivedP = derive(poly);

	if(derivedP.size() == 0) {
		cerr << "FATAL ERROR : a polynomial has an empty derivative in berlekamp." << endl;
		exit(-1);
	}
	
	Univariate pgcd = gcd(poly, derivedP);
	
	if(pgcd.size() > 1) {
		return pgcd;
	} 

	Matrix<Mod> matrice(poly.size() - 1, poly.size() - 1, 0);
	
	Univariate x_pow_ip({1});
	
	for (size_t i = 0;i + 1 < poly.size();i++) {
		if(i != 0)
			x_pow_ip = x_pow_ip << modulo;
		x_pow_ip = x_pow_ip % poly;
		
		matrice.coeffs[i] = MatrixRow<Mod>(x_pow_ip.coeffs);
	}
	
	matrice = matrice - identity<Mod>(poly.size() - 1);
	Matrix<Mod> basis = kernel_basis(matrice);
	
	if(basis.nbRows() == 1) {
		return poly;
	} 
	else if(basis.nbRows() == 0) {
		cerr << "FATAL ERROR : unexpected behavior in berlekamp (1)." << endl;
		exit(-1);
	}
	
	int index = 0;
	Univariate g(basis.coeffs[0]);
	
	while(g.size() <= 1) {
		index++;
		g = Univariate(basis.coeffs[index]);
	}
	
	for (int i = 0;i < modulo;i++) {
		g.setCoeff(0, i);
		
		Univariate pgcd = gcd(poly, g);
		
		if(pgcd.size() > 1 && pgcd.size() != poly.size()) {
			return pgcd;
		}
	}
	
	cerr << "FATAL ERROR : unexpected behavior in berlekamp (2)." << endl;
	exit(-1);
}

vector<Univariate> decompose(Univariate poly) {
 	vector<Univariate> decomp;
 	
 	vector<Univariate> factors = {poly};

	while(!factors.empty()) {
		Univariate factor = factors.back();
		factors.pop_back();
		
		Univariate div = berlekamp(factor);
		
		if(div == factor) {
			decomp.push_back(factor);
			continue;
		}
		
		factors.push_back(div);
		factors.push_back(factor / div);
	}
	
	return decomp;
}
#endif
