#include "polynomial.h"
#include "matrix.h"

ModularP berlekamp(ModularP poly) {
	ModularP derivedP = derive(poly);

	if(derivedP.size() == 0) {
		cerr << "FATAL ERROR : a polynomial has an empty derivative in berlekamp." << endl;
		exit(-1);
	}
	
	ModularP pgcd = gcd(poly, derivedP);
	
	if(pgcd.size() > 1) {
		return pgcd;
	} 

	Matrix<Mod> matrice(poly.size() - 1, poly.size() - 1, 0);
	
	ModularP x_pow_ip({1});
	
	for (size_t i = 0;i + 1 < poly.size();i++) {
		if(i != 0)
			x_pow_ip = x_pow_ip << modulo;
		x_pow_ip = x_pow_ip % poly;
		
		for(size_t j = 0;j < x_pow_ip.size();j++) {
			matrice.coeffs[i][j] = x_pow_ip.getCoeff(j);
	 	}
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
	ModularP g(basis.coeffs[0]);
	
	while(g.size() <= 1) {
		index++;
		g = ModularP(basis.coeffs[index]);
	}
	
	for (int i = 0;i < modulo;i++) {
		g.setCoeff(0, i);
		
		ModularP pgcd = gcd(poly, g);
		
		if(pgcd.size() > 1 && pgcd.size() != poly.size()) {
			return pgcd;
		}
	}
	
	cerr << "FATAL ERROR : unexpected behavior in berlekamp (2)." << endl;
	exit(-1);
}

vector<ModularP> decompose(ModularP poly) {
 	vector<ModularP> decomp;
 	
 	vector<ModularP> factors = {poly};

	while(!factors.empty()) {
		ModularP factor = factors.back();
		factors.pop_back();
		
		ModularP div = berlekamp(factor);
		
		if(div == factor) {
			decomp.push_back(factor);
			continue;
		}
		
		factors.push_back(div);
		factors.push_back(factor / div);
	}
	
	return decomp;
}
