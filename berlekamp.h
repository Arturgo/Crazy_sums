#include "polynomial.h"


IntegralsP remainder(IntegralsP num, IntegralsP div, BigInt p, BigInt invertTable[]) {
  
  int headD = static_cast<int>(div.getCoeff(div.size() - 1));
  BigInt invD = invertTable[headD];

  for (int iCoeff = 0; iCoeff < (int)div.size(); iCoeff++) {
    div.setCoeff(iCoeff, (div.getCoeff(iCoeff) * invD) % p);
  }
  div.reduce();
  for (int iCoeff = 0; iCoeff < (int)num.size(); iCoeff++) {
    num.setCoeff(iCoeff, (p + (num.getCoeff(iCoeff) * invD) % p) % p);
  }
  num.reduce();
  int c = 0;
  //cout << "entered " << toString(num, "X") << " by " << toString(div, "X") << endl;
  while (num.size() >= div.size()) {
    BigInt lead = num.getCoeff(num.size() - 1);
    num = num - (lead * div << (num.size() - div.size()));
    for (int iCoeff = 0; iCoeff < (int)num.size(); iCoeff++) {
      num.setCoeff(iCoeff, (p + (num.getCoeff(iCoeff) * invD) % p) % p);
    }
    //cout << "num' : " << toString(num, "X")  << " " << toString(div, "X")<< endl;
    assert(c++<1000);
  }
  //cout << "left" << endl;

  return num;
}

IntegralsP division(IntegralsP num, IntegralsP div, BigInt p, BigInt invertTable[]) {
  vector<BigInt>vect;
  
  int headD = static_cast<int>((p + div.getCoeff(div.size() - 1) % p ) % p);
  BigInt invD = invertTable[headD];

  for (int iCoeff = 0; iCoeff < (int)div.size(); iCoeff++) {
    div.setCoeff(iCoeff, (div.getCoeff(iCoeff) * invD) % p);
  }
  int c = 0;
  while (num.size() >= div.size()) {
    BigInt lead = num.getCoeff(num.size() - 1);
    vect.push_back(lead);
    for (int i = 0; i < (int)div.size(); i++) {
        num.setCoeff(i + num.size() - div.size(), (p + (num.getCoeff(i + num.size() - div.size()) - lead * div.getCoeff(i))%p )%p);
    }
    //num = num - (lead * div << (num.size() - div.size()));
   // cout << toString(num, "X") << endl;
    assert(c++<100);
  }

  vector<BigInt>vect2;
  if (num.size() == 0) {
    for (int i = vect.size() - 1; i >= 0; i--) {
      vect2.push_back(vect[i]);
    }
  }
  //cout << "quotient " << toString(IntegralsP(vect2), "X")  << " " << toString(num, "X") << endl; 
  return IntegralsP(vect2);
}

template<typename T>
Polynomial<T> berlekamp(Polynomial<T> polynomial) {
  Polynomial<T> derivedP = derive(polynomial);
  //cout << "computed derivation" << endl;
  derivedP.reduce();
  bool t = true;
  for (int iCoeff = 0; iCoeff < (int)derivedP.size(); iCoeff++) {
    t = t && (derivedP.getCoeff(iCoeff) == T(0));
  }
  //cout << "here" << endl;
  if (t) {
    //TODO
    cout << "case 1" << endl;
    cout << "derived : " << toString(derivedP, "X") << endl;
    cout << toString(polynomial, "X") << endl;
    assert(false);
  } else {
    Polynomial<T> gcd = modular_gcd(polynomial, derivedP);
    //cout << "computed gcd" << endl;
    if (gcd.size() > 1) {
      //cout << "case 2 " << toString(gcd, "X") << " = " << toString(polynomial, "X") << " ^ " << toString(derivedP, "X") << endl;
      return gcd;
    } else {
      //cout << "case 3 " << toString(gcd, "X") << " = " << toString(polynomial, "X") << " ^ " << toString(derivedP, "X") << endl;
      Matrix<T> matrice(polynomial.size()-1, polynomial.size()-1, T(0));
      vector<T> v;
      v.push_back(T(1));
      //cout << "l.91" << endl;
      Polynomial<T> x_pow_ip(v);
      for (int i = 0; i < (int)polynomial.size()-1; i++) {
        if (i!=0) x_pow_ip = x_pow_ip << (int)modulo.value;
        x_pow_ip = remainder(x_pow_ip, polynomial);
        //cout << toString(x_pow_ip, "X") << endl;
        for (int j = 0; j < (int)x_pow_ip.size(); j++) {
          matrice.coeffs[i][j] = x_pow_ip.getCoeff(j);
        }
      }

      //cout << "here" << endl;

      for (int i = 0; i < (int)polynomial.size()-1; i++) {
        matrice.coeffs[i][i] = matrice.coeffs[i][i] - T(1);
      }

      Matrix<T> basis = kernel_basis(matrice);

      if (basis.coeffs.size() == 1) {
        return polynomial;
      } else if (basis.coeffs.size() == 0) {
        assert(false);
      } else {
        //cout << "basis size : " << basis.coeffs.size() << endl;
        //cout << "case 3 " << toString(gcd, "X") << " = " << toString(polynomial, "X") << " ^ " << toString(derivedP, "X") << " " << basis.coeffs.size() << endl;
        int index = 0;
        Polynomial<T> g(basis.coeffs[index]);
        g.reduce();
        while (g.size() <= 1) {
            assert(++index < basis.coeffs.size());
            g = Polynomial<T>(basis.coeffs[index]);
            g.reduce();
        }
        int p2 = (int)modulo.value;
        for (int i = 0; i < p2; i++) {
            g.setCoeff(0, T(i));
            //IntegralsP g2(basis.coeffs[index]);
            //cout << "entered gcd" << endl;
            gcd = modular_gcd(polynomial, g);
            //cout << "left gcd" << endl; 
            if (gcd.size() > 1) {
                //cout << "divisor found " << toString(gcd, "X") << " : " << toString(gcd, "X") << endl;
                //cout <<  i << " -> " << toString(g, "X") << endl;
                return gcd;
            }
            //cout << "here " << toString(gcd, "X") << " " << toString(g2, "X") << " " << toString(polynomial, "X") << endl;
        }
        cout << "P " << toString(polynomial, "X") << endl;
        for (vector<T> row : basis.coeffs) {
          for(T x : row)
            cout << toString(x) << " ";
          cout << endl;
        }

        cout << "Matrix : " << endl;
        for (vector<T> row : matrice.coeffs) {
          for(T x : row)
            cout << toString(x) << " ";
          cout << endl;
        }
        assert(false);
      }
    }
  }
}

IntegralsP berlekamp(IntegralsP polynomial, BigInt& p, BigInt invertTable[]) {
  IntegralsP derivedP = derive(polynomial, p);
  //cout << "computed derivation" << endl;
  derivedP.reduce();
  bool t = true;
  for (int iCoeff = 0; iCoeff < (int)derivedP.size(); iCoeff++) {
    BigInt coeff = derivedP.getCoeff(iCoeff);
    BigInt coeffRed = coeff % p;
    t = t && (coeffRed == BigInt(0));
  }
  //cout << "here" << endl;
  if (t) {
    //TODO
    cout << "case 1" << endl;
    cout << "derived : " << toString(derivedP, "X") << endl;
    cout << toString(polynomial, "X") << endl;
    assert(false);
  } else {
    IntegralsP gcd = integrals_gcd(polynomial, derivedP, p, invertTable);
    //cout << "computed gcd" << endl;
    if (gcd.size() > 1) {
      //cout << "case 2 " << toString(gcd, "X") << " = " << toString(polynomial, "X") << " ^ " << toString(derivedP, "X") << endl;
      return gcd;
    } else {
      //cout << "case 3 " << toString(gcd, "X") << " = " << toString(polynomial, "X") << " ^ " << toString(derivedP, "X") << endl;
      Matrix<BigInt> matrice(polynomial.size()-1, polynomial.size()-1, BigInt(0));
      vector<BigInt> v;
      v.push_back(BigInt(1));
      IntegralsP x_pow_ip(v);
      for (int i = 0; i < (int)polynomial.size()-1; i++) {
        if (i!=0) x_pow_ip = x_pow_ip << (int)p;
        x_pow_ip = remainder(x_pow_ip, polynomial, p, invertTable);
        for (int j = 0; j < (int)x_pow_ip.size(); j++) {
          matrice.coeffs[i][j] = x_pow_ip.getCoeff(j);
        }
      }

      //cout << "here" << endl;

      for (int i = 0; i < (int)polynomial.size()-1; i++) {
        matrice.coeffs[i][i] = (matrice.coeffs[i][i] - BigInt(1) + p)%p;
      }
      //cout << "l.100" << endl;
      Matrix<BigInt> basis = kernel_basis_finite_set(matrice, p, invertTable);
      //cout << "l.101" << endl;
      if (basis.coeffs.size() == 1) {
        return polynomial;
      } else if (basis.coeffs.size() == 0) {
        assert(false);
      } else {
        //cout << "case 3 " << toString(gcd, "X") << " = " << toString(polynomial, "X") << " ^ " << toString(derivedP, "X") << " " << basis.coeffs.size() << endl;
        int index = -1;
        while (basis.coeffs[++index].size() <= 1) {
            assert(index < basis.coeffs.size());
        }
        IntegralsP g(basis.coeffs[index]);
        int p2 = (int)p;
        for (int i = 0; i < p2; i++) {
            g.setCoeff(0, BigInt(i));
            //IntegralsP g2(basis.coeffs[index]);
            //cout << "entered gcd" << endl;
            gcd = integrals_gcd(polynomial, g, p, invertTable);
            //cout << "left gcd" << endl; 
            if (gcd.size() > 1) {
                //cout << "divisor found " << toString(gcd, "X") << " : " << toString(gcd, "X") << endl;
                //cout <<  i << " -> " << toString(g, "X") << endl;
                return gcd;
            }
            //cout << "here " << toString(gcd, "X") << " " << toString(g2, "X") << " " << toString(polynomial, "X") << endl;
        }
        cout << "P " << toString(polynomial, "X") << endl;
        for (vector<BigInt> row : basis.coeffs) {
          for(BigInt x : row)
            cout << x << " ";
          cout << endl;
        }

        cout << "Matrix : " << endl;
        for (vector<BigInt> row : matrice.coeffs) {
          for(BigInt x : row)
            cout << x << " ";
          cout << endl;
        }
        assert(false);
      }
    }
  }
}

vector<IntegralsP> decompose(Univariate b_raw, BigInt& p, BigInt invertTable[]) {
  vector<IntegralsP> v;
  IntegralsP b = toBigInt(b_raw);
  //cout << "b = " << toString(b, "X") << endl;
  //cout << "b'= " << toString(b, "X") << endl;
  while (b.size() > 1) {
    //cout << "enter kerlekamp" << endl;
    for (int iCoeff = 0; iCoeff < (int)b.size(); iCoeff++) {
        b.setCoeff(iCoeff, (p + b.getCoeff(iCoeff)%p )%p);
    }
    //cout << "begins computing" << endl;
    IntegralsP div = berlekamp(b, p, invertTable);
    //cout << "finished computing" << endl;
    for (int iCoeff = 0; iCoeff < (int)div.size(); iCoeff++) {
        div.setCoeff(iCoeff, (p + div.getCoeff(iCoeff)%p )%p);
    }
    //cout << "berlekamp ended with : "<< toString(b, "X") << " " << toString(div, "X") << endl;
    v.push_back(div);
    bool t = true;
    while (t) {
      IntegralsP quotient = division(b, div, p, invertTable);
//      cout << toString(quotient, "X") << endl;
      if (quotient.size() == 0) {
        t = false;
      } else {
        b = quotient;
      }
    }
  }
  return v;
}