#pragma once
#include <boost/multiprecision/gmp.hpp>
#include <string>
#include <vector>
using namespace boost::multiprecision;
using namespace std;
using BigInt = mpz_int;

BigInt gcd(const BigInt& a, const BigInt& b) {
  if(b == 0)
    return max(a, -a);
  return gcd(b, a % b);
}

BigInt normalFactor(const BigInt& a, const BigInt& b) {
  if(b >= 0)
    return gcd(a, b);
  else
    return -gcd(a, b);;
}

string toString(const BigInt& a) {
  return a.str();
}



class Mod {
public:
  Mod(int val = 0);
  int value;
};

int modulo;
vector<Mod> inverseTable;

void precomputeInverses(int PRIME) {
	modulo = PRIME;
	
   inverseTable.push_back(Mod(0));
   for(int i = 1;i < PRIME;i++) {
   	int j = 1;
      for(;(i * j) % PRIME != 1;j++);
      inverseTable.push_back(Mod(j));
   }
}

string toString(const Mod& a) {
	if(a.value > modulo / 2)
		return to_string(a.value - modulo);
	return to_string(a.value);
}

Mod::Mod(int _constant) {
	value = ((_constant % modulo) + modulo) % modulo;
}

Mod operator + (const Mod& a, const Mod& b) {
	return Mod(a.value + b.value);
}

Mod operator * (const Mod& a, const Mod& b) {
	return Mod(a.value * b.value);
}

Mod operator - (const Mod& a, const Mod& b) {
	return Mod(a.value - b.value);
}

Mod inverse(const Mod& a) {
	return inverseTable[a.value];
}

Mod operator / (const Mod&a, const Mod& b) {
	return a * inverse(b);
}

bool operator == (const Mod& a, const Mod& b) {
	return a.value == b.value;
}
ostream& operator << (ostream& out, const Mod &r) {
   out << toString(r);
   return out;
}

/*
NON-SENS:
Modular gcd(const Modular& a, const Modular& b) {
  return Modular(gcd(a.value, b.value));
}*/
