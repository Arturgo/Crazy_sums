#include <vector>
#include <boost/multiprecision/gmp.hpp>
using namespace boost::multiprecision;
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


template<typename T>
class Mod {
public:
  Mod(int64_t _constant);
  Mod(T val);
  T value;
};

typedef Mod<BigInt> Modular;

string toString(const Modular& a) {
  return a.value.str();
}

Modular modulo(1);
vector<BigInt> inverseTable;

template<typename T>
Mod<T>::Mod(int64_t _constant) {
  value = T(_constant);
}

template<typename T>
Mod<T>::Mod(T val) {
  value = val;
}

template<typename T>
Mod<T> operator + (const Mod<T>& a, const Mod<T>& b) {
  return Mod<T>((a.value + b.value) % modulo.value);
}

template<typename T>
Mod<T> operator * (const Mod<T>& a, const Mod<T>& b) {
  return Mod<T>((modulo.value + (a.value * b.value) % modulo.value) % modulo.value);
}

template<typename T>
Mod<T> operator - (const Mod<T>& a, const Mod<T>& b) {
  return Mod<T>((modulo.value + a.value - b.value) % modulo.value);
}

template<typename T>
Mod<T> operator / (const Mod<T>&a, const Mod<T>& b) {
  return Mod<T>((a.value * inverseTable[(long)b.value]) % modulo.value);
}

template<typename T>
bool operator == (const Mod<T>& a, const Mod<T>& b) {
  return a.value == b.value;
}

Modular gcd(const Modular& a, const Modular& b) {
  return Modular(gcd(a.value, b.value));
}


Modular inverse (const Modular a) {
  return Modular(inverseTable[(int)a.value]);
}