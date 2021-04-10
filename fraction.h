#pragma once
#include "bigint.h"

template<typename T>
class Fraction {
public:
  Fraction(const T& _numerator, const T& _denominator, bool simplify=true);
  Fraction(const T& _numerator);
  Fraction(int64_t _constant = 0);
  T getNumerator() const;
  T getDenominator() const;
  void operator += (const Fraction<T>& a);
  Fraction<T> operator + (const Fraction<T>& a) const;
  template<class U>
  friend std::ostream& operator << (std::ostream& out,const Fraction<U>& a);
private:
  T numerator, denominator;
};

template<typename T>
Fraction<T>::Fraction(const T& _numerator, const T& _denominator, bool simplify) {
  numerator = _numerator;
  denominator = _denominator;
  
  if (simplify) {
    T factor = normalFactor(numerator, denominator);
    if (normalFactorCanReduce(factor)) {
      numerator = numerator / factor;
      denominator = denominator / factor;
    }
  }
}

template<typename T>
Fraction<T>::Fraction(const T& _numerator) {
  numerator = _numerator;
  denominator = T(1);
}

template<typename T>
Fraction<T>::Fraction(int64_t _constant) {
  numerator = T(_constant);
  denominator = T(1);
}

template<typename T>
T Fraction<T>::getNumerator() const {
  return numerator;
}

template<typename T>
T Fraction<T>::getDenominator() const {
  return denominator;
}

template<typename T>
void Fraction<T>::operator += (const Fraction<T>& a) {
  if (0) {
    /* With pre-simplification */
    T factor = gcd(denominator, a.denominator);
    T x = denominator;
    T y = a.denominator;
    if (normalFactorCanReduce(factor)) {
      x = x / factor;
      y = y / factor;
    }
    numerator = numerator * y + a.numerator * x;
    denominator = x * y;
    assert(!normalFactorCanReduce(gcd(numerator, denominator)));
  } else {
    /* With post-simplification */
    numerator = numerator * a.denominator + denominator * a.numerator;
    denominator = denominator * a.denominator;

    T factor = normalFactor(numerator, denominator);
    if (normalFactorCanReduce(factor)) {
      numerator = numerator / factor;
      denominator = denominator / factor;
    }
  }
}

template<typename T>
Fraction<T> Fraction<T>::operator + (const Fraction<T>& a) const {
  Fraction<T> res = *this;
  res += a;
  return res;
}

template<class T>
std::ostream& operator << (std::ostream& out,const Fraction<T>& a) {
  out << toString(a.getNumerator(), "x") << KGRN "/" KRST
      << toString(a.getDenominator(), "x");
  return out;
}

template<typename T>
bool operator == (const Fraction<T>& a, const Fraction<T>& b) {
  return (a.getNumerator() == b.getNumerator()) && (a.getDenominator() == b.getDenominator());
}

template<typename T>
bool operator != (const Fraction<T>& a, const Fraction<T>& b) {
  return !(a == b);
}

template<typename T>
Fraction<T> operator - (const Fraction<T>& a) {
  return Fraction<T>(-a.getNumerator(), a.getDenominator(), false);
}

template<typename T>
Fraction<T> operator - (const Fraction<T>& a, const Fraction<T>& b) {
  return a + (-b);
}

template<typename T>
Fraction<T> operator * (const Fraction<T>& a, const Fraction<T>& b) {
  return Fraction<T>(
    a.getNumerator() * b.getNumerator(),
    a.getDenominator() * b.getDenominator()
  );
}

template<typename T>
Fraction<T> inverse(const Fraction<T>& a) {
  return Fraction<T>(a.getDenominator(), a.getNumerator(), false);
}

template<typename T>
Fraction<T> operator / (const Fraction<T>& a, const Fraction<T>& b) {
  return a * inverse(b);
}

template<typename T>
string toString(const Fraction<T>& a) {
  if(a.getNumerator() == T(0))
    return "0";
  if(a.getDenominator() == T(1))
    return toString(a.getNumerator());
  return toString(a.getNumerator()) + "/" + toString(a.getDenominator());
}

template<typename T>
bool is_positive(const Fraction<T>& a) {
  assert(a.getDenominator() > T(0));
  return a.getNumerator() > T(0);
}

template<typename T>
Fraction<T> abs(const Fraction<T>& a) {
  return is_positive(a) ? a : -a;
}

template<typename T>
bool operator < (const Fraction<T>& a, const Fraction<T>& b) {
  return is_positive(b - a);
}

template<typename T>
bool operator > (const Fraction<T>& a, const Fraction<T>& b) {
  return b < a;
}

template<typename T>
bool is_integer(const Fraction<T>& a) {
  return a.getDenominator() == T(1);
}

typedef Fraction<SomeInt> Rational;

ostream& operator << (ostream& out, const Rational &r) {
   out << toString(r);
   return out;
}

Rational round(Rational a) {
   return Rational((a.getNumerator() + a.getDenominator() / 2) / a.getDenominator());
}
