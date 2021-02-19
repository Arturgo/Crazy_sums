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
