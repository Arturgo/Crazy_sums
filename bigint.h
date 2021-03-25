#pragma once
#include <boost/multiprecision/gmp.hpp>
#include <string>
#include <vector>
using namespace std;

class SmallInt {
    typedef int_fast64_t smallint;
protected:
    smallint n;

public:
    SmallInt() {
        n = 0;
    }

    SmallInt(smallint val) {
        n = val;
    }

    string str() const {
        return std::to_string(n);
    }

    SmallInt operator ++ (void) {
        return n++;
    }

    SmallInt operator ++ (int) {
        return n++;
    }

    SmallInt operator -- (void) {
        return n--;
    }

    SmallInt operator -- (int) {
        return n--;
    }

    SmallInt operator + (const SmallInt& other) const {
        return n + other.n;
    }

    SmallInt operator - (const SmallInt& other) const {
        return n - other.n;
    }

    SmallInt operator - (void) const {
        return -n;
    }

    SmallInt operator * (const SmallInt& other) const {
        return n * other.n;
    }

    SmallInt operator / (const SmallInt& other) const {
        return n / other.n;
    }

    SmallInt operator % (const SmallInt& other) const {
        return n % other.n;
    }

    bool operator == (const SmallInt& other) const {
        return n == other.n;
    }

    void operator += (const SmallInt& other) {
        n += other.n;
    }

    void operator -= (const SmallInt& other) {
        n -= other.n;
    }

    bool operator != (const SmallInt& other) const {
        return n != other.n;
    }

    bool operator < (const SmallInt& other) const {
        return n < other.n;
    }

    bool operator > (const SmallInt& other) const {
        return n > other.n;
    }

    bool operator <= (const SmallInt& other) const {
        return n <= other.n;
    }

    bool operator >= (const SmallInt& other) const {
        return n >= other.n;
    }
};

#if 0
using SomeInt = boost::multiprecision::mpz_int;
#else
using SomeInt = SmallInt;
#endif

SomeInt gcd(const SomeInt& a, const SomeInt& b) {
  if(b == SomeInt(0))
    return max(a, -a);
  return gcd(b, a % b);
}

SomeInt normalFactor(const SomeInt& a, const SomeInt& b) {
  if(b >= SomeInt(0))
    return gcd(a, b);
  else
    return -gcd(a, b);;
}

string toString(const SomeInt& a) {
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
