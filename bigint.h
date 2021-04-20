#pragma once
#include <string>
#include <vector>

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

    std::string str() const {
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

    bool operator != (const SmallInt& other) const {
        return n != other.n;
    }

    void operator += (const SmallInt& other) {
        n += other.n;
    }

    void operator -= (const SmallInt& other) {
        n -= other.n;
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

    int to_int() {
        return n;
    }
};

#if 0
#include <boost/multiprecision/gmp.hpp>
using SomeInt = boost::multiprecision::mpz_int;
#else
using SomeInt = SmallInt;
#endif

SomeInt gcd(const SomeInt& a, const SomeInt& b) {
  if(b == SomeInt(0))
    return std::max(a, -a);
  return gcd(b, a % b);
}

SomeInt normalFactor(const SomeInt& a, const SomeInt& b) {
  if(b >= SomeInt(0))
    return gcd(a, b);
  else
    return -gcd(a, b);
}

bool normalFactorCanReduce(const SomeInt& a) {
  return a != SomeInt(1);
}

std::string toString(const SomeInt& a) {
  return a.str();
}

constexpr int modulo = PRIME_MODULO;

class Mod {
public: /* TODO: Refactor inverse someday */
    int value; /* !!! Invariant: 0 <= value < modulo !!! */

public:
    Mod(int _constant) {
        value = ((_constant % modulo) + modulo) % modulo;
    }

    bool operator == (const Mod& m) const {
        return value == m.value;
    }

    bool operator != (const Mod& other) const {
        return value != other.value;
    }

    bool operator < (const Mod& m) const {
        return value < m.value;
    }

    void operator += (const Mod& m) {
        int v = value + m.value;
        if (v >= modulo) {
            v -= modulo;
        }
        value = v;
    }

    void operator -= (const Mod& m) {
        int v = modulo + value - m.value;
        if (v >= modulo) {
            v -= modulo;
        }
        value = v;
    }

    void operator *= (const Mod& m) {
        value = (value * m.value) % modulo;
    }

    friend std::string toString(const Mod& a) {
        if(a.value > modulo / 2)
            return std::to_string(a.value - modulo);
        return std::to_string(a.value);
    }

    friend std::ostream& operator << (std::ostream& out, const Mod &r) {
        out << toString(r);
        return out;
    }
};

std::vector<Mod> inverseTable;

void precomputeInverses() {
    inverseTable.push_back(Mod(0));
    for(int i = 1;i < modulo;i++) {
        int j = 1;
        for(;(i * j) % modulo != 1;j++);
        inverseTable.push_back(Mod(j));
    }
}

Mod inverse(const Mod& a) {
    return inverseTable[a.value];
}

Mod operator + (const Mod& a, const Mod& b) {
    Mod sum = a;
    sum += b;
    return sum;
}

Mod operator - (const Mod& a, const Mod& b) {
    Mod sum = a;
    sum -= b;
    return sum;
}

Mod operator - (const Mod& a) {
    Mod sum = 0 - a;
    return sum;
}

Mod operator * (const Mod& a, const Mod& b) {
    Mod sum = a;
    sum *= b;
    return sum;
}

Mod operator / (const Mod&a, const Mod& b) {
    return a * inverse(b);
}
