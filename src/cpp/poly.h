#include <gmpxx.h>
#include <vector>
#include <iostream>
#include <sstream>
#include <string>

template <typename T>
class Polynomial {
    public:
    std::vector<T> coeffs; 
    int length;
    int degree;
    bool monic;
    bool is_zero;

    // Standard constructor
    Polynomial(std::vector<T> c) {
        init(c);
    };

    // Construct zero polynomial
    Polynomial(int len) {
        std::vector<T> c(len);
        length = len;
        for (int i = 0; i < length; i++)
            c[i] = 0;
        init(c);
    };

    // Copy constructor
    Polynomial(const Polynomial &p){
        init(p.coeffs);
    }

    std::string str();
    void print();
    void pad(int len);
    void negate();

    private:
    void init(std::vector<T> c) {
        coeffs = c;
        length = coeffs.size();
        int i = length -1;
        while(coeffs[i] == 0 && i >= 0)
            i--;
        degree = i;
        monic = (coeffs[degree] == 1);
        is_zero = (degree == -1);
    }
};

// Pad a polynomial length to a specified length
template <typename T>
void Polynomial<T>::pad(int len) {
    T zero = 0;
    if (len > length){
        for (int i = length; i < len; i++)
            coeffs.push_back (zero);
    }
    init(coeffs);
}

// Negate polynomial
template <typename T>
void Polynomial<T>::negate() {
    for (std::size_t i = 0; i != length; i++) {
        coeffs[i] *= -1;
    }
    init(coeffs);
}

// Print to stdout
template <typename T>
void Polynomial<T>::print() {
    std::string polystr = str();
    std::cout << polystr << std::endl;
}

// Credit Jarod42 on SO
template <typename T>
std::string Polynomial<T>::str() {
    const char* plus = "";
    const char* space = "";
    
    std::stringstream buff;
    for (std::size_t i = 0; i != length; i++) {
        if (coeffs[i] == 0) {
            continue;
        }
        if (coeffs[i] > 0) {
            buff << space << plus << space;
        } else {
            buff << space << "-" << space;
        }
        plus = "+";
        space = " ";
        if (i == 0 || abs(coeffs[i]) != 1) { // to avoid to print 1x^k
            buff << abs(coeffs[i]);
        }
        if (i == 0) {
            continue;
        }
        buff << "x";
        if (i == 1) { // not print x^1
            continue;
        }
        buff << "^" << i;
    }

    std::string str = buff.str();
    return str;
}

Polynomial<mpf_class> mpz_to_mpf(Polynomial<mpz_class> p, int precision) {
    std::vector<mpf_class> c(p.length);
    for (std::size_t i = 0; i != p.length; i++) {
        c[i] = mpf_class(p.coeffs[i],precision);
    };
    Polynomial<mpf_class> q(c);
    return q;
}

Polynomial<mpq_class> mpz_to_mpq(Polynomial<mpz_class> p) {
    std::vector<mpq_class> c(p.length);
    for (std::size_t i = 0; i != p.length; i++) {
        c[i] = mpq_class(p.coeffs[i]);
    };
    Polynomial<mpq_class> q(c);
    return q;
}

Polynomial<mpz_class> mpq_to_mpz(Polynomial<mpq_class> p) {
    std::vector<mpz_class> c(p.length);
    for (std::size_t i = 0; i != p.length; i++) {
        if (mpz_class(p.coeffs[i].get_den_mpz_t()) == 1) {
            c[i] = mpz_class(p.coeffs[i]);
        } else {
            std::string polystr = p.str();
            std::string err("Could not coerce polynomial ");
            err += polystr;
            err += " to Z[x]";
            throw std::invalid_argument( err );
        }
    };
    Polynomial<mpz_class> q(c);
    return q;

}

template <typename T>
Polynomial<T> derivative(Polynomial<T> p){
    Polynomial<T> q(p.length);
    for (int i = 0; i != p.length-1; i++) {
        q.coeffs[i] = (i+1) * p.coeffs[i+1];
    }
    return q;
}

// Polynomial division: p = q*d + r. Find d and r given p and q. 
// Returns 0 if r = 0, 1 if r != 0.
// Valid for mpz_class and mpq_class Polynomials only (might produce unexpected behavior otherwise).
template <typename T>
int polydivide(const Polynomial<T> p, Polynomial<T> d, Polynomial<T> &q, Polynomial<T> &r){
    int i,j;
    int offset = p.length - d.degree;
    T leading_coeff = d.coeffs[d.degree];

    // Ensure they are all big enough
    d.pad(p.length);
    q.pad(p.length);
    r.pad(p.length);

    // Handle monic, depending on base ring
    if constexpr (std::is_same_v<T, mpz_class>){
        if (!d.monic)
            throw std::invalid_argument( "Cannot divide integer polynomial by non-monic integer polynomial." );
    }

    // Special case when d is degree 0
    if (d.degree == 0) {
        for(i = 0; i < q.length; i++){
            q.coeffs[i] = p.coeffs[i] / leading_coeff;
        }
        return 0;
    }

    // Initialize q to be reverse(p) and r to be zero
    for (i = 0; i < p.length; i++) {
        q.coeffs[i] = p.coeffs[p.length-i-1];
        r.coeffs[i] = 0;
    }

    // Synthetic division
    for (i = 0; i < offset; i++){
        if (leading_coeff != 1)
            q.coeffs[i] = q.coeffs[i] / leading_coeff;
        if (q.coeffs[i] != 0) {
            for (j = offset; j < d.length; j++){
                q.coeffs[i+j-offset+1] = q.coeffs[i+j-offset+1] - d.coeffs[d.length-j-1]*q.coeffs[i];
            }
        }
    }

    // set remainder and remove it from q ([reverse(q),reverse(r)] -> [reverse(q),0])
    for (i = 0; i < d.degree; i++) {
        r.coeffs[i] = q.coeffs[q.length - i - 1];
        q.coeffs[q.length - i - 1] = 0;
    }

    // fix quotient ([reverse(q), 0] -> [q,0])
    int start = 0;
    int end = q.length - d.degree - 1;
    T swap;
    while (start < end) {
        swap = q.coeffs[start];
        q.coeffs[start] = q.coeffs[end];
        q.coeffs[end] = swap;
        start++;
        end--;
    }
    
    return !r.is_zero;
}

// Compute gcd of two polynomials. If monic = true, then return a monic polynomial by
// dividing by the leading coefficient.
// Only valid for mpz_class and mpq_class types.
template <typename T>
Polynomial<T> gcd(const Polynomial<T> p1, const Polynomial<T> p2, bool monic = false){
    int i, count_max = p1.length;
    T leading_coeff;

    Polynomial<mpq_class> a(p1.length);
    Polynomial<mpq_class> b(p1.length);
    Polynomial<mpq_class> q(p1.length);
    Polynomial<mpq_class> r(p1.length);

    if constexpr (std::is_same_v<T, mpz_class>){
        a = mpz_to_mpq(p1);
        b = mpz_to_mpq(p2);
    } else if constexpr (std::is_same_v<T, mpq_class>){
        a = p1;
        b = p2;
    } else {
        throw std::invalid_argument( "Only mpz and mpq Polynomials are supported for `gcd` operation." );
    }

    int count = 0;
    while( polydivide(a,b,q,r) && count < count_max) { // While remainder of a/b is nonzero
        a = b;
        b = r;
    }
    if(count == count_max)
        std::cerr << "GCD error -- Euclidean algorithm did not terminate\n" << std::endl; 
    
    // set gcd to be last nonzero remainder
    Polynomial<mpq_class> gcd(p1.length);
    leading_coeff = b.coeffs[b.degree];
    for (i = 0; i < p1.length; i++) {
        gcd.coeffs[i] = b.coeffs[i];
        if (monic)
            gcd.coeffs[i] /= leading_coeff;
    }

    // Return gcd, coercing to appropriate ring.
    if constexpr (std::is_same_v<T, mpz_class>){
        return mpq_to_mpz(gcd);
    } else {
        return gcd;
    }
}
