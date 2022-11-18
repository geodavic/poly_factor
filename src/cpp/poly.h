#include "gcd.h"
#include <gmpxx.h>
#include <iostream>
#include <regex>
#include <sstream>
#include <string>
#include <vector>

// TODO: handle errors here
template <typename T, typename R> T assign_value(R v, int precision) {
  if constexpr (std::is_same_v<T, mpz_class>)
    return mpz_class(v);
  if constexpr (std::is_same_v<T, mpq_class>) {
    T val;
    val = mpq_class(v);
    val.canonicalize();
    return val;
  }
  if constexpr (std::is_same_v<T, mpf_class>)
    return mpf_class(v, precision);

  throw std::invalid_argument(
      "Only mpz, mpq, and mpf class assignments are supported.");
}

template <typename T> class Polynomial {
public:
  std::vector<T> coeffs;

  // Standard constructor
  Polynomial(std::vector<T> c) { init(c); };

  // Construct zero polynomial
  Polynomial(int len) {
    std::vector<T> c(len);
    for (int i = 0; i < len; i++)
      c[i] = 0;
    init(c);
  };

  // Copy constructor
  Polynomial(const Polynomial<T> &p) { init(p.coeffs); };

  // String constructor
  Polynomial(std::string, int = 64);

  // Comparison
  bool operator==(Polynomial<T> &p) {
    if (degree() != p.degree())
      return false;
    int len = std::min(length(), p.length());
    for (int i = 0; i < len; i++) {
      if (coeffs[i] != p.coeffs[i]) {
        return false;
      }
    }
    return true;
  }

  // Properties and simple transformations
  void pad(int len);
  void negate();
  void print() const;
  std::string str() const;

  int length() const { return coeffs.size(); }

  int degree() const {
    int i = length() - 1;
    while (coeffs[i] == 0 && i >= 0)
      i--;
    return i;
  }

  bool monic() const { return (coeffs[degree()] == 1); }

  bool is_zero() const { return (degree() == -1); }

private:
  void init(std::vector<T> c) {
    if constexpr (!(std::is_same_v<T, mpz_class> |
                    std::is_same_v<T, mpq_class> |
                    std::is_same_v<T, mpf_class>))
      throw std::invalid_argument(
          "Only mpz, mpq, and mpf class Polynomials are supported.");
    coeffs = c;
  }
};

template <typename T>
Polynomial<T>::Polynomial(std::string str, int precision) {
  str = std::regex_replace(str, std::regex("\\s"), "");
  str = std::regex_replace(str, std::regex("-"), "+-");

  std::regex monomial_re(R"([x|X](\^[0-9]*)?$)");

  // Split string on '+'
  std::string segment;
  std::istringstream stream(str);
  std::vector<std::string> L;
  while (std::getline(stream, segment, '+'))
    L.push_back(segment);

  int _degree = 0;
  int d = 0;

  // Regularize terms
  for (int i = 0; i < L.size(); i++) {
    std::string term = L[i];
    if (term.size() == 0)
      continue;
    term = std::regex_replace(term, std::regex("-x"), "-1*x");
    if (term.at(0) == 'x')
      term = std::regex_replace(term, std::regex("x"), "1*x");
    if (term.back() == 'x')
      term = std::regex_replace(term, std::regex("x"), "x^1");
    if (term.find('*') == std::string::npos)
      term = std::regex_replace(term, std::regex("x"), "*x");

    L[i] = term;
    // get degree of polynomial
    d = 0;

    if (term.find('^') != std::string::npos) {
      std::smatch match;
      if (std::regex_search(term, match, monomial_re)) {
        std::string power;
        std::string m(match[0]);
        power = std::regex_replace(m, std::regex(R"(x\^)"), "");
        d = std::atoi(power.c_str());

      } else {
        std::string err("Improperly formatted polynomial term `");
        err = err + term + "`";
        throw std::invalid_argument(err);
      }
    }

    if (d > _degree)
      _degree = d;
  }

  // Create coefficients vector

  std::vector<T> coeffs(_degree + 1);
  for (int i = 0; i < coeffs.size(); i++)
    coeffs[i] = assign_value<T, int>(0, precision);

  for (int i = 0; i < L.size(); i++) {
    std::string term = L[i];
    if (term.size() == 0)
      continue;
    if (term.find('^') != std::string::npos) {
      std::smatch match;
      if (std::regex_search(term, match, monomial_re)) {
        std::string power;
        std::string m(match[0]);
        power = std::regex_replace(m, std::regex(R"(x\^)"), "");
        d = std::atoi(power.c_str());

        std::string coef = std::regex_replace(term, monomial_re, "");
        coef = std::regex_replace(coef, std::regex("\\*"), "");
        coeffs[d] += assign_value<T, std::string>(coef, precision);
      }
    } else {
      coeffs[0] += assign_value<T, std::string>(term, precision);
    }
  }

  init(coeffs);
}

// Pad a polynomial length to a specified length
template <typename T> void Polynomial<T>::pad(int len) {
  T zero = 0;
  if (len > length()) {
    for (int i = length(); i < len; i++)
      coeffs.push_back(zero);
  }
  init(coeffs);
}

// Negate polynomial
template <typename T> void Polynomial<T>::negate() {
  for (std::size_t i = 0; i != length(); i++) {
    coeffs[i] *= -1;
  }
  init(coeffs);
}

// Print to stdout
template <typename T> void Polynomial<T>::print() const {
  std::string polystr = str();
  std::cout << polystr << std::endl;
  // for (int i = 0; i != length(); i++)
  //     std::cout << coeffs[i] << ",";
  // std::cout<<std::endl;
}

// Credit Jarod42 on SO
template <typename T> std::string Polynomial<T>::str() const {
  const char *plus = "";
  const char *space = "";

  if (degree() == -1)
    return "0";

  std::stringstream buff;
  for (std::size_t i = 0; i != length(); i++) {
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
  std::vector<mpf_class> c(p.length());
  for (std::size_t i = 0; i != p.length(); i++) {
    c[i] = mpf_class(p.coeffs[i], precision);
  };
  Polynomial<mpf_class> q(c);
  return q;
}

Polynomial<mpq_class> mpz_to_mpq(Polynomial<mpz_class> p) {
  std::vector<mpq_class> c(p.length());
  for (std::size_t i = 0; i != p.length(); i++) {
    c[i] = mpq_class(p.coeffs[i]);
  };
  Polynomial<mpq_class> q(c);
  return q;
}

Polynomial<mpz_class> mpq_to_mpz(Polynomial<mpq_class> p) {
  std::vector<mpz_class> c(p.length());
  for (std::size_t i = 0; i != p.length(); i++) {
    p.coeffs[i].canonicalize();
    if (mpz_class(p.coeffs[i].get_den_mpz_t()) == 1) {
      c[i] = mpz_class(p.coeffs[i]);
    } else {
      std::string polystr = p.str();
      std::string err("Could not coerce polynomial ");
      err += polystr;
      err += " to Z[x]";
      throw std::invalid_argument(err);
    }
  };
  Polynomial<mpz_class> q(c);
  return q;
}

template <typename T> Polynomial<T> derivative(Polynomial<T> p) {
  Polynomial<T> q(p.length());
  for (int i = 0; i != p.length() - 1; i++) {
    q.coeffs[i] = (i + 1) * p.coeffs[i + 1];
  }
  return q;
}

// Polynomial division: p = q*d + r. Find d and r given p and q.
// NOTE: if over Z, this assumes p is monic.
// Returns 0 if r = 0, 1 if r != 0.
// Valid for mpz_class and mpq_class Polynomials only (might produce unexpected
// behavior otherwise).
// TODO: handle p = 0 or dd = 0
template <typename T>
int polydivide(const Polynomial<T> p, const Polynomial<T> dd, Polynomial<T> &q,
               Polynomial<T> &r) {
  int i, j;
  int offset = p.length() - dd.degree();
  T leading_coeff = dd.coeffs[dd.degree()];
  Polynomial<T> d(dd.coeffs);

  // If deg d > deg p, r = d and q = 0
  if (d.degree() > p.degree()) {
    r = d;
    q = Polynomial<T>(q.length());
    return 1;
  }

  // Ensure they are all big enough
  int len = p.length();
  d.pad(len);
  q.pad(len);
  r.pad(len);

  // Handle monic, depending on base ring
  if constexpr (std::is_same_v<T, mpz_class>) {
    if (!d.monic() || !p.monic())
      throw std::invalid_argument(
          "Division of non-monic integer polynomials is not supported.");
  }

  // Special case when d is degree 0
  if (d.degree() == 0) {
    for (i = 0; i < len; i++) {
      q.coeffs[i] = p.coeffs[i] / leading_coeff;
    }
    return 0;
  }

  // Initialize q to be reverse(p) and r to be zero
  for (i = 0; i < len; i++) {
    q.coeffs[i] = p.coeffs[len - i - 1];
    r.coeffs[i] = 0;
  }

  // Synthetic division
  for (i = 0; i < offset; i++) {
    if (leading_coeff != 1)
      q.coeffs[i] = q.coeffs[i] / leading_coeff;
    if (q.coeffs[i] != 0) {
      for (j = offset; j < len; j++) {
        q.coeffs[i + j - offset + 1] =
            q.coeffs[i + j - offset + 1] - d.coeffs[len - j - 1] * q.coeffs[i];
      }
    }
  }

  // set remainder and remove it from q ([reverse(q),reverse(r)] ->
  // [reverse(q),0])
  for (i = 0; i < d.degree(); i++) {
    r.coeffs[i] = q.coeffs[len - i - 1];
    q.coeffs[len - i - 1] = 0;
  }

  // fix quotient ([reverse(q), 0] -> [q,0])
  int start = 0;
  int end = len - d.degree() - 1;
  T swap;
  while (start < end) {
    swap = q.coeffs[start];
    q.coeffs[start] = q.coeffs[end];
    q.coeffs[end] = swap;
    start++;
    end--;
  }

  return !r.is_zero();
}

// Compute gcd of two polynomials.
// Only valid for mpz_class and mpq_class types.
template <typename T>
Polynomial<T> gcd(const Polynomial<T> p1, const Polynomial<T> p2) {
  int i, p_len = std::max(p1.length(), p2.length());

  Polynomial<mpq_class> a(p_len);
  Polynomial<mpq_class> b(p_len);
  Polynomial<mpq_class> q(p_len);
  Polynomial<mpq_class> r(p_len);

  // Make sure nonzero
  if (p1.degree() < 0 || p2.degree() < 0)
    throw std::invalid_argument("Cannot compute gcd of zero polynomials.");

  // Cast to rationals
  if constexpr (std::is_same_v<T, mpz_class>) {
    a = mpz_to_mpq(p1);
    b = mpz_to_mpq(p2);
    r = b;
  } else if constexpr (std::is_same_v<T, mpq_class>) {
    a = p1;
    b = p2;
    r = b;
  } else {
    throw std::invalid_argument(
        "Only mpz and mpq Polynomials are supported for `gcd` operation.");
  }

  // Divide out leading coefficients (will be added at end)
  b.pad(p_len);
  a.pad(p_len);
  for (i = 0; i < p_len; i++) {
    a.coeffs[i] = a.coeffs[i] / a.coeffs[a.degree()];
    b.coeffs[i] = b.coeffs[i] / b.coeffs[b.degree()];
  }

  // If a < b, swap them so a > b
  if (a.degree() < b.degree()) {
    r = a;
    a = b;
    b = r;
  }

  int count = 0;
  while (polydivide(a, b, q, r) &&
         count < p_len) { // While remainder of a/b is nonzero
    count++;
    a = b;
    b = r;
  }
  if (count == p_len)
    std::cerr << "GCD error -- Euclidean algorithm did not terminate\n"
              << std::endl;

  // set gcd to be last nonzero remainder
  Polynomial<mpq_class> g(p_len);
  mpq_class leading_coeff = b.coeffs[b.degree()];
  for (i = 0; i < p_len; i++)
    g.coeffs[i] =
        b.coeffs[i] / leading_coeff; // Force to be monic, this is okay since a
                                     // and b are assumed to be monic

  // If over Z, clear denominators of gcd and multiply back in gcd of
  // coefficients of p1 and p2
  if constexpr (std::is_same_v<T, mpz_class>) {
    // Find lcm of denominators of g
    std::vector<mpz_class> denom(g.length());
    for (i = 0; i < g.length(); i++)
      denom[i] = mpz_class(g.coeffs[i].get_den_mpz_t());
    mpq_class lcm_denom = mpq_class(lcm(denom));

    // Find gcd of p1 and p2 coeffs
    std::vector<mpz_class> all_coeffs;
    for (i = 0; i < p1.length(); i++) {
      if (p1.coeffs[i] != 0)
        all_coeffs.push_back(abs(p1.coeffs[i]));
    }
    for (i = 0; i < p2.length(); i++) {
      if (p2.coeffs[i] != 0)
        all_coeffs.push_back(abs(p2.coeffs[i]));
    }
    mpq_class gcd_coeffs = mpq_class(integer_gcd_v(all_coeffs));

    // Clear denominators and multiply gcd back in
    for (i = 0; i < g.length(); i++)
      g.coeffs[i] *= gcd_coeffs * lcm_denom;
  }

  // Return gcd, coercing to appropriate ring.
  if constexpr (std::is_same_v<T, mpz_class>) {
    return mpq_to_mpz(g);
  } else {
    return g;
  }
}
