#include "poly.h"
#include <cassert>
#include <iostream>

Polynomial<mpz_class> P1("1");
Polynomial<mpz_class> P2("x - 1");
Polynomial<mpz_class> P3("x^5 - 1");
Polynomial<mpz_class> P4("x^5 + 3x - 4");
Polynomial<mpz_class> P5("8*x^3 - 6*x + 2");
Polynomial<mpz_class> P6("24x^2 - 6");
Polynomial<mpz_class> P10("20x^3 + 40");
Polynomial<mpz_class> P11("15x^2 + 30");
Polynomial<mpz_class>
    P12("x^8 + x^7 - 2 x^6 + 6 x^5 - 18 x^4 - 4 x^3 - 19 x^2 - 23 x - 6");
Polynomial<mpz_class>
    P13("-23 - 38 x - 12 x^2 - 72 x^3 + 30 x^4 - 12 x^5 + 7 x^6 + 8 x^7");
Polynomial<mpz_class>
    P14("204 x^8 + 17 x^7 - 83 x^6 - 189 x^5 - 37 x^4 + 163 x^3 + 9 x - 48");
Polynomial<mpz_class> P15("51 x^5 + 17 x^4 - 8 x^3 - 45 x^2 - 19 x + 16");

Polynomial<mpq_class> Q1("3/4x^5 - 6/7");
Polynomial<mpq_class> Q2("2/3x^2 - 3/2x + 1/5");

Polynomial<mpf_class> R1("1.1x^5 - 2.3");

bool test_string_and_constructor();
bool test_derivative();
bool test_divide();
bool test_gcd();

int main() {

  bool result = true;

  std::cout << " ---------------------- Poly tests ---------------------------"
            << std::endl;

  result &= test_string_and_constructor();
  result &= test_derivative();
  result &= test_divide();
  result &= test_gcd();

  if (!result) {
    std::cout
        << "---------------------- Poly tests FAILED ----------------------"
        << std::endl;
    return 1;
  }
}

bool test_string_and_constructor() {
  std::cout << "test_string_and_constructor() ... ";

  Polynomial<mpz_class> P1copy(P1.str());
  Polynomial<mpz_class> P2copy(P2.str());
  Polynomial<mpz_class> P5copy(P5.str());
  Polynomial<mpq_class> Q1copy(Q1.str());
  Polynomial<mpf_class> R1copy(R1.str());

  assert(P1 == P1copy);
  assert(P2 == P2copy);
  assert(P5 == P5copy);
  assert(Q1 == Q1copy);
  assert(R1 == R1copy);

  std::cout << "passed" << std::endl;
  return true;
}

bool test_derivative() {
  std::cout << "test_derivative() ... ";

  Polynomial<mpz_class> P5prime = derivative(P5);
  Polynomial<mpz_class> P12prime = derivative(P12);

  assert(P5prime == P6);
  assert(P12prime == P13);

  std::cout << "passed" << std::endl;
  return true;
}

// TODO: add division of smaller by bigger test case
bool test_divide() {
  std::cout << "test_divide() ... ";

  Polynomial<mpz_class> r(P3.length());
  Polynomial<mpz_class> q(P3.length());
  Polynomial<mpz_class> quotient3("x^4+x^3+x^2+x+1");
  Polynomial<mpz_class> remainder3("0");
  int rzero1 = polydivide(P3, P2, q, r);
  assert(quotient3 == q);
  assert(quotient3 == q);
  assert(remainder3 == r);
  assert(rzero1 == 0);

  Polynomial<mpq_class> rq(Q1.length());
  Polynomial<mpq_class> qq(Q1.length());
  Polynomial<mpq_class> quotient7(
      "9/8x^3 + 81/32x^2 + 3429/640*x + 28917/2560");
  Polynomial<mpq_class> remainder7("406323/25600*x - 279219/89600");
  int rzero2 = polydivide(Q1, Q2, qq, rq);
  assert(quotient7 == qq);
  assert(remainder7 == rq);
  assert(rzero2 == 1);

  std::cout << "passed" << std::endl;
  return true;
}

bool test_gcd() {
  std::cout << "test_gcd() ... ";

  Polynomial<mpz_class> G1 = gcd(P5, P6);
  Polynomial<mpz_class> GG1("4x-2");
  assert(G1 == GG1);

  Polynomial<mpz_class> G2 = gcd(P6, P5);
  Polynomial<mpz_class> GG2("4x-2");
  assert(G2 == GG2);

  Polynomial<mpz_class> G3 = gcd(P5, P3);
  Polynomial<mpz_class> GG3("1");
  assert(G3 == GG3);

  Polynomial<mpq_class> G4 = gcd(mpz_to_mpq(P5), mpz_to_mpq(P6));
  Polynomial<mpq_class> GG4("x-1/2");
  assert(G4 == GG4);

  Polynomial<mpz_class> G5 = gcd(P10, P11);
  Polynomial<mpz_class> GG5("5");
  assert(G5 == GG5);

  Polynomial<mpz_class> G6 = gcd(P3, P4);
  Polynomial<mpz_class> GG6("x-1");
  assert(G6 == GG6);

  Polynomial<mpz_class> G7 = gcd(P12, P13);
  Polynomial<mpz_class> GG7("x^3+2x+1");
  assert(G7 == GG7);

  Polynomial<mpz_class> G8 = gcd(P14, P15);
  Polynomial<mpz_class> GG8("17x^3+3x-16");
  assert(G8 == GG8);

  std::cout << "passed" << std::endl;
  return true;
}
