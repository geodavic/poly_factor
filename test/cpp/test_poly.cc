#include "poly.h"
#include <iostream>
#include <cassert>

Polynomial<mpz_class> P1("1");
Polynomial<mpz_class> P2("x - 1");
Polynomial<mpz_class> P3("x^5 - 1");
Polynomial<mpz_class> P5("8x^3 - 6x + 2");
Polynomial<mpz_class> P6("24x^2 - 6");
Polynomial<mpq_class> P7("3/3x^2 - 6/7");
Polynomial<mpf_class> P8("1.1x^2 - 2.3");

bool test_string_and_constructor();

int main() {
    
    bool result = true;
    
    std::cout << " ---------------------- Poly tests ---------------------------" << std::endl;

    result &= test_string_and_constructor();

    if (!result) {
        std::cout << "----------------- Poly tests FAILED -          ------" << std::endl;
        return 1;
    }
}

bool test_string_and_constructor() {
    std::cout << "test_string_and_constructor() ... ";
    bool result = true; 

    Polynomial<mpz_class> P1copy(P1.str());
    Polynomial<mpz_class> P2copy(P2.str());
    Polynomial<mpz_class> P5copy(P5.str());
    Polynomial<mpq_class> P7copy(P7.str());
    Polynomial<mpf_class> P8copy(P8.str());

    result &= (P1 ==P1copy);
    result &= (P2 ==P2copy);
    result &= (P5 ==P5copy);
    result &= (P7 ==P7copy);
    result &= (P8 ==P8copy);

    if (result)
        std::cout << "passed" << std::endl;
    return result;
}

