#include "poly.h"
#include <iostream>

int main() {
   
    std::vector<mpz_class> c1 {1,4,6,4,1};
    Polynomial<mpz_class> p(c1);
    Polynomial<mpz_class> pp = derivative(p);

    p.print();
    pp.print();

    bool monic = true;
    Polynomial<mpz_class> g = gcd(p,pp,monic=monic);
    g.print();

    return 0; 
}
