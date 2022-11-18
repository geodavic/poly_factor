#include "poly.h"
#include <iostream>

int main() {

    std::vector<mpz_class> nums{20,60,16,24};
    mpz_class g = lcm(nums);

    std::cout << g << std::endl;

    return 0; 
}
