#include <gmpxx.h>
#include <vector>

// GCD over integers (Euclidean alg)
mpz_class integer_gcd(const mpz_class p, const mpz_class d) {
    mpz_class q;
    mpz_class r;
    mpz_class a = p;
    mpz_class b = d;

    if (p < d){
        r = a;
        a = b;
        b = r;
    }
    r = b;

    // Special case
    if (b<2)
        return mpz_class(1);

    while (true) {
        q = 0;
        while (q*b <= a)
            q++;

        if (a == (q-1)*b)
            break;

        r = a - (q-1)*b;
        a = b;
        b = r;
    }
    return r;
}

// gcd of list of positive integers
mpz_class integer_gcd_v(const std::vector<mpz_class> nums) {
    mpz_class g = nums[0];
    int pointer = 1;
    while (pointer < nums.size()) {
        g = integer_gcd(g,nums[pointer]); 
        pointer ++;
    }
    return g;
}

// lcm of list of positive integers
mpz_class lcm(const std::vector<mpz_class> nums) {
    mpz_class min, lcm = 0;
    bool all_equal = true;
    int i, min_idx;
    std::vector<mpz_class> nums_c = nums;

    for (i=1; i<nums_c.size(); i++)
        all_equal &= (nums_c[i]==nums_c[i-1]); 

    while (!all_equal) {
        min = nums_c[0];
        min_idx = 0;
        for (i=0; i<nums_c.size(); i++){
            if (nums_c[i] < min){
                min = nums_c[i];
                min_idx = i;
            }
        }
        nums_c[min_idx] += nums[min_idx];

        all_equal = true;
        for (i=1; i<nums_c.size(); i++)
            all_equal &= (nums_c[i]==nums_c[i-1]); 
    }
    return nums_c[0];
}
