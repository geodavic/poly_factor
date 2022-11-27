#include "poly.h"
#include <iostream>

int main() {

		Polynomial<mpz_class> p("x^10-1");
		int precision = 200;
		mpc_t x; mpc_init2(x,precision);
		mpc_t y; mpc_init2(y,precision);

		mpc_set_d_d(x,0.13,-1.023,MPC_RNDNN);
		rootfind(p,x,y,57);

		mpc_out_str(stdout,10,0,y,MPC_RNDNN);

    return 0; 
}
