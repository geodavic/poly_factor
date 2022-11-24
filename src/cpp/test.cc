#include "poly.h"
#include <iostream>

int main() {

		Polynomial<mpq_class> p("1+x^3+x^6");
		int precision = 200;
		mpc_t x; mpc_init2(x,precision);
		mpc_t y; mpc_init2(y,10);

		mpc_set_d_d(x,0.13,-1.023,MPC_RNDNN);
		evaluate(p,x,&y);

		mpc_out_str(stdout,10,0,y,MPC_RNDNN);

    return 0; 
}
