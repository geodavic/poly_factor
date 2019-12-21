#include <stdio.h>
#include <stdlib.h>
#include "../../mpz_algebraic.h"
#include "../poly_functions.h"

int main(){
	int i,poly_len=6;
	mpz_t *p=malloc(poly_len*sizeof(mpz_t)); //polynomial
	mpz_t *pp=malloc(poly_len*sizeof(mpz_t)); //other polynomial
	mpz_t *r=malloc(poly_len*sizeof(mpz_t)); //gcd
	
	for(i=0;i<poly_len;i++){
		mpz_init(p[i]);
		mpz_init(pp[i]);
		mpz_init(r[i]);
		mpz_set_ui(r[i],0);
	}

	mpz_set_si(p[0],-1);
	mpz_set_si(p[1],0);
	mpz_set_si(p[2],0);
	mpz_set_si(p[3],0);
	mpz_set_si(p[4],0);
	mpz_set_si(p[5],1);

	mpz_set_si(pp[0],0);
	mpz_set_si(pp[1],0);
	mpz_set_si(pp[2],0);
	mpz_set_si(pp[3],0);
	mpz_set_si(pp[4],5);
	mpz_set_si(pp[5],0);
		
	gcd(pp,p,r,poly_len);
	
	print_poly(poly_len,r,1);

	return 0;
}

