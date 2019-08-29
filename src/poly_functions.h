//Function library for all polynomial factorization files (these files must also include mpz_algebraic.h, which this references extensively)

//----major changes------------//
//GDT 01.2018
//    05.2018
//    07.2018
//---------about---------------//


//------------TODO------------------//
//implement irreducibility tests
//			- Eisenstein (plus shifts?)
//      - Perron's criterion
//      - Murdy, Osada, Bauer criteria
//      - reverting
//      - Newton polygons?

void evaluate(mpz_t *p, int len,const mpf_t input, mpf_t output,int PRECISION);
void evaluate_cx(mpz_t *p, int len, const mpc_t input, mpc_t output, int PRECISION);
int degree(mpz_t *p, int len);
int degree_q(mpq_t *p, int len);
int rootfind(mpz_t *p, int len, mpf_t start,mpf_t root,int log10_thresh, int PRECISION);
int rootfind_cx(mpz_t *p, int len, mpc_t start,mpc_t root,int log10_thresh, int PRECISION);
int polydivide(mpz_t *p,mpz_t *d,mpz_t *out,int len);
int polydivide_r(mpq_t *p,mpq_t *d,mpq_t *r,int len);
void gcd(mpz_t *poly1, mpz_t *poly2, mpz_t *gcd, int poly_len);
int find_factor(mpz_t *poly, mpz_t *d, mpz_t *q, int poly_len,int PRECISION);
int find_factor_cx(mpz_t *poly, mpz_t *d, mpz_t *q, int poly_len,int PRECISION,int verbosity,double d_delta);
int factorize(mpz_t *poly,int poly_len,int PRECISION,mpz_t *factors,int verbosity,double delta);
int factorize_full(mpz_t *poly,int poly_len,int PRECISION,mpz_t *factors,int verbosity, double delta);
int monic_slide(int len, mpz_t *p);
int monic_slide_dont_multiply(int len, mpz_t *p);
void derivative(mpz_t *p,mpz_t *pp,int poly_len);

//evaluate p(x) at input using Horner's method
void evaluate(mpz_t *p, int len,const mpf_t input, mpf_t output, int PRECISION){
	mpf_t bi; mpf_init2(bi,PRECISION);
	mpf_t dummy; mpf_init2(dummy,PRECISION);
	int i;	
	mpf_set_z(bi,p[len-1]); //bi=an 
	for(i=len-2;i>=0;i--){ 
		//set bi=ai+a(i+1)*input
		mpf_mul(bi,bi,input); 
		mpf_set_z(dummy,p[i]);
		mpf_add(bi,dummy,bi);
	}
	mpf_set(output,bi); //f(input) = b0
	mpf_clear(bi);
	mpf_clear(dummy);
}

//evaluate p(x), where the input is a complex number. Also Horner's method
void evaluate_cx(mpz_t *p, int len, const mpc_t input, mpc_t output, int PRECISION){
	mpc_t bi; mpc_init2(bi,PRECISION);
	mpc_t dummy; mpc_init2(dummy,PRECISION);
	int i;
	mpc_set_z(bi,p[len-1],MPC_RNDNN); //valgrind takes issue with this
	for(i=len-2;i>=0;i--){
		//set bi=ai+a(i+1)*input
		mpc_mul(bi,bi,input,MPC_RNDNN);
		mpc_set_z(dummy,p[i],MPC_RNDNN); //or this
		mpc_add(bi,dummy,bi,MPC_RNDNN);
	}
	mpc_set(output,bi,MPC_RNDNN);
	mpc_clear(dummy);
	mpc_clear(bi);
}

//degree of p(x) (in Z[x])
//considers the zero polynomial to have degree -1
int degree(mpz_t *p, int len){
	int i=len-1;
	while(mpz_sgn(p[i])==0&&i>=0)
		i--;
	return i;
}

//degree of p(x) (in Q[x])
//considers the zero polynomial to have degree -1
int degree_q(mpq_t *p, int len){
	int i=len-1;
	while(mpq_sgn(p[i])==0&&i>=0)
		i--;
	return i;
}

//find one real root of p(x) using second order Newton's method (Halley's method). 
//cubic convergence, stop when log10(|xn-x(n+1)|)<-log10_thresh 
//returns 1 on success and 0 on failure (i.e. exceeding a certain amount of iterations without getting within the threshold of zero).
int rootfind(mpz_t *p, int len, mpf_t start,mpf_t root,int log10_thresh, int PRECISION){
	int i;
	int max_iterates=50; //maximum number of iterations allowed. 50 is generous
	mpf_t diff; mpf_init2(diff,PRECISION);
	mpf_t quot; mpf_init2(quot,PRECISION);
	mpf_t dummy; mpf_init2(dummy,PRECISION);
	mpf_t thresh; mpf_init2(thresh,PRECISION);
	mpz_t *pp; pp=malloc((len)*sizeof(mpz_t)); //p'(x)
	mpz_t *ppp; ppp=malloc((len)*sizeof(mpz_t)); //p''(x)
	mpf_t eval_p; mpf_init2(eval_p,PRECISION);
	mpf_t eval_pp; mpf_init2(eval_pp,PRECISION);
	mpf_t eval_ppp; mpf_init2(eval_ppp,PRECISION);

	//threshold for stopping
	mpf_set_ui(thresh,10);
	mpf_pow_ui(thresh,thresh,log10_thresh);
	mpf_ui_div(thresh,1,thresh);


	//compute coefficients of p'(x) and p''(x)
	for(i=0;i<(len);i++){
		mpz_init(pp[i]); 
		mpz_init(ppp[i]); 
	}
	for(i=0;i<(len-1);i++)
		mpz_mul_ui(pp[i],p[i+1],i+1);
	for(i=0;i<(len-2);i++)
		mpz_mul_ui(ppp[i],pp[i+1],i+1);

	mpf_set_ui(diff,1);
	mpf_set(root,start);
	int c=0; //counter for loop
	//main loop
	while(mpf_cmp(diff,thresh)>0&&(c<max_iterates)){
		c++;
		evaluate(p,len,root,eval_p,PRECISION);//evaluate p(xn)
		evaluate(pp,len-1,root,eval_pp,PRECISION); //evaluate pp(xn)
		evaluate(ppp,len-2,root,eval_ppp,PRECISION);//evaluate ppp(xn)

		mpf_set_ui(quot,2);
		mpf_mul(quot,quot,eval_pp);
		mpf_mul(quot,quot,eval_pp);
		mpf_mul(dummy,eval_p,eval_ppp);
		mpf_sub(quot,quot,dummy); //denominator of halley quotient
		mpf_set_ui(dummy,2);
		mpf_mul(dummy,dummy,eval_p);
		mpf_mul(dummy,dummy,eval_pp); //numerator
		mpf_div(quot,dummy,quot); //halley quotient

		mpf_abs(diff,quot); //diff = |quot|
		mpf_sub(root,root,quot); //x(n+1) = xn - quot

		//progress prints
		//mpf_out_str(stdout,10,0,root);
		//printf("\n");
	}

	//clear variables
	for(i=0;i<(len);i++){
		mpz_clear(pp[i]);
		mpz_clear(ppp[i]);
	}
	free(pp);
	free(ppp);
	mpf_clear(diff);
	mpf_clear(eval_p);
	mpf_clear(eval_pp);
	mpf_clear(eval_ppp);
	mpf_clear(quot);
	mpf_clear(dummy);
	mpf_clear(thresh);
	if(c==max_iterates) //didn't find root 
		return 0;
	else //found root
		return 1;
}

//same as above, except for complex numbers
//note: start value should not be totally real, since this iteration sends reals to reals
//    : make sure p has no repeated roots (otherwise this isn't guaranteed to converge)
//TODO: make exit condition depend on abs(Re(x)) and abs(Im(x))
int rootfind_cx(mpz_t *p, int len, mpc_t start,mpc_t root,int log10_thresh, int PRECISION){
	int i;
	int max_iterates=100;
	mpfr_t diff; mpfr_init2(diff,PRECISION);
	mpc_t quot; mpc_init2(quot,PRECISION);
	mpc_t dummy; mpc_init2(dummy,PRECISION);
	mpfr_t thresh; mpfr_init2(thresh,PRECISION);
	mpz_t *pp; pp=malloc((len)*sizeof(mpz_t));
	mpz_t *ppp; ppp=malloc((len)*sizeof(mpz_t));
	mpc_t eval_p; mpc_init2(eval_p,PRECISION);
	mpc_t eval_pp; mpc_init2(eval_pp,PRECISION);
	mpc_t eval_ppp; mpc_init2(eval_ppp,PRECISION);

	mpfr_set_ui(thresh,10,MPFR_RNDN);
	mpfr_pow_ui(thresh,thresh,log10_thresh,MPFR_RNDN);
	mpfr_ui_div(thresh,1,thresh,MPFR_RNDN);

	for(i=0;i<(len);i++){
		mpz_init(pp[i]); 
		mpz_init(ppp[i]); 
	}
	for(i=0;i<(len-1);i++)
		mpz_mul_ui(pp[i],p[i+1],i+1);
	for(i=0;i<(len-2);i++)
		mpz_mul_ui(ppp[i],pp[i+1],i+1);
	
	mpfr_set_ui(diff,1,MPFR_RNDN);
	mpc_set(root,start,MPC_RNDNN);
	int c=0;
	while(mpfr_cmp(diff,thresh)>=0&&(c<max_iterates)){
		c++;
		evaluate_cx(p,len,root,eval_p,PRECISION);
		evaluate_cx(pp,len-1,root,eval_pp,PRECISION);
		evaluate_cx(ppp,len-2,root,eval_ppp,PRECISION);

		mpc_set_ui(quot,2,MPC_RNDNN);
		mpc_mul(quot,quot,eval_pp,MPC_RNDNN);
		mpc_mul(quot,quot,eval_pp,MPC_RNDNN);
		mpc_mul(dummy,eval_p,eval_ppp,MPC_RNDNN);
		mpc_sub(quot,quot,dummy,MPC_RNDNN);
		mpc_set_ui(dummy,2,MPC_RNDNN);
		mpc_mul(dummy,dummy,eval_p,MPC_RNDNN);
		mpc_mul(dummy,dummy,eval_pp,MPC_RNDNN);
		mpc_div(quot,dummy,quot,MPC_RNDNN);

		mpc_abs(diff,quot,MPC_RNDNN);
		mpc_sub(root,root,quot,MPC_RNDNN);
		//mpc_out_str(stdout,10,0,root,MPC_RNDNN);
		//printf("\n");
	}

	for(i=0;i<len;i++){
		mpz_clear(pp[i]);
		mpz_clear(ppp[i]);
	}
	free(pp);
	free(ppp);
	mpfr_clear(diff);
	mpc_clear(eval_p);
	mpc_clear(eval_pp);
	mpc_clear(eval_ppp);
	mpc_clear(quot);
	mpc_clear(dummy);
	mpfr_clear(thresh);
	if(c==max_iterates)
		return 0;
	else
		return 1;
}


//divide polynomial p by polynomial d, i.e. compute quotient q in p=d*q+r
//everything monic
//each poly should be allocated to 'len' size
//returns 0 if r=0, 1 otherwise
int polydivide(mpz_t *p,mpz_t *d,mpz_t *out,int len){
	int i,j,offset,deg_d,return_val;
	mpz_t dummy;mpz_init(dummy);

	//populate out. starts off as reverse(p), becomes [reverse(q),reverse(r)] after division. 
	//reverse because it makes synthetic division easier
	for(i=0;i<len;i++)
		mpz_set(out[i],p[len-i-1]);
	

	//offset index in d, (also the max index of loop below)
	deg_d=degree(d,len);
	offset=len-deg_d;

	//synthetic division loop
	for(i=0;i<offset;i++){
		mpz_set(dummy,out[i]);
		if(mpz_sgn(dummy)!=0){
			for(j=offset;j<len;j++)
				mpz_submul(out[i+j-offset+1],d[len-j-1],dummy); //out[i+j-offset+1] -= reverse(d)[j]*dummy
		}
	}

	//set return value
	i=0;
	while(mpz_sgn(out[len-1-i])==0&&i<=deg_d)
		i++;
	if(i==deg_d) //if r=0
		return_val=0;
	else
		return_val=1;

	//turn out into q (forgetting r): [reverse(q),reverse(r)] -> [q,0]
	for(i=len-deg_d;i<len;i++)
		mpz_set_ui(out[i],0); //set r=0
	int start=0;
	int end=len-deg_d-1;	
	while(start<end){ //reverse reverse(q)
		mpz_set(dummy,out[start]);
		mpz_set(out[start],out[end]);
		mpz_set(out[end],dummy);
		start++;
		end--;
	}

	//clear variables
	mpz_clear(dummy);
	return return_val;
}

//same as above, except returns remainder instead of q
// p = d*q + r  (r is in Q[x])
// all polynomials are rational, since none are assumed to be monic
// return 0 if r=0, 1 otherwise
int polydivide_r(mpq_t *p,mpq_t *d,mpq_t *r,int len){
	int i,j,offset,deg_d,return_val;
	mpq_t dummy;mpq_init(dummy);
	mpq_t dummy2;mpq_init(dummy2);
	mpq_t lc;mpq_init(lc);
	mpq_t *outq=malloc(len*sizeof(mpq_t)); 

	//populate out. starts off as reverse(p), becomes [reverse(q),reverse(r)] after division. 
	//reverse because it makes synthetic division easier
	for(i=0;i<len;i++){
		mpq_init(outq[i]);
		mpq_set(outq[i],p[len-i-1]);
		mpq_set_ui(r[i],0,1); //make sure r=0 to start
	}
	

	//offset index in d, (also the max index of loop below)
	deg_d=degree_q(d,len);
	offset=len-deg_d;

	//edge case (if d is degree zero)
	if(deg_d==0){
		for(i=0;i<len;i++)
			mpq_clear(outq[i]);
		free(outq);
		return 0; //can always divide by degree zero divisor
	//edge case (if d=0)
	}
	if(deg_d==-1){
		printf("warning! tried to divide by zero polynomial.\n");
		for(i=0;i<len;i++)
			mpq_clear(outq[i]);
		free(outq);
		return 0; //can't divide by zero
	}

	//leading coefficient of d
	mpq_set(lc,d[deg_d]);
	
	//synthetic division loop
	for(i=0;i<offset;i++){
		mpq_div(outq[i],outq[i],lc);
		mpq_set(dummy,outq[i]);
		if(mpq_sgn(dummy)!=0){
			for(j=offset;j<len;j++){
				//out[i+j-offset+1] -= reverse(d)[j]*dummy
				mpq_mul(dummy2,d[len-j-1],dummy);
				mpq_sub(outq[i+j-offset+1],outq[i+j-offset+1],dummy2);
				//mpz_submul(out[i+j-offset+1],d[len-j-1],dummy);
			}
		}
	}

	//set r
	for(i=0;i<deg_d;i++){
		mpq_set(r[i],outq[len-i-1]);
	}
	
	//set return value
	if(degree_q(r,len)>=0) //if r is nonzero
		return_val=1;
	else
		return_val=0;

	//clear variables
	for(i=0;i<len;i++){
		mpq_clear(outq[i]);
	}
	free(outq);	
	mpq_clear(dummy);
	mpq_clear(dummy2);
	mpq_clear(lc);

	return return_val;
}

//compute gcd using synthetic division
//this algorithm coerces poly1,poly2 into Q[x] and then does a gcd 
//algorithm over the field Q. What results is generically a non-integer gcd; this
//can then be coerced into Z[x] by an appropriate multiple
//for now, this assumes that the gcd will be monic, based on how it is used in factorize()
//hence the appropriate multiple is the inverse of the leading coefficient
void gcd(mpz_t *poly1, mpz_t *poly2, mpz_t *gcd, int poly_len){
	int i,count_max=poly_len;
	mpq_t lc; mpq_init(lc);
	mpq_t *a=malloc(poly_len*sizeof(mpq_t));//temporary polynomials
	mpq_t *b=malloc(poly_len*sizeof(mpq_t));
	mpq_t *r=malloc(poly_len*sizeof(mpq_t));
	for(i=0;i<poly_len;i++){
		mpq_init(a[i]);
		mpq_init(b[i]);
		mpq_init(r[i]);
		mpq_set_z(a[i],poly1[i]); 
		mpq_set_z(b[i],poly2[i]);
	}

	int count=0;
	while(polydivide_r(a,b,r,poly_len)&&count<count_max){//while remainder of a/b is nonzero
		//set a <- b
		//    b <- r
		count++;
		for(i=0;i<poly_len;i++){
			mpq_set(a[i],b[i]);
			mpq_set(b[i],r[i]);
		}
	}
	if(count==count_max)
		printf("gcd error! Euclidean division did not terminate\n");
	

	//set gcd to last nonzero remainder
	//(make monic first)
	int degr=degree_q(b,poly_len);
	mpq_set(lc,b[degr]); //leading coefficient of gcd, to be divided out
	for(i=0;i<poly_len;i++){// set gcd = b[i]/lc
		mpq_div(b[i],b[i],lc);
		mpq_get_num(gcd[i],b[i]);
	}

	//clear variables
	for(i=0;i<poly_len;i++){
		mpq_clear(a[i]);
		mpq_clear(b[i]);
		mpq_clear(r[i]);
	}
	mpq_clear(lc);
	free(a);
	free(b);
	free(r);

}


//find irreducible factor of poly, poly=d*q. Return zero if no factors found, return 1 if factor is found
//Warning: this sets d,q to zero upon failure.
//only finds real roots (for a slight speedup if that's all that is needed)- see below for more general version
//This function is now slightly outdated. See find_factor_cx for most recent parameters and tweaks
int find_factor(mpz_t *poly, mpz_t *d, mpz_t *q, int poly_len,int PRECISION){
	int i,j;
	int sig_digits,deg,input_degree=poly_len-1;
	int LLL_found_divisor=0;
	int log10thresh=(int)(PRECISION*log10(2.0)); //closest we can get to root with given PRECISION
	mpf_t input; mpf_init2(input,PRECISION);
	mpf_t output; mpf_init2(output,PRECISION);
	mpf_t delta; mpf_init2(delta,PRECISION);mpf_set_d(delta,0.75);//LLL parameter
	mpf_t thresh; mpf_init2(thresh,PRECISION); //10^(-log10thresh)
	mpf_t dummy; mpf_init2(dummy,PRECISION); //dummy variables
	mpz_t dummy_z; mpz_init(dummy_z);

	//make sure d,q are zeroed out
	for(i=0;i<poly_len;i++){
		mpz_set_ui(d[i],0);
		mpz_set_ui(q[i],0);
	}

	printf("Finding a factor of:\n");
	print_poly(poly_len,poly);

	//find a root
	mpf_set_d(input,3.1415926); //starting value for rootfind
	if(!rootfind(poly,poly_len,input,output,log10thresh,PRECISION)){
		printf("No root found\n");
		//clear variables and exit
		mpf_clear(input);
		mpf_clear(output);
		mpf_clear(delta);
		mpf_clear(thresh);
		mpf_clear(dummy);
		mpz_clear(dummy_z);
		return 0;
	}

	printf("root chosen: ");
	mpf_out_str(stdout,10,15,output);

	//check if its an integer. If it is, we're done
	mpf_set_ui(thresh,10);
	mpf_pow_ui(thresh,thresh,log10thresh);
	mpf_ui_div(thresh,1,thresh); //set value of thresh from log10thresh
	Mpf_round_f(dummy,output,PRECISION);
	mpf_sub(dummy,dummy,output); 
	mpf_abs(dummy,dummy);//dummy = |difference between output and nearest integer|
	if(mpf_cmp(dummy,thresh)<=0){//integer check
		printf(" (integer)");
		Mpf_round(dummy_z,output,PRECISION); //round output to mpz
		mpz_mul_si(d[0],dummy_z,-1); //set d(x) = x-output, since output is integer
		mpz_set_si(d[1],1);
		//check that d is a divisor
		if(polydivide(poly,d,q,poly_len)==0){
			printf("Factor:\n");
			print_poly(poly_len,d);
			printf("Quotient:\n");
			print_poly(poly_len,q);
			//clear variables
			mpf_clear(input);
			mpf_clear(output);
			mpf_clear(delta);
			mpf_clear(thresh);
			mpf_clear(dummy);
			mpz_clear(dummy_z);
			return 1;
		}
	}
	printf("\n");

	//otherwise run LLL on each degree less than input degree to find minimal polynomial
	mpz_t *basis;//initialize basis
	basis=malloc((input_degree+1)*(input_degree+2)*sizeof(mpz_t));
	for(i=0;i<(input_degree+1)*(input_degree+2);i++)
		mpz_init(basis[i]);

	for(deg=2;deg<=input_degree;deg++){//loop on degrees
		printf("      LLL on degree %d\n",deg);
		sig_digits=sig_mpf(output,deg,PRECISION);
		//find irreducible polynomial for chosen root
		create_basis(basis,output,deg,sig_digits,PRECISION);
		LLL(deg+2,deg+1,basis,delta,PRECISION);
		for(j=0;j<deg+1;j++)
			mpz_set(d[j],basis[j]); //set first vector of reduced basis to divisor d
		
		//TODO: pick shortest vector instead of first one
		
		//check if p is monic and perform synthetic division
		if((monic_slide(deg+1,d)>=0)&&(polydivide(poly,d,q,poly_len)==0)){
			LLL_found_divisor=1;
			printf("Factor:\n");
			print_poly(poly_len,d);
			printf("Quotient:\n");
			print_poly(poly_len,q);
			break; //quit once you've found lowest degree divisor
		}
	}

	if(LLL_found_divisor==0){
		//no divisor found by LLL, clear variables and exit
		mpf_clear(input);
		mpf_clear(output);
		mpf_clear(delta);
		mpf_clear(thresh);
		mpf_clear(dummy);
		mpz_clear(dummy_z);
		for(i=0;i<(input_degree+1)*(input_degree+2);i++)
			mpz_clear(basis[i]);
		free(basis);
		return 0;
	}
	
	//clear variables
	mpf_clear(input);
	mpf_clear(output);
	mpf_clear(delta);
	mpf_clear(thresh);
	mpf_clear(dummy);
	mpz_clear(dummy_z);
	for(i=0;i<(input_degree+1)*(input_degree+2);i++)
		mpz_clear(basis[i]);
	free(basis);
	return 1;
}

//same as above, but allows for complex roots (more general, maybe requires /slightly/ more precision)
//notes: - might be able to reduce down to at most one dummy variable of each data type
int find_factor_cx(mpz_t *poly, mpz_t *d, mpz_t *q, int poly_len,int PRECISION,int verbosity, double d_delta){
	int i,j,iter=0,iter_max=3;
	int sig_digits,deg,input_degree=poly_len-1;
	int LLL_found_divisor=0;
	int log10thresh=(int)(PRECISION*log10(2.0)); //closest we can get to root with given PRECISION
	gmp_randstate_t seed; gmp_randinit_default(seed);//seed for random starting value of rootfind
	mpc_t input; mpc_init2(input,PRECISION);
	mpc_t output; mpc_init2(output,PRECISION);
	mpfr_t output_r; mpfr_init2(output_r,PRECISION);
	mpfr_t output_i; mpfr_init2(output_i,PRECISION);
	mpf_t delta; mpf_init2(delta,PRECISION);mpf_set_d(delta,d_delta);//LLL parameter
	mpfr_t thresh; mpfr_init2(thresh,PRECISION); //10^(-log10thresh)
	mpfr_t dummy; mpfr_init2(dummy,PRECISION); //dummy variables
	mpfr_t dummy2; mpfr_init2(dummy2,PRECISION); //dummy variables 
	mpz_t dummy_z; mpz_init(dummy_z);

	//make sure d,q are zeroed out
	for(i=0;i<poly_len;i++){
		mpz_set_ui(d[i],0);
		mpz_set_ui(q[i],0);
	}

	if(verbosity){
	printf("Finding a factor of:\n");
	print_poly(poly_len,poly);}

	//find a root
	mpc_set_d_d(input,0.13,-1.023,MPC_RNDNN); //starting value for rootfind
	//while rootfind doesn't succeed from specified starting point, pick a new starting point
	while(!rootfind_cx(poly,poly_len,input,output,log10thresh,PRECISION)&&iter<iter_max){//find root
		gmp_randseed_ui(seed,rand()); //set random seed
		mpfr_urandomb(dummy,seed); //real part of starting value in [0,1)
		mpfr_urandomb(dummy2,seed); //im part of starting value in [0,1)
		mpfr_mul_si(dummy,dummy,(2*(rand()%2)-1)*2,MPFR_RNDN); //force them to be in (-2,2)
		mpfr_mul_si(dummy2,dummy2,(2*(rand()%2)-1)*2,MPFR_RNDN);
		mpc_set_fr_fr(input,dummy,dummy2,MPC_RNDNN); //starting value for rootfind
		iter++;
	}
	if(iter==iter_max){//failed to find a root after 50 tries
		if(verbosity){
		printf("Failed to find a root.\n");}
		//clear variables and exit
		gmp_randclear(seed);
		mpc_clear(input);
		mpc_clear(output);
		mpfr_clear(output_r);
		mpfr_clear(output_i);
		mpf_clear(delta);
		mpfr_clear(thresh);
		mpfr_clear(dummy);
		mpfr_clear(dummy2);
		mpz_clear(dummy_z);
		return 0;
	}

	if(verbosity){
	printf("root chosen: ");
	mpc_out_str(stdout,10,0,output,MPC_RNDNN);}

	//get real and im parts
	mpc_real(output_r,output,MPC_RNDNN);
	mpc_imag(output_i,output,MPC_RNDNN);

	//check if its an integer. If it is, we're done
	mpfr_set_ui(thresh,10,MPFR_RNDN);
	mpfr_pow_ui(thresh,thresh,log10thresh,MPFR_RNDN);
	mpfr_ui_div(thresh,1,thresh,MPFR_RNDN); //set value of thresh from log10thresh

	mpfr_round(dummy,output_r); //round real and imaginary parts to integer
	mpfr_round(dummy2,output_i); 
	mpfr_sub(dummy,dummy,output_r,MPFR_RNDN); 
	mpfr_abs(dummy,dummy,MPFR_RNDN);//dummy = |difference between output_r and nearest integer|
	mpfr_abs(dummy2,dummy2,MPFR_RNDN);//dummy2 = |imaginary part|

	if((mpfr_cmp(dummy,thresh)<=0)&&(mpfr_cmp(dummy2,thresh)<=0)){//integer check
		if(verbosity){
		printf(" (integer)");}
		mpfr_get_z(dummy_z,output_r,MPFR_RNDN); //round output_r to mpz
		mpz_mul_si(d[0],dummy_z,-1); //set d(x) = x-output, since output is integer
		mpz_set_si(d[1],1);
		//check that d is a divisor
		if(polydivide(poly,d,q,poly_len)==0){
			if(verbosity){
			printf("\nFactor:\n");
			print_poly(poly_len,d);
			printf("Quotient:\n");
			print_poly(poly_len,q);
			printf("\n");}
			//clear variables
			gmp_randclear(seed);
			mpc_clear(input);
			mpc_clear(output);
			mpfr_clear(output_r);
			mpfr_clear(output_i);
			mpf_clear(delta);
			mpfr_clear(thresh);
			mpfr_clear(dummy);
			mpfr_clear(dummy2);
			mpz_clear(dummy_z);
			return 1;
		}
	}
	if(verbosity){
	printf("\n");}


	//otherwise run LLL on each degree less than input degree to find minimal polynomial
	mpz_t *basis;//initialize basis
	basis=malloc((input_degree+1)*(input_degree+3)*sizeof(mpz_t));
	for(i=0;i<(input_degree+1)*(input_degree+3);i++)
		mpz_init(basis[i]);

	for(deg=2;deg<=input_degree;deg++){//loop on degrees
		if(verbosity){
		printf("      LLL on degree %d\n",deg);}
		sig_digits=sig_mpc(output,deg,PRECISION);
		//find irreducible polynomial for chosen root
		create_basis_cx(basis,output,deg,sig_digits,PRECISION);
		LLL(deg+3,deg+1,basis,delta,PRECISION); //use passed PRECISION value
		//LLL(deg+3,deg+1,basis,delta,(3*PRECISION)/4); //use fraction of passed PRECISION value
		for(j=0;j<deg+1;j++)
			mpz_set(d[j],basis[j]); //set first vector of reduced basis to divisor d
		
		//TODO: pick shortest vector instead of first one

		//synthetic division to check it actually divides and to find both factors
		if((monic_slide(deg+1,d)>=0)&&(polydivide(poly,d,q,poly_len)==0)){
			LLL_found_divisor=1;
			if(verbosity){
			printf("Factor:\n");
			print_poly(poly_len,d);
			printf("Quotient:\n");
			print_poly(poly_len,q);}
			break; //quit once you've found lowest degree divisor
		}
	}

	if(LLL_found_divisor==0){
		//no divisor found by LLL, clear variables and exit
		if(verbosity){
		printf("No factor found, increase precision.\n");}
		gmp_randclear(seed);
		mpc_clear(input);
		mpc_clear(output);
		mpfr_clear(output_r);
		mpfr_clear(output_i);
		mpf_clear(delta);
		mpfr_clear(thresh);
		mpfr_clear(dummy);
		mpfr_clear(dummy2);
		mpz_clear(dummy_z);
		for(i=0;i<(input_degree+1)*(input_degree+3);i++)
			mpz_clear(basis[i]);
		free(basis);
		return 0;
	}
	
	//clear variables
	gmp_randclear(seed);
	mpc_clear(input);
	mpc_clear(output);
	mpfr_clear(output_r);
	mpfr_clear(output_i);
	mpf_clear(delta);
	mpfr_clear(thresh);
	mpfr_clear(dummy);
	mpfr_clear(dummy2);
	mpz_clear(dummy_z);
	for(i=0;i<(input_degree+1)*(input_degree+3);i++)
		mpz_clear(basis[i]);
	free(basis);
	if(verbosity){printf("\n");}
	return 1;
}

//factorize poly, store factors in mpz_t factors.
//returns the number of factors (0 if failed, 1 if irreducible, etc)
//not guaranteed to work if poly has factors of higher multiplicity (due 
//to Halley's method rounding). Consequently, one should pass poly/gcd(poly,poly')
//and keep track of the gcd separately
int factorize(mpz_t *poly,int poly_len,int PRECISION,mpz_t *factors,int verbosity, double delta){
	int i;
	int is_reducible=1;
	int factor_counter=0;
	int degree_q=poly_len-1;
	mpz_t *d;
	mpz_t *q;

	//check if poly is monic (must be)
	int degree_poly=degree(poly,poly_len);
	if(mpz_cmp_ui(poly[degree_poly],1)!=0){
		printf("Polynomial not monic, cannot divide\n");
		return 0;
	}

	d=malloc(poly_len*sizeof(mpz_t));//divisor for intermediate step
	q=malloc(poly_len*sizeof(mpz_t));//quotient for intermediate step
	for(i=0;i<poly_len;i++){
		mpz_init(d[i]);
		mpz_init(q[i]);
		mpz_set_ui(d[i],0);
		mpz_set_ui(q[i],0);
	}

	if(verbosity){
	printf("-----------------------------------------------------------------------------\n\n");}
	while(is_reducible&&degree_q>0){
		//find a factor
		is_reducible=find_factor_cx(poly,d,q,degree_q+1,PRECISION,verbosity,delta);
		degree_q=degree(q,poly_len);
		//copy factor d to factor bank
		if(is_reducible){
			for(i=0;i<poly_len;i++)
				mpz_set(factors[factor_counter*poly_len+i],d[i]);
			factor_counter++;
		}
		else{
			for(i=0;i<poly_len;i++){
				mpz_clear(d[i]);
				mpz_clear(q[i]);
			}
			free(d);
			free(q);
			return 0;
		}
		//zero out d (not sure why this is necessary, since it is done in find_factor_cx)
		for(i=0;i<poly_len;i++){mpz_set_ui(d[i],0);}
		//set poly to quotient by factor
		for(i=0;i<poly_len;i++)
			mpz_set(poly[i],q[i]);
	}
	if(verbosity){
	printf("-----------------------------------------------------------------------------\n");}

	//clear variables
	for(i=0;i<poly_len;i++){
		mpz_clear(d[i]);
		mpz_clear(q[i]);
	}
	free(d);
	free(q);
	return factor_counter;
}

//factorizes poly even if it has repeated factors 
//essentially calls factorize() on poly/gcd(poly,poly') and gcd(poly,poly')
//iteratively, to mitigate factoring polynomials with repeated factors
int factorize_full(mpz_t *poly,int poly_len,int PRECISION,mpz_t *factors,int verbosity, double delta){
	//idea: strip off poly/gcd, call factorize(poly/gcd), set poly <- gcd. Repeat while gcd!=1
	int i,total_factors=0,exp_count=0;

	//if it is degree 1 or less: we are done
	if(degree(poly,poly_len)<=1){
		for(i=0;i<poly_len;i++)
			mpz_set(factors[i],poly[i]);
		return 1;
	}

	mpz_t *p=malloc(poly_len*sizeof(mpz_t));
	mpz_t *pp=malloc(poly_len*sizeof(mpz_t));
	mpz_t *gcd_p=malloc(poly_len*sizeof(mpz_t));
	mpz_t *stripped=malloc(poly_len*sizeof(mpz_t));

	for(i=0;i<poly_len;i++){
		mpz_init(p[i]);
		mpz_init(pp[i]);
		mpz_init(stripped[i]);
		mpz_init(gcd_p[i]);
		mpz_set(p[i],poly[i]);
	}
		
	derivative(p,pp,poly_len);
	gcd(p,pp,gcd_p,poly_len);

	while(degree(gcd_p,poly_len)>1){//while there are repeated factors
		derivative(p,pp,poly_len);
		gcd(p,pp,gcd_p,poly_len);		
		if(polydivide(p,gcd_p,stripped,poly_len)!=0){//strip off pice with no repeated factors 
			printf("gcd wasn't a divisor!\n");
			for(i=0;i<poly_len;i++){
				mpz_clear(p[i]);
				mpz_clear(pp[i]);
				mpz_clear(gcd_p[i]);
				mpz_clear(stripped[i]);
			}
			free(p);
			free(pp);
			free(gcd_p);
			free(stripped);

			return 0;
		}
		exp_count++;
		if(verbosity){printf("Factoring exponent-free part #%d\n",exp_count);};
		total_factors=total_factors+factorize(stripped,poly_len,PRECISION,&factors[total_factors*poly_len],verbosity,delta);
		for(i=0;i<poly_len;i++) //set p <- gcd and repeat
			mpz_set(p[i],gcd_p[i]);
	}
	
	if(degree(p,poly_len)>0){//factor remaining part
		total_factors=total_factors+factorize(p,poly_len,PRECISION,&factors[total_factors*poly_len],verbosity,delta);
	}


	//clear variables
	for(i=0;i<poly_len;i++){
		mpz_clear(p[i]);
		mpz_clear(pp[i]);
		mpz_clear(gcd_p[i]);
		mpz_clear(stripped[i]);
	}
	free(p);
	free(pp);
	free(gcd_p);
	free(stripped);

	return total_factors;
}

//make sure p is monic and divide out by highest power of x dividing p
//if the leading coefficient is -1, multiply through by -1 and proceed
//return highest power of x dividing p, -1 if not monic
int monic_slide(int len, mpz_t *p){
	int first_nz=0;
	int last_nz=len-1;
	int i;

	for(i=0;i<len;i++){
		if(mpz_sgn(p[i])!=0){
			first_nz=i;
			break;
		}
	}
	for(i=len-1;i>=0;i--){
		if(mpz_sgn(p[i])!=0){
			last_nz=i;
			break;
		}
	}
	//divide out by highest power of x dividing p
	for(i=first_nz;i<len;i++)
		mpz_set(p[i-first_nz],p[i]);
	for(i=last_nz-first_nz+1;i<len;i++)
		mpz_set_ui(p[i],0);
	
	//monic check - if the leading term is -1, multiply through by -1 to make monic
	//				if not unit, return 0
	if(mpz_cmp_si(p[last_nz-first_nz],-1)==0){
		for(i=0;i<last_nz-first_nz+1;i++)
			mpz_mul_si(p[i],p[i],-1);
	}
	else if(mpz_cmp_si(p[last_nz-first_nz],1)!=0) 
		return -1;

	return first_nz;
}

//same as above, except doesn't multiply through by -1 if that is the leading coef
int monic_slide_dont_multiply(int len, mpz_t *p){
	int first_nz=0;
	int last_nz=len-1;
	int i;

	for(i=0;i<len;i++){
		if(mpz_sgn(p[i])!=0){
			first_nz=i;
			break;
		}
	}
	for(i=len-1;i>=0;i--){
		if(mpz_sgn(p[i])!=0){
			last_nz=i;
			break;
		}
	}
	//divide out by highest power of x dividing p
	for(i=first_nz;i<len;i++)
		mpz_set(p[i-first_nz],p[i]);
	for(i=last_nz-first_nz+1;i<len;i++)
		mpz_set_ui(p[i],0);
	
	if(mpz_cmp_si(p[last_nz-first_nz],1)!=0) 
		return -1;

	return first_nz;
}

//compute the derivative of p, save as pp
void derivative(mpz_t *p,mpz_t *pp,int poly_len){
	int i;
	mpz_set_ui(pp[poly_len-1],0);
	for(i=0;i<(poly_len-1);i++)
		mpz_mul_ui(pp[i],p[i+1],i+1);
}

