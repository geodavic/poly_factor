//Function library for all polynomial factorization files (these files must also include lll_gs.h, which this references extensively)
//see lll_factor.c for more detailed documentation

//----major changes------------//
//GDT 01.2018
//    05.2018
//    07.2018

//------------TODO------------------//
//allow for non-monic polynomials
//implement irreducibility tests
//			- Eisenstein (plus shifts?)
//      - Perron's criterion
//      - Murdy, Osada, Bauer criteria
//      - reverting
//      - Newton polygons?


void print_poly(int len,const  mpz_t *x,int newline);
void print_factors(mpz_t *factors,int *multiplicities, int num_factors, int poly_len,int trivial_power,int newline);
void evaluate_cx(mpz_t *p, int len, const mpc_t input, mpc_t output, int PRECISION);
int degree(mpz_t *p, int len);
int degree_q(mpq_t *p, int len);
int rootfind_cx(mpz_t *p, int len, mpc_t start,mpc_t root,int log10_thresh, int PRECISION);
int polydivide(mpz_t *p,mpz_t *d,mpz_t *out,int len);
int polydivide_r(mpq_t *p,mpq_t *d,mpq_t *r,int len);
void gcd(mpz_t *poly1, mpz_t *poly2, mpz_t *gcd, int poly_len);
int find_factor_cx(mpz_t *poly, mpz_t *d, mpz_t *q, int poly_len,int PRECISION,int verbosity, double d_delta, int stop_deg);
int factorize(mpz_t *poly,int poly_len,int PRECISION,mpz_t *factors,int verbosity, double delta, int stop_deg);
int factorize_full(mpz_t *poly,int poly_len,int PRECISION,mpz_t *factors, int *multiplicities,int verbosity, double delta, int stop_deg);
int monic_slide(int len, mpz_t *p);
int monic_slide_dont_multiply(int len, mpz_t *p);
void derivative(mpz_t *p,mpz_t *pp,int poly_len);

//print polynomial with coefficient list x. (pass newline=1 if \n is needed)
void print_poly(int len,const  mpz_t *x,int newline){
    int i;
    int sgn;
    int first_done=0;
    mpz_t abs; mpz_init(abs);
    if(mpz_sgn(x[0])!=0&&len>0){
        mpz_out_str(stdout,10,x[0]);
        first_done=1;
    }
    for(i=1;i<len;i++){
        sgn=mpz_sgn(x[i]);
        if(sgn > 0){
            if(mpz_cmp_ui(x[i],1)==0){
                if(first_done){
                    if(i!=1)
                        printf(" + x^%d",i);
                    else
                        printf(" + x");
                }
                else{
                    if(i!=1)
                        printf("x^%d",i);
                    else
                        printf("x");
                    first_done=1;
                }
            }
            else{
                if(first_done){
                    printf(" + ");
                    mpz_out_str(stdout,10,x[i]);
                    if(i!=1)
                        printf("x^%d",i);
                    else
                        printf("x");
                }
                else{
                    mpz_out_str(stdout,10,x[i]);
                    if(i!=1)
                        printf("x^%d",i);
                    else
                        printf("x");
                    first_done=1;
                }
            }
        }
        if(sgn < 0){
            if(mpz_cmp_d(x[i],-1.0)==0){
                if(first_done){
                    if(i!=1)
                        printf(" - x^%d",i);
                    else
                        printf(" - x");
                }
                else{
                    if(i!=1)
                        printf("-x^%d",i);
                    else
                        printf("-x");
                    first_done=1;
                }

            }
            else{
                if(first_done){
                    printf(" - ");
                    mpz_neg(abs,x[i]);
                    mpz_out_str(stdout,10,abs);
                    if(i!=1)
                        printf("x^%d",i);
                    else
                        printf("x");
                }
                else{
                    mpz_out_str(stdout,10,x[i]);
                    if(i!=1)
                        printf("x^%d",i);
                    else
                        printf("x");
                    first_done=1;
                }
            }
        }
    }
    if(newline)
        printf("\n");
    mpz_clear(abs);
}

//print list of factors of a polynomial
//trivial_power is the highest power of x dividing the polynomial
void print_factors(mpz_t *factors,int *multiplicities, int num_factors, int poly_len,int trivial_power,int newline){
    int j;
    if(trivial_power==1)
        printf("x");
    else if(trivial_power>1)
        printf("x^%d",trivial_power);
    if(newline)
        printf("\n");
    if(num_factors>0){
        for(j=0;j<num_factors;j++){
            printf("(");  
            print_poly(poly_len,&factors[j*poly_len],0);
            printf(")");
            if(multiplicities[j]>1)
                printf("^%d",multiplicities[j]);
            if(newline)
                printf("\n");
        } 
    }
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


//find one complex root of p(x) using second order Newton's method (Halley's method). 
//cubic convergence, stop when log10(|xn-x(n+1)|)<-log10_thresh 
//returns 1 on success and 0 on failure (i.e. exceeding a certain amount of iterations without getting within the threshold of zero).
//note: start value should not be totally real, since this iteration sends reals to reals
//    : make sure p has no repeated roots (otherwise this isn't guaranteed to converge)
//TODO: make exit condition depend on abs(Re(x)) and abs(Im(x))
int rootfind_cx(mpz_t *p, int len, mpc_t start,mpc_t root,int log10_thresh, int PRECISION){

    // If p is linear, the root is already known: -p[0]
    if(len<=2){
        mpc_set_z(root,p[0],MPC_RNDNN);
        mpc_ui_sub(root,0,root,MPC_RNDNN); //root = -root
        return 1;
    }

    int i;
    int max_iterates=100;
    mpfr_t diff; mpfr_init2(diff,PRECISION);
    mpc_t quot; mpc_init2(quot,PRECISION);
    mpc_t dummy; mpc_init2(dummy,PRECISION);
    mpfr_t thresh; mpfr_init2(thresh,PRECISION);
    mpz_t *pp; pp=malloc((len)*sizeof(mpz_t)); //valgrind sees an issue here
    mpz_t *ppp; ppp=malloc((len)*sizeof(mpz_t)); //or here
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
        fprintf(stderr,"warning! tried to divide by zero polynomial.\n");
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
//notes: - might be able to reduce down to at most one dummy variable of each data type
int find_factor_cx(mpz_t *poly, mpz_t *d, mpz_t *q, int poly_len,int PRECISION,int verbosity, double d_delta, int stop_deg){
    int i,j,iter=0,iter_max=3;
    int sig_digits,deg,input_degree=poly_len-1;
    int LLL_found_divisor=0;
    int LLL_hit_cap=0;
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
        print_poly(poly_len,poly,1);}

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
        fprintf(stderr,"Failed to find a root.\n");
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
                printf("\nFactor found:\n");
                printf("--> ");
                print_poly(poly_len,d,0);
                printf(" <--\n");
                printf("Quotient:\n");
                print_poly(poly_len,q,1);
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

        // break if cap hit
        if(deg==stop_deg){
            LLL_hit_cap=1;
            break; 
        }

        if(verbosity){
            printf("      LLL searching for factor of degree %d...",deg);}
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
                printf("\nFactor found:\n");
                printf("--> ");
                print_poly(poly_len,d,0);
                printf(" <--\n");
                printf("Quotient:\n");
                print_poly(poly_len,q,1);}
            break; //quit once you've found lowest degree divisor
        }
        if(verbosity){
            printf(" none\n");
        }
    }

    if(!LLL_found_divisor||LLL_hit_cap){
        //no divisor found by LLL, clear variables and exit
        if(LLL_hit_cap){
            fprintf(stderr,"Maximum allowed degree %d hit.\n",stop_deg);
        }
        else{
            fprintf(stderr,"No factor found, increase precision or delta parameter.\n");
        }

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
int factorize(mpz_t *poly,int poly_len,int PRECISION,mpz_t *factors,int verbosity, double delta, int stop_deg){
    int i;
    int is_reducible=1;
    int factor_counter=0;
    int degree_q=poly_len-1;
    mpz_t *d;
    mpz_t *q;

    //check if poly is monic (must be)
    int degree_poly=degree(poly,poly_len);
    if(mpz_cmp_ui(poly[degree_poly],1)!=0){
        fprintf(stderr,"Polynomial not monic, cannot divide\n");
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


    while(is_reducible&&degree_q>0){
        //find a factor
        is_reducible=find_factor_cx(poly,d,q,degree_q+1,PRECISION,verbosity,delta,stop_deg);
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
//If p = f_1^{n_1} * ... * f_k^{n_k}, then this first finds f_1,...,f_k then n_1,...,n_k
//the f_i are stored in factors and the n_i are stored in multiplicities
//finds the largest degree, square free factor. Factors that and then find the multiplicities of those factors
int factorize_full(mpz_t *poly,int poly_len,int PRECISION,mpz_t *factors, int *multiplicities,int verbosity, double delta, int stop_deg){
    int i,j,mult,new_factors=0;

    //if it is degree 1 or less: we are done
    if(degree(poly,poly_len)<=1){
        for(i=0;i<poly_len;i++)
            mpz_set(factors[i],poly[i]);
        return 1;
    }

    //allocate auxillary polynomials used
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

    //check that it is monic and divide out by highest power of x dividing it
    int trivial_power=monic_slide_dont_multiply(poly_len,p);
    if(trivial_power<0){
        fprintf(stderr,"Polynomial not monic. Unable to divide.\n");
        return 0;
    }
    else if(verbosity && trivial_power>0){
        printf("Trivial factor found:\n");
        printf("x^%d\n\n",trivial_power);
    }

    //compute gcd of p and p'
    derivative(p,pp,poly_len);
    gcd(p,pp,gcd_p,poly_len);

    int gcd_deg=degree(gcd_p,poly_len);
    if(gcd_deg>0){//if higher multiplicity factors exist (p is not square-free)
        if(polydivide(p,gcd_p,stripped,poly_len)!=0){//strip off pice with no repeated factors 
            fprintf(stderr,"gcd wasn't a divisor!\n");
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
        if(verbosity)
            printf("**Higher multiplicity factors detected.**\n\n");

        //set p to square-free part of p
        for(i=0;i<poly_len;i++)
            mpz_set(p[i],stripped[i]);
    }

    if(degree(p,poly_len)>0){//factor square-free part
        new_factors=factorize(p,poly_len,PRECISION,&factors[0],verbosity,delta,stop_deg);
        if(new_factors==0){
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
    }

    //ensure multiplicities are all zero to start
    for(i=0;i<poly_len;i++)
        multiplicities[i]=0;


    //compute multiplicites of each factor
    if(gcd_deg>0){
        if(verbosity){printf("Counting multiplicities:\n");}
        for(j=0;j<new_factors;j++){
            //set p back to original polynomial
            for(i=0;i<poly_len;i++)
                mpz_set(p[i],gcd_p[i]);
            mult=1;

            //set pp to be the jth factor
            for(i=0;i<poly_len;i++)
                mpz_set(pp[i],factors[j*poly_len+i]);

            //divide by pp until you can't
            while(polydivide(p,pp,stripped,poly_len)==0){
                mult++;
                for(i=0;i<poly_len;i++)
                    mpz_set(p[i],stripped[i]);

            }
            multiplicities[j]=mult;
            if(verbosity){
                printf("factor: ");
                print_poly(poly_len,pp,0);
                printf("\nmultiplicity: %d\n",mult);
            }
        }
    }

    //add x^i to list of factors 
    if(trivial_power>0){
        for(i=0;i<poly_len;i++)
            mpz_set_ui(factors[new_factors*poly_len+i],0);
        mpz_set_ui(factors[new_factors*poly_len+trivial_power],1);
        multiplicities[new_factors]=1;
        new_factors++;
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

    return new_factors;
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

