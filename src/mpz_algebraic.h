//GDT 12.2017

//---------about-----------------//

//function library for mpz_algebraic and all other related files
//essentially LLL + gram schmidt + a few other basic functions

//---------notes-----------------//

//speedup notes: - removed ACC condition on project(), seems to save a bit of time. 
//               - tried to start g.s. at k instead of 0, since starting at 0 does some unnecessary computation. This seems to break it. not sure why. EDIT: fixed.
//				 			 - tried replacing Mpf_round with mpfr_get_z inside LLL for perhaps a small speedup.
//               - looks like gram schmidt is the bottleneck of LLL (~93% of cpu time)

//-----------TODO-----------------//
// re-evaluate mpc/mpf_sig(). Not sure I actually used the correct heuristics.
// make fast version of LLL (or gram schmidt) using low-level gmp functions
// maybe also look into more efficient LLL modification


//arbitrary precision and math libraries
#include <gmp.h>
#include <mpc.h>
#include <mpfr.h>
#include <math.h>

#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define STR_MAX 8192 //Maximum string size for alpha. Pretty generous. Be sure this is the same in all files that reference this one. Only used in sig(), which has been replaced by sig_mpf().


void Mpf_round(mpz_t rop, const mpf_t op, int PRECISION);
void Mpf_round_f(mpf_t rop, const mpf_t op, int PRECISION);
void print_matrix(int dim, int nvec, mpf_t *mat);
void print_matrix_i(int dim, int nvec, int base, mpz_t *mat);
void print_vector(int dim, mpf_t *x);
void print_vector_z(int dim, mpz_t *x);
void create_basis(mpz_t *mat,mpf_t alpha,int deg,int sig_digits, int PRECISION);//
void project(int len, mpf_t *x, mpf_t *y, mpf_t *proj, int PRECISION);
void project_zf(int len, mpz_t *x, mpf_t *y, mpf_t *proj, int PRECISION);
void normalize(mpf_t *x, int len, int PRECISION, int ACC);
void sq_norm(mpf_t *x, int len, mpf_t norm, int PRECISION);
void sq_norm_z(mpz_t *x, int len, mpz_t norm);
int sig(char inputStr[STR_MAX],int deg);
int sig_mpf(mpf_t alpha,int deg,int PRECISION); //
int sig_mpc(mpc_t alpha,int deg,int PRECISION);
void gram_schmidt(int start,int dim, int nvec, mpz_t *basis, mpf_t *obasis, int PRECISION);
void gram_coef(int dim,int k,int j, mpz_t *basis, mpf_t *obasis, mpf_t g_coef, int PRECISION);
void LLL(int dim, int nvec, mpz_t *basis, mpf_t delta, int PRECISION);
int shortest_vec(int dim, int nvec, mpz_t *basis);

//variables used throughout:
//PRECISION: working precision level for all floats in bits.
//ACC: accuracy parameter. any vector of square length less than 10^(-ACC) is considered zero. not sure how it should compare to PRECISION, or if it matters for this application. Only used in normalize() at the moment, which isn't part of LLL.



//nearest integer round - store result as mpz
void Mpf_round(mpz_t rop, const mpf_t op, int PRECISION){
    mpf_t dummy; mpf_init2(dummy,PRECISION);
    int compare_int;
    mpf_floor(dummy,op);
    mpf_sub(dummy,op,dummy);
    compare_int=mpf_cmp_d(dummy,0.5);
    mpf_floor(dummy,op);
    if(compare_int>0){
        mpf_add_ui(dummy,dummy,1);
        mpz_set_f(rop,dummy);
    }
    else
        mpz_set_f(rop,dummy);

    mpf_clear(dummy);
}

//nearest integer round - store result as mpf
void Mpf_round_f(mpf_t rop, const mpf_t op, int PRECISION){
    mpf_t dummy; mpf_init2(dummy,PRECISION);
    int compare_int;
    mpf_floor(dummy,op);
    mpf_sub(dummy,op,dummy);
    compare_int=mpf_cmp_d(dummy,0.5);
    mpf_floor(dummy,op);
    if(compare_int>0){
        mpf_add_ui(dummy,dummy,1);
        mpf_set(rop,dummy);
    }
    else
        mpf_set(rop,dummy);

    mpf_clear(dummy);
}



//print float matrix
void print_matrix(int dim, int nvec, mpf_t *mat){
    int i,j;
    for(i=0;i<nvec;i++){
        for(j=0;j<dim;j++){
            mpf_out_str(stdout,10,6,mat[i*dim+j]);
            printf(", ");
        }
        printf("\n");
    }
    printf("\n");
}

//print integer matrix
void print_matrix_i(int dim, int nvec,int base, mpz_t *mat){
    int i,j;
    for(i=0;i<nvec;i++){
        for(j=0;j<dim;j++){
            mpz_out_str(stdout,base,mat[i*dim+j]);
            printf(", ");
        }
        printf("\n");
    }
    printf("\n");
}

//print float vector
void print_vector(int dim, mpf_t *x){
    int i;
    for(i=0;i<dim;i++){
        mpf_out_str(stdout,10,6,x[i]);
        printf(", ");
    }
    printf("\n");
}

//print integer vector
void print_vector_z(int dim, mpz_t *x){
    int i;
    for(i=0;i<dim;i++){
        mpz_out_str(stdout,10,x[i]);
        printf(", ");
    }
    printf("\n");
}

//create lattice basis for given alpha
void create_basis(mpz_t *mat,mpf_t alpha,int deg,int sig_digits, int PRECISION){
    int i,j;
    mpf_t ten_power; //10^sig_digits
    mpf_t alpha_power; //alpha^i
    mpz_t alpha_round; // alpha^i*ten_power rounded to int

    for(i=0;i<deg+1;i++){
        for(j=0;j<deg+2;j++)
            mpz_set_ui(mat[i*(deg+2)+j],0);
    }

    mpf_init2(alpha_power,PRECISION);
    mpf_set_ui(alpha_power,1);

    mpf_init2(ten_power,PRECISION);
    mpf_set_ui(ten_power,10);
    mpf_pow_ui(ten_power,ten_power,sig_digits);

    mpz_init(alpha_round);
    mpz_set_ui(alpha_round,1);

    for(i=0;i<deg+1;i++){
        mpf_pow_ui(alpha_power,alpha,i);//compute alpha^i
        mpf_mul(alpha_power,alpha_power,ten_power);//compute alpha^i*10^sig_digits
        Mpf_round(alpha_round,alpha_power,PRECISION); //round it
        //mpz_set_f(alpha_round,alpha_power);  //old rounding (floor)
        mpz_set(mat[i*(deg+2)+deg+1],alpha_round);  //set it to mat[i][deg+1]
        mpz_set_ui(mat[i*(deg+2)+i],1); //set mat[i][i] to 1
    }
    mpz_clear(alpha_round);
    mpf_clear(ten_power);
    mpf_clear(alpha_power);
}

//same as above, execpt alpha is complex (adding extra column for imaginary part)
//make sure mat is allocated the extra column, i.e. has dimension (deg+1)x(deg+3)
void create_basis_cx(mpz_t *mat,mpc_t alpha,int deg,int sig_digits, int PRECISION){
    int i,j;
    mpfr_t ten_power; //10^sig_digits
    mpc_t alpha_power; //alpha^i*ten_power
    mpfr_t alpha_r; //real part of alpha^i*ten_power
    mpfr_t alpha_i; //im part of alpha^i*ten_power
    mpz_t alpha_round_r; // alpha^i*ten_power rounded to int (real part)
    mpz_t alpha_round_i; // alpha^i*ten_power rounded to int (im. part)

    for(i=0;i<deg+1;i++){
        for(j=0;j<deg+3;j++)
            mpz_set_ui(mat[i*(deg+3)+j],0);
    }

    mpc_init2(alpha_power,PRECISION);

    mpfr_init2(ten_power,PRECISION);
    mpfr_set_ui(ten_power,10,MPFR_RNDN);
    mpfr_pow_ui(ten_power,ten_power,sig_digits,MPFR_RNDN);

    mpfr_init2(alpha_r,PRECISION);
    mpfr_init2(alpha_i,PRECISION);
    mpz_init(alpha_round_r);
    mpz_init(alpha_round_i);

    for(i=0;i<deg+1;i++){
        mpc_pow_ui(alpha_power,alpha,i,MPC_RNDNN);//compute alpha^i
        mpc_mul_fr(alpha_power,alpha_power,ten_power,MPC_RNDNN);//compute alpha^i*10^sig_digits
        mpc_real(alpha_r,alpha_power,MPC_RNDNN); //get real and imag parts
        mpc_imag(alpha_i,alpha_power,MPC_RNDNN);
        mpfr_get_z(alpha_round_r,alpha_r,MPFR_RNDN); //round them
        mpfr_get_z(alpha_round_i,alpha_i,MPFR_RNDN); //round them
        mpz_set(mat[i*(deg+3)+deg+1],alpha_round_r);  //set them to mat[i][deg+1], mat[i][deg+2]
        mpz_set(mat[i*(deg+3)+deg+2],alpha_round_i); 

        mpz_set_ui(mat[i*(deg+3)+i],1); //set mat[i][i] to 1
    }

    mpc_clear(alpha_power);
    mpfr_clear(ten_power);
    mpfr_clear(alpha_r);
    mpfr_clear(alpha_i);
    mpz_clear(alpha_round_r);
    mpz_clear(alpha_round_i);
}




//perform a projection of x onto y. len is the length of x,y. proj must be initialized beforehand.
void project(int len, mpf_t *x, mpf_t *y, mpf_t *proj, int PRECISION){
    mpf_t d; mpf_init2(d,PRECISION); mpf_set_ui(d,0);
    mpf_t s; mpf_init2(s,PRECISION); mpf_set_ui(s,0);
    mpf_t mu; mpf_init2(mu,PRECISION);
    mpf_t dummy; mpf_init2(dummy,PRECISION);
    int i;

    // compute x \cdot y and y\cdot y = |y|^2
    for(i=0;i<len;i++){
        mpf_mul(dummy,x[i],y[i]);
        mpf_add(d,d,dummy); //d+= x[i]*y[i]

        mpf_mul(dummy,y[i],y[i]);
        mpf_add(s,s,dummy); //s += y[i]*y[i]
    }

    /* not sure if this check is necessary if precision is high. uncomment this if div by 0 errors occur 
       mpf_set_ui(dummy,10);
       mpf_pow_ui(dummy,dummy,ACC); 
       mpf_ui_div(dummy,1,dummy);
       is_zero=mpf_cmp(s,dummy); // testing if |x|^2 < 10^(-ACC), which is considered zero. 

       if(is_zero<0){ //projection is zero
       for(i=0;i<len;i++)
       mpf_set_ui(proj[i],0);
       }
       else*/
    mpf_div(mu,d,s); //mu=d/s;
    for(i=0;i<len;i++)
        mpf_mul(proj[i],mu,y[i]); //proj[i]=mu*y[i];

    mpf_clear(d);
    mpf_clear(s);
    mpf_clear(mu);
    mpf_clear(dummy);
}
//same as above, except x is an integer array
void project_zf(int len, mpz_t *x, mpf_t *y, mpf_t *proj, int PRECISION){
    mpf_t d; mpf_init2(d,PRECISION); mpf_set_ui(d,0);
    mpf_t s; mpf_init2(s,PRECISION); mpf_set_ui(s,0);
    mpf_t mu; mpf_init2(mu,PRECISION);
    mpf_t dummy; mpf_init2(dummy,PRECISION);
    mpf_t dummy2; mpf_init2(dummy2,PRECISION);
    int i;

    // compute x \cdot y and y\cdot y = |y|^2
    for(i=0;i<len;i++){
        mpf_set_z(dummy2,x[i]);	
        mpf_mul(dummy,dummy2,y[i]);
        mpf_add(d,d,dummy); //d+= x[i]*y[i]

        mpf_mul(dummy,y[i],y[i]);
        mpf_add(s,s,dummy); //s += y[i]*y[i]
    }

    /* not sure if this check is necessary if precision is high. uncomment this if div by 0 errors occur 
       mpf_set_ui(dummy,10);
       mpf_pow_ui(dummy,dummy,ACC); 
       mpf_ui_div(dummy,1,dummy);
       is_zero=mpf_cmp(s,dummy); // testing if |x|^2 < 10^(-ACC), which is considered zero. 

       if(is_zero<0){ //projection is zero
       for(i=0;i<len;i++)
       mpf_set_ui(proj[i],0);
       }
       else*/
    mpf_div(mu,d,s); //mu=d/s;
    for(i=0;i<len;i++)
        mpf_mul(proj[i],mu,y[i]); //proj[i]=mu*y[i];

    mpf_clear(d);
    mpf_clear(s);
    mpf_clear(mu);
    mpf_clear(dummy);
    mpf_clear(dummy2);
}



//normalize x to unit length. 
void normalize(mpf_t *x, int len, int PRECISION, int ACC){
    mpf_t s; mpf_init2(s,PRECISION); mpf_set_ui(s,0);
    mpf_t dummy; mpf_init2(dummy,PRECISION);
    int i,is_zero;

    for(i=0;i<len;i++){
        mpf_mul(dummy,x[i],x[i]);
        mpf_add(s,s,dummy); //s += x[i]*x[i]
    }

    mpf_set_ui(dummy,10);
    mpf_pow_ui(dummy,dummy,ACC); 
    mpf_ui_div(dummy,1,dummy);
    is_zero=mpf_cmp(s,dummy); // testing if |x|^2 < 10^(-ACC), which is considered zero. 

    if(is_zero>0){ 
        mpf_sqrt(s,s); //s -> sqrt(s)
        for(i=0;i<len;i++)
            mpf_div(x[i],x[i],s); // x -> x/|x|
    }
    else
        printf("Error: cannot normalize zero vector\n");

    mpf_clear(s);
    mpf_clear(dummy);
}


//return |x|^2
void sq_norm(mpf_t *x, int len, mpf_t norm, int PRECISION){
    mpf_t s; mpf_init2(s,PRECISION); mpf_set_ui(s,0);
    mpf_t dummy; mpf_init2(dummy,PRECISION);
    int i;

    for(i=0;i<len;i++){
        mpf_mul(dummy,x[i],x[i]);
        mpf_add(s,s,dummy); //s += x[i]*x[i]
    }
    mpf_set(norm,s);

    mpf_clear(s);
    mpf_clear(dummy);
}

//same as above, except x[] is mpz_t
void sq_norm_z(mpz_t *x, int len, mpz_t norm){
    mpz_t s; mpz_init(s); mpz_set_ui(s,0);
    mpz_t dummy; mpz_init(dummy);
    int i;

    for(i=0;i<len;i++){
        mpz_mul(dummy,x[i],x[i]);
        mpz_add(s,s,dummy); //s += x[i]*x[i]
    }
    mpz_set(norm,s);

    mpz_clear(s);
    mpz_clear(dummy);
}


//outdated; see function below
int sig(char inputStr[STR_MAX],int deg){
    int MAX_PRECISION=(int)STR_MAX*log2(10); //maximum precision of a number based on input buffer size
    mpf_t alpha; mpf_init2(alpha,MAX_PRECISION);
    mpf_set_str(alpha,inputStr,10);
    int i,sig_digits=0,past_dot=0;
    mpf_abs(alpha,alpha);

    //probably don't need to read alpha into an mpf to test if |alpha|<1, but whatever
    if(mpf_cmp_ui(alpha,1)<0){ //if |alpha|<1
        //add (deg+1)*(-log10(alpha)) to sig_digits
        //-log10(alpha) ~= 1+#of zeros after decimal
        sig_digits=sig_digits+deg+1;
        for(i=0;inputStr[i]!='\0';i++){
            if(past_dot==1&&inputStr[i]!='0')
                break;
            else{
                if(past_dot==1)
                    sig_digits=sig_digits+deg+1;
            }
            if(inputStr[i]=='.')
                past_dot=1;
        }
    }
    past_dot=0;
    for(i=0;inputStr[i]!='\0';i++){
        if(past_dot==1)
            sig_digits++;
        if(inputStr[i]=='.')
            past_dot=1;
    }
    mpf_clear(alpha);
    return sig_digits-1; //subtract a little 
}


//get number of necessary significant digits to accurately represent \{alpha,...,\alpha^deg\}. 
//if |alpha|>1, this returns the number of digits after the decimal in alpha. if |alpha|<1, this returns the number of digits after the decimal in alpha plus (deg+1)*(-log10(alpha)), since alpha^{deg} will require more digits than alpha in this case. 
//PRECISION value passed is that of alpha
//NOTE: not actually sure the |alpha|<1 case is the correct heuristic
//NOTE: it is very important that this function does not return too high a value. This causes the tail ends of the powers of alpha to be in the lattice, which is bad because the ends have error.
int sig_mpf(mpf_t alpha,int deg,int PRECISION){
    int sig_digits=0;
    mpfr_t alpha2;mpfr_init2(alpha2,PRECISION);
    mpfr_set_f(alpha2,alpha,MPFR_RNDN);
    mpfr_abs(alpha2,alpha2,MPFR_RNDN);
    if(mpfr_cmp_ui(alpha2,1)<0){//if |alpha|<1
        mpfr_log10(alpha2,alpha2,MPFR_RNDN); 
        //set sig_digis+= (deg+1)*(|log10(alpha2)|): (rounding toward -infinity)
        sig_digits=sig_digits+(deg+1)*labs(mpfr_get_si(alpha2,MPFR_RNDU));
    }
    sig_digits=sig_digits+(int)(PRECISION*log10(2.0));
    mpfr_clear(alpha2);
    return sig_digits-1; 
}


//same as sig_mpf except with complex alpha.
int sig_mpc(mpc_t alpha,int deg,int PRECISION){
    int sig_digits=0;
    mpfr_t norm;mpfr_init2(norm,PRECISION);
    mpc_abs(norm,alpha,MPC_RNDNN);
    if(mpfr_cmp_ui(norm,1)<0){//if |alpha|<1
        mpfr_log10(norm,norm,MPFR_RNDN);
        sig_digits=sig_digits+(deg+1)*labs(mpfr_get_si(norm,MPFR_RNDU));
    }
    sig_digits=sig_digits+(int)(PRECISION*log10(2.0));
    mpfr_clear(norm);
    return sig_digits-1;
}



//gram schmidt reduce ``basis" (save as ``obasis"). Do not normalize.
//both must be dim by nvec sized arrays
//probably not numerically stable. Householders are better for stability, but not sure if that will affect the LLL stage
void gram_schmidt(int start,int dim, int nvec, mpz_t *basis, mpf_t *obasis, int PRECISION){
    int i,j,k,s;	
    mpf_t *dummy;
    dummy=malloc(dim*sizeof(mpf_t));

    for(i=0;i<dim;i++)
        mpf_init2(dummy[i],PRECISION);

    for(k=start;k<nvec;k++){
        //start obasis[k] out as as basis[k]
        for(s=0;s<dim;s++)
            mpf_set_z(obasis[k*dim+s],basis[k*dim+s]);
        //substract off projections
        for(j=0;j<k;j++){
            project(dim,&obasis[k*dim],&obasis[j*dim],dummy,PRECISION);
            for(s=0;s<dim;s++){
                mpf_sub(obasis[k*dim+s],obasis[k*dim+s],dummy[s]); //obasis[k][s]=obasis[k][s]-dummy[s];
            }
        }
    }
    //F R E E  H I M
    for(i=0;i<dim;i++)
        mpf_clear(dummy[i]);
    free(dummy);
}

//calculate the (k,j) gram coefficient of the basis
void gram_coef(int dim,int k,int j, mpz_t *basis, mpf_t *obasis, mpf_t g_coef, int PRECISION){
    mpf_t mu_n; mpf_init2(mu_n,PRECISION); mpf_set_ui(mu_n,0);
    mpf_t mu_d; mpf_init2(mu_d,PRECISION); mpf_set_ui(mu_d,0);
    mpf_t dummy; mpf_init2(dummy,PRECISION);

    int s;
    for(s=0;s<dim;s++){
        mpf_set_z(dummy,basis[k*dim+s]); //cast basis[k][s] to a float
        mpf_mul(dummy,dummy,obasis[j*dim+s]); //dummy = basis[k][s]*obasis[j][s]
        mpf_add(mu_n,mu_n,dummy); //mu_n=mu_n+basis[k][s]*obasis[j][s];

        mpf_set(dummy,obasis[j*dim+s]);
        mpf_mul(dummy,dummy,obasis[j*dim+s]);
        mpf_add(mu_d,mu_d,dummy); //mu_d=mu_d+obasis[j][s]*obasis[j][s];
    }
    mpf_div(g_coef,mu_n,mu_d);// g_coef = mu_n/mu_d

    mpf_clear(mu_n);
    mpf_clear(mu_d);
    mpf_clear(dummy);
}


//perform LLL reduction on basis
//uses standard version of LLL, no improvements as of yet. 
void LLL(int dim, int nvec, mpz_t *basis, mpf_t delta, int PRECISION){
    int i,j,k,s;
    int compare_int;

    //initialize g.s. basis
    mpf_t *obasis;
    obasis=malloc(dim*nvec*sizeof(mpf_t));
    for(i=0;i<(dim*nvec);i++)
        mpf_init2(obasis[i],PRECISION);

    gram_schmidt(0,dim,nvec,basis,obasis,PRECISION);
    k=1;

    //set dummy variables and intermediate variables
    mpz_t rnd; mpz_init(rnd); mpz_set_ui(rnd,0);
    mpz_t dummyz; mpz_init(dummyz);
    mpf_t dummy; mpf_init2(dummy,PRECISION);
    mpf_t dummy2; mpf_init2(dummy2,PRECISION);//probably could be done with only one dummy var, but whatever. its only memory.
    mpf_t mu; mpf_init2(mu,PRECISION); mpf_set_ui(mu,0);
    mpfr_t mu2; mpfr_init2(mu2,PRECISION);

    //LLL loop
    while(k<nvec){
        for(j=k-1;j>=0;j--){
            gram_coef(dim,k,j,basis,obasis,mu,PRECISION);
            mpf_abs(dummy,mu);
            compare_int=mpf_cmp_d(dummy,0.5); //check if abs(mu)>0.5
            if(compare_int>0){
                mpfr_set_f(mu2,mu,MPFR_RNDN); //set mu2=mu. (they are different data types)
                mpfr_get_z(rnd,mu2,MPFR_RNDN); //round using mpfr (new rounding)
                //Mpf_round(rnd,mu,PRECISION);  //set rnd=round(mu). (old rounding)
                for(s=0;s<dim;s++){
                    //set basis[k][s]=basis[k][s]-rnd*basis[j][s];
                    mpz_mul(dummyz,rnd,basis[j*dim+s]);
                    mpz_sub(basis[k*dim+s],basis[k*dim+s],dummyz);
                }
                //update obasis (starting from basis[k])
                gram_schmidt(k,dim,nvec,basis,obasis,PRECISION);
            }
        }
        //compute sq_norm(obasis[k],dim), save as dummy
        gram_coef(dim,k,k-1,basis,obasis,mu,PRECISION); 
        sq_norm(&obasis[k*dim],dim,dummy,PRECISION);

        //compute (delta- mu*mu)*sq_norm(obasis[k-1],dim), save as dummy2
        sq_norm(&obasis[(k-1)*dim],dim,dummy2,PRECISION);
        mpf_mul(mu,mu,mu);
        mpf_sub(mu,delta,mu);
        mpf_mul(dummy2,dummy2,mu); 

        //lovasz condition
        compare_int=mpf_cmp(dummy,dummy2);
        if( compare_int>=0 ){
            k++;
        }
        else{
            //swap basis[k] and basis[k-1]
            for(s=0;s<dim;s++){
                mpz_swap(basis[k*dim+s],basis[(k-1)*dim+s]);
            }
            //update obasis (starting from basis[k-1])
            gram_schmidt(k-1,dim,nvec,basis,obasis,PRECISION);
            k=MAX(k-1,1);
        }
    }

    for(i=0;i<(dim*nvec);i++)
        mpf_clear(obasis[i]);
    free(obasis);
    mpz_clear(rnd);
    mpz_clear(dummyz);
    mpf_clear(dummy);
    mpf_clear(dummy2);
    mpf_clear(mu);
    mpfr_clear(mu2);
}

//return index of vector in basis with shortest l2 length
int shortest_vec(int dim, int nvec, mpz_t *basis){
    int j,index_shortest=0;
    mpz_t norm; mpz_init(norm);
    mpz_t shortest_norm; mpz_init(shortest_norm);

    sq_norm_z(&basis[0],dim,shortest_norm); //set shortest_norm to norm of basis[0]
    for(j=1;j<nvec;j++){
        sq_norm_z(&basis[j*dim],dim,norm);
        if(mpz_cmp(shortest_norm,norm)>0){
            index_shortest=j;
            mpz_set(shortest_norm,norm);
        }
    }
    mpz_clear(norm);
    mpz_clear(shortest_norm);
    return index_shortest;
}








