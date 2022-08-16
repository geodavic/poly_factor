/*

   George D. Torres, 2017.
   Polynomial factorization using LLL. 

 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <ctype.h>
#include <regex.h>
#include "lll_gs.h" //includes gmp.h, mpfr.h, mpc.h, math.h
#include "lll_functions.h" //function library for polynomials

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

//----major-changes-------------//
//GDT 01.2018
//    05.2018
//    07.2018

//---------about----------------//
//polynomial facorization algorithm: find root of f(x) using Newton's method -> use LLL to find algebraic relation on root -> get irreducible divisor of f(x) -> repeat on quotient

//----------notes---------------//
//not a guaranteed algorithm, since LLL method won't always find minimal polynomial, but is nearly guaranteed. Higher precision -> higher likelihood of finding a factor
//LLL is polynomial time, but method requires high precision arithmetic, so it doesn't scale well to higher degree polynomials (see complexity below)
//compleixity for one factorization is approximately O(n^6*A^3), where n is the degree of the input polynomial and A = n+PRECISION*log(2); for fixed (and reasonable) precision, this is O(n^9).
//there is a trade-off between the parameters PRECISION and delta (the LLL parameter). Lower values of delta make LLL run faster, but be less accurate. Raising PRECISION a certain amount will fix this, but at the cost of slower computations. Raising delta makes LLL slower but allows for lower PRECISION, hence faster gmp computations. Not sure where the sweet spot is yet.


int lll_factor(int argc, char *argv[]);
int read_csv(char *polystr,mpz_t *poly,int poly_len);
int csv_len(char *str);
int parameter_set(int argc, char *argv[],int *PRECISION, int *verbosity,int *timer,int *newline,double *delta,int *poly_len, int *stop_deg);


int main(int argc,char *argv[]){
    return !lll_factor(argc,argv);
}

// factor polynomial from command line input (argv,argc)
// return 0 if failed
int lll_factor(int argc, char *argv[]) {
    srand(1); //initialize random
    int i,j;
    int PRECISION=128; //bits of precision (for floats) - default 64
    int poly_len=0; //length of polynomial 
    int factor_counter; //number of factors 
    int verbosity=0; //verbosity bool
    int newline=0; //newline bool (for printing)
    int timer=0; //timer bool
    int stop_deg=0;// stop degree for LLL
    double delta=0.5;//LLL parameter default
    mpz_t *poly; //polynomial coefficients
    mpz_t *allfactors; //factors list
    int *multiplicities; //multiplicities of factors

    //read command line parameters
    if(parameter_set(argc,argv,&PRECISION,&verbosity,&timer,&newline,&delta,&poly_len,&stop_deg)==0)
        return 0;
    if(stop_deg==0){stop_deg=poly_len;}

    //print parameters
    if(verbosity){printf("Working precision: %d\n",PRECISION);printf("LLL parameter: %lf\n",delta);}


    //allocate input polynomial, its list of factors, and their multiplicities
    poly=malloc(poly_len*sizeof(mpz_t));
    allfactors=malloc(poly_len*(poly_len-1)*sizeof(mpz_t));
    multiplicities=malloc(poly_len*sizeof(int));
    for(i=0;i<poly_len;i++){
        mpz_init(poly[i]);
        mpz_set_ui(poly[i],0);
    }
    for(i=0;i<poly_len*(poly_len-1);i++)
        mpz_init(allfactors[i]);


    //read polynomial and print it
    if(!read_csv(argv[1],poly,poly_len)){
        fprintf(stderr,"Input error. Expect <polynomial in csv> <-v for verbosity>\n");
        for(i=0;i<poly_len;i++){
            mpz_clear(poly[i]);
            for(j=0;j<poly_len-1;j++)
                mpz_clear(allfactors[j*poly_len+i]);
        }
        free(allfactors);
        free(poly);
        free(multiplicities);

        return 0;
    }
    if(verbosity){
        printf("Polynomial input: ");
        print_poly(poly_len,poly,1);
    }
    if(verbosity){
        printf("\n");
        printf("================================================================================================\n");
        printf("Begin factorization\n");
        printf("================================================================================================\n\n");
    }

    //factor it
    clock_t start=clock(),diff;
    factor_counter=factorize_full(poly,poly_len,PRECISION,allfactors,multiplicities,verbosity,delta,stop_deg);	
    diff=clock()-start;
    int msec_time=diff*1000/CLOCKS_PER_SEC;

    if(verbosity){
        printf("================================================================================================\n");
        printf("End\n");
        printf("================================================================================================\n\n");
    }

    //print factors
    if(factor_counter>0){
        if(verbosity){printf("Factorization:\n");}
        print_factors(allfactors,multiplicities,factor_counter,poly_len,0,newline);
        printf("\n");
    }
    else{
        fprintf(stderr,"Factorization failed\n");
        for(i=0;i<poly_len;i++){
            mpz_clear(poly[i]);
            for(j=0;j<poly_len-1;j++)
                mpz_clear(allfactors[j*poly_len+i]);
        }
        free(allfactors);
        free(poly);
        free(multiplicities);

        return 0;
    }

    //print time taken
    if(timer)
        printf("Time: %dms\n",msec_time);


    //clear variables
    for(i=0;i<poly_len;i++){
        mpz_clear(poly[i]);
        for(j=0;j<poly_len-1;j++)
            mpz_clear(allfactors[j*poly_len+i]);
    }
    free(allfactors);
    free(poly);
    free(multiplicities);
    return 1;
}


//read csv to polynomial
//return zero if failed, nonzero value if succeeded. 
int read_csv(char *polystr,mpz_t *poly,int poly_len){
    char *p=polystr;
    int current_deg=0;

    // split on commas
    char *token = strtok(p,",");
    while( token != NULL ){
        if(current_deg>=poly_len){
            fprintf(stderr,"Error: polynomial not allocated enough memory\n");
            return 0;
        }
        if(mpz_set_str(poly[current_deg],token,10)){
            fprintf(stderr,"Invalid integer coefficient: %s\n",token);
            return 0;
        };
        token = strtok(NULL,",");
        current_deg++;
    }
    return 1;	
}

//return length of csv string, counted as the number of entries between commas
int csv_len(char *str){
    char *p=str;
    int len=1;
    while(*p){
        if(!strncmp(p,",",1))
            len++;
        p++;
    }
    return len;
}

//parse command line input and set the relevant parameters
int parameter_set(int argc, char *argv[],int *PRECISION, int *verbosity,int *timer,int *newline,double *delta,int *poly_len, int* stop_deg){
    int i;

    //no arguments passed
    if(argc==1){
        printf("Input is a monic polynomial in Z[x], written without spaces (e.g. x^2-x+2)\nFormat: <polynomial> <OPTS>\n        OPTS: -v: verbosity\n              -t: timer\n              -p: precision in bits (e.g. -p 150). Default is 64, minimum of 32.\n              -d: LLL parameter (0.25<d<1). Default is 0.5.\n              -newline: print each factor on a new line.\n              -stop: Stop degree for LLL algorithm. Default is infinity.\n");
        return 0;
    }
    //get putative polynomial length and set options
    else{
        *poly_len=csv_len(argv[1]);
        for(i=2;i<argc;i++){
            if(strcmp(argv[i],"-v")==0)
                *verbosity=1;
            else if(strcmp(argv[i],"-t")==0)
                *timer=1;
            else if(strcmp(argv[i],"-vt")==0||strcmp(argv[i],"-tv")==0){
                *verbosity=1;
                *timer=1;
            }
            else if(strcmp(argv[i],"-newline")==0)
                *newline=1;
            else if(strcmp(argv[i],"-p")==0){
                i++;
                if(i==argc){
                    fprintf(stderr,"Precision indicated not an integer.\n");
                    return 0;
                }
                *PRECISION=strtol(argv[i],&argv[i],10);
                if(*PRECISION<=0){
                    fprintf(stderr,"Precision indicated not a positive integer.\n");
                    return 0;
                }
                *PRECISION=MAX(*PRECISION,32); //Assume at least 32 bits of precision
            }
            else if(strcmp(argv[i],"-d")==0){
                i++;
                if(i==argc){
                    fprintf(stderr,"Delta parameter indicated not recognized.\n");
                    return 0;
                }
                *delta=strtod(argv[i],&argv[i]);
                if(*delta==0.0){
                    fprintf(stderr,"Delta parameter indicated not a double.\n");
                    return 0;
                }
                if(*delta<=0.25||*delta>=1){
                    fprintf(stderr,"Delta parameter must be in the range 0.25 < d < 1\n");
                    return 0;
                }

            }
            else if(strcmp(argv[i],"-stop")==0){
                i++;
                if(i==argc){
                    fprintf(stderr,"Stop parameter not recognized.\n"); 
                    return 0;
                }
                *stop_deg=strtol(argv[i],&argv[i],10);
            }
            else{
                fprintf(stderr,"Input error. Unrecognized options (ensure polynomial input has no spaces).\n");
                return 0;
            }
        }
    }
    return 1;
}

