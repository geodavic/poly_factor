#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "lll_gs.h" //homemade libraries
#include "poly_functions.h" //for polynomial print

#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define STR_MAX 8192 //maximum number of digits of alpha

//NOT WORKING

int parameter_set(int argc, char *argv[],int *PRECISION, int *verbosity,int *timer,double *delta,int *deg);


int main(int argc, char *argv[]){
    int deg,i,input_fail;
    int PRECISION;//,ACC;
    int sig_digits=0;
    char inputStr[STR_MAX];
    int verbose=1;
    double delta_=0.5;
    int timer=0;

    mpf_t alpha; //input float
    mpf_t delta; //LLL parameter

    //read command line parameters
    if(parameter_set(argc,argv,&PRECISION,&verbose,&timer,&delta_,&deg)==0)
        return 0;


    //check if its an integer
    int past_dot=0;
    for(i=0;inputStr[i]!='\0';i++){
        if(inputStr[i]=='.'){
            past_dot=1;
            break;
        }
    }
    if(past_dot==0){
        if verbose{
            printf("Integer entered! This is clearly algebraic\n");	
        }
        printf("x - %s",inputStr);
        return 0;
    }

    //get degree of alpha
    printf("Enter degree: ");
    printf("(enter exact degree if known)\n");
    scanf("%d",&deg);
    printf("\n");
    if(deg<=0){
        printf("Degree must be positive.\n");
        return 0;
    }

    //get sig_digits
    sig_digits=sig(inputStr,deg);
    //printf("sig_digits: %d\n",sig_digits);

    //set precision level and ACC
    int precision_given=(int)(sig_digits*log2(10));
    PRECISION=precision_given+3;
    //ACC=MIN(900,9*PRECISION/10);


    //initialize alpha and delta with specified bits of precision
    mpf_init2(alpha,PRECISION); 
    mpf_set_ui(alpha,0);

    mpf_init2(delta,PRECISION);
    mpf_set_d(delta,delta_); //LLL parameter

    //read input string into alpha
    input_fail=mpf_set_str(alpha,inputStr,10); 
    if(input_fail){
        printf("read error! Input must only contain 0-9 and at most one '.' and 'e'\n");
        return 0;
    }

    //print precision given
    if(!input_fail&&verbose){
        printf("precision needed: %d bits (operating with %d bits)",precision_given,PRECISION);
        printf("\n\n");
    }

    //create and print basis matrix
    mpz_t *basis;
    basis=malloc((deg+1)*(deg+2)*sizeof(mpz_t));
    for(i=0;i<(deg+1)*(deg+2);i++){
        mpz_init(basis[i]);	
    }
    create_basis(basis,alpha,deg,sig_digits,PRECISION);
    if(verbose){
        printf("Initial lattice basis:\n");
        print_matrix_i(deg+2,deg+1,10,basis);
    }

    if(verbose){
        printf("---------------------------------------------\n\n");
        printf("searching for minimal polynomial...\n\n");
        printf("---------------------------------------------\n\n");
    }

    //perform LLL
    clock_t start=clock(),diff;
    LLL(deg+2,deg+1,basis,delta,PRECISION);
    diff=clock()-start;
    if(verbose){
        printf("\nReduced lattice basis:\n");
        print_matrix_i(deg+2,deg+1,10,basis);
    }

    int msec_time= diff *1000/CLOCKS_PER_SEC;

    //pick out shortest vector, divide by largest power of x dividing the polynomial it represents
    printf("Likely minimal polynomial:\n");

    int poly_len=deg+1;
    i=0;
    while(mpz_sgn(basis[i])==0){
        i++;
        poly_len--;
    }
    print_poly(poly_len,&basis[i],1);
    if(verbose)
        printf("\nTime: \n%d.%ds\n",msec_time/1000,msec_time%1000);

    //clear mp variables
    mpf_clear(alpha);
    mpf_clear(delta);
    for(i=0;i<(deg+1)*(deg+2);i++)
        mpz_clear(basis[i]);
    free(basis);
    return 1;
}

//parse command line input and set the relevant parameters
int parameter_set(int argc, char *argv[],int *PRECISION, int *verbosity,int *timer,double *delta,int *deg){
    int i;

    //no arguments passed
    if(argc==1){
        printf("Input is a monic polynomial in Z[x], written without spaces (e.g. x^2-x+2)\nFormat: <polynomial> <OPTS>\n        OPTS: -v: verbosity\n              -t: timer\n              -p: precision in bits (e.g. -p 150). Default is 64, minimum of 32.\n              -d: LLL parameter (0.25<d<1). Default is 0.5.\n              -newline: print each factor on a new line.\n");
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
                    printf("Precision indicated not an integer.\n");
                    return 0;
                }
                *PRECISION=strtol(argv[i],&argv[i],10);
                if(*PRECISION<=0){
                    printf("Precision indicated not a positive integer.\n");
                    return 0;
                }
                *PRECISION=MAX(*PRECISION,32); //Assume at least 32 bits of precision
            }
            else if(strcmp(argv[i],"-d")==0){
                i++;
                if(i==argc){
                    printf("Delta parameter indicated not recognized.\n");
                    return 0;
                }
                *delta=strtod(argv[i],&argv[i]);
                if(*delta==0.0){
                    printf("Delta parameter indicated not a double.\n");
                    return 0;
                }
                if(*delta<=0.25||*delta>=1){
                    printf("Delta parameter must be in the range 0.25 < d < 1\n");
                    return 0;
                }

            }
            else{
                printf("Input error. Unrecognized options (ensure polynomial input has no spaces).\n");
                return 0;
            }
        }
    }
    return 1;
}

