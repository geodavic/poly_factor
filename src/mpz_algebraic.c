#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "lll_gs.h" //homemade libraries
#include "lll_functions.h" //for polynomial print

//GDT 12.2017
//---------about-----------------//

//This determines minimal polynomial of an algebraic integer alpha, given as input. Result is a 'best guess' based on lattice reduction. Not guaranteed to work, but throwing enough precision at it will usually work.
//compile with 'gcc -Wall -Wextra -o mpz_algebraic mpz_algebraic.c -lgmp -lmpfr'
//first attempt with GMP libraries
//First working version: 12/25/17
//Mostly spaghetti code at the moment

//----------notes---------------//

//this can also be thought of as an algorithm that finds an integer polynomial that has a root which is very close to the given number. Of course, there is always a trivial such polynomial: p(x)=xq-p, where alpha = p/q. This finds more interesting (higher degree) algebraic relations on alpha which might be causal (i.e. have very small coefficients)
//there is a subtlety to the degree chosen. if you choose it below deg(alpha), you won't find the integer relation. if you choose it too high, the arithmetic seems to garble the relation. ideally, you input the exact degree of alpha, no more.
//this was able to find the min poly of 3*sqrt(2)+sqrt(3)+sqrt(5)+sqrt(7), which is degree 16, with about 750 bits of alpha in about 25 seconds. Not a simple minimal polynomial (huge coefficients!)
//maybe try householder reflections instead of gram schmidt for LLL. not sure if LLL requires gs or not. 

//speedup ideas:
//be smarter about gs inside lll- only update based on vectors that change. EDIT: made this change, about ~20-30% speed up. Less than I would have thought (was expecting closer to %50);


//------------------------------//

#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define STR_MAX 8192 //maximum number of digits of alpha

int main(int argc, char *argv[]){
    int deg,i,input_fail,past_dot,precision_given;
    int PRECISION;//,ACC;
    int sig_digits=0;
    char inputStr[STR_MAX];
    int verbose=1;

    mpf_t alpha; //input float
    mpf_t delta; //LLL parameter

    if(argc==2){
        if(strcmp(argv[1],"-v")==0)
            verbose=1; //verbosity option 
    }

    //print statements and variable reading
    printf("_____________________________________________\n\n");
    printf("|     Check that a number is algebraic      |\n");
    printf("| -Now with arbitrary precision             |\n");
    printf("|                                           |\n");
    printf("|                         GDT December 2017 |\n");
    printf("---------------------------------------------\n\n");
    printf("Enter number: ");
    printf("(more digits preferable)\n");
    scanf("%s",inputStr);

    //check if its an integer
    past_dot=0;
    for(i=0;inputStr[i]!='\0';i++){
        if(inputStr[i]=='.'){
            past_dot=1;
            break;
        }
    }
    if(past_dot==0){
        printf("Integer entered! This is clearly algebraic\n");	
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
    precision_given=(int)(sig_digits*log2(10));
    PRECISION=precision_given+3;
    //ACC=MIN(900,9*PRECISION/10);


    //initialize alpha and delta with specified bits of precision
    mpf_init2(alpha,PRECISION); 
    mpf_set_ui(alpha,0);

    mpf_init2(delta,PRECISION);
    mpf_set_d(delta,0.75); //LLL parameter

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


