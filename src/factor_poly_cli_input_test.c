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
#include "mpz_algebraic.h" //includes gmp.h, mpfr.h, mpc.h, math.h
#include "poly_functions.h" //function library for polynomials

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

//----major-changes-------------//
//GDT 01.2018
//    05.2018
//    07.2018

//---------about----------------//

//polynomial facorization algorithm: find root of f(x) using Newton's method -> use LLL to find algebraic relation on root -> get irreducible divisor of f(x) -> repeat on quotient
//assume f(x) monic for now
//compile with 'gcc -Wall -Wextra -o factor_poly factor_poly.c -lgmp -lmpfr -lmpc' (or type make)
//this is the command line input version
//syntax: ./factor_poly <polynomial in csv> <-v for verbose mode> <-t for time display>
//      csv polynomial example: x^2 - x + 3 is encoded as 3,-1,1
//		if both verbosity and time are required, can use -vt or -tv
//there is a python script in msc_code that will parse a polynomial into csv format

//----------notes---------------//

//not a guaranteed algorithm, since LLL method won't always find minimal polynomial, but is nearly guaranteed. Higher precision -> higher likelihood of finding a factor
//LLL is polynomial time, but method requires high precision arithmetic, so it doesn't scale well to higher degree polynomials (see complexity below)
//compleixity for one factorization is approximately O(n^6*A^3), where n is the degree of the input polynomial and A = n+PRECISION*log(2); for fixed (and reasonable) precision, this is O(n^9).
//there is a trade-off between the parameters PRECISION and delta (the LLL parameter). Lower values of delta make LLL run faster, but be less accurate. Raising PRECISION a certain amount will fix this, but at the cost of slower computations. Raising delta makes LLL slower but allows for lower PRECISION, hence faster gmp computations. Not sure where the sweet spot is yet.


//----------TODO----------------//
//replace python script with c script to parse polynomial input

//It seems like calling factorize_full many times leads to segfaults sometimes. This is likely due to memory not being cleared somewhere. Need to find where the buildup is. UPDATE: rand memory leak detection with valgrind and fixed a few malloc bugs. This didn't fix the issue, however. 
//continue messing with valgrind. Suspect possible bug in GMP?

//multi-thread the starting root. Execute same alg for different roots so that roots of low degree are more likely to be found first.


int read_poly_better(char *polystr,char *csvstr);
char *str_replace(char *orig, char *rep, char *with);
int read_poly(char *polystr,mpz_t *poly,int poly_len);
int read_degree(char *polystr);
int read_csv(char *polystr,mpz_t *poly,int poly_len);
int csv_len(char *str);
int parameter_set(int argc, char *argv[],int *PRECISION, int *verbosity,int *timer,double *delta,int *poly_len);


int main(int argc,char *argv[]){

	//variables used
	srand(1); //initialize random
	int i,j;
	int PRECISION=64; //bits of precision (for floats) - default 64
	int poly_len; //length of polynomial 
	int factor_counter; //number of factors 
	int verbosity=0; //verbosity bool
	int timer=0; //timer bool
	double delta=0.3;//LLL parameter default
	mpz_t *poly; //polynomial coefficients
	mpz_t *allfactors; //factors list

	int a;
	char *test_str="x^3+x^2+3x-2";
	char *csvstr=NULL;
	printf("%s\n",test_str);
	a=read_poly_better(test_str,csvstr);
	printf("%d\n",a);

	//read command line parameters
	if(parameter_set(argc,argv,&PRECISION,&verbosity,&timer,&delta,&poly_len)==0)
		return 0;
	
	//print parameters
	if(verbosity){printf("Working precision: %d\n",PRECISION);printf("LLL parameter: %lf\n",delta);}


	//allocate input polynomial and its list of factors
	poly=malloc(poly_len*sizeof(mpz_t));
	allfactors=malloc(poly_len*(poly_len-1)*sizeof(mpz_t));
	for(i=0;i<poly_len;i++){
		mpz_init(poly[i]);
		mpz_set_ui(poly[i],0);
	}
	for(i=0;i<poly_len*(poly_len-1);i++)
		mpz_init(allfactors[i]);
	
	
	//read polynomial and print it
	if(!read_csv(argv[1],poly,poly_len)){
		printf("Input error. Expect <polynomial in csv> <-v for verbosity>\n");
		return 0;
	}
	printf("Polynomial input:\n");
	print_poly(poly_len,poly,1);
	if(verbosity){
		printf("\n");
	};


	//check that it is monic and divide out by highest power of x dividing it
	int trivial_power=monic_slide_dont_multiply(poly_len,poly);
	if(trivial_power<0){
		printf("Polynomial not monic. Unable to divide.\n");
		return 0;
	}


	//factor it
	clock_t start=clock(),diff;
	factor_counter=factorize_full(poly,poly_len,PRECISION,allfactors,verbosity,delta);	
	diff=clock()-start;
	int msec_time=diff*1000/CLOCKS_PER_SEC;


	//print factors
	if(factor_counter>0){printf("Factorization:\n");
	if(trivial_power>0)
		printf("x^%d\n",trivial_power);
	for(j=0;j<factor_counter;j++)
		print_poly(poly_len,&allfactors[j*poly_len],1);
	}


	//print time taken
	if(timer)
		printf("\nTime: %dms\n",msec_time);


	//clear variables
	for(i=0;i<poly_len;i++){
		mpz_clear(poly[i]);
		for(j=0;j<poly_len-1;j++)
			mpz_clear(allfactors[j*poly_len+i]);
	}
	free(allfactors);
	free(poly);
	return 0;
}

// credit: jmucchiello, StackOverflow
char *str_replace(char *orig, char *rep, char *with) {
    char *result; // the return string
    char *ins;    // the next insert point
    char *tmp;    // varies
    int len_rep;  // length of rep (the string to remove)
    int len_with; // length of with (the string to replace rep with)
    int len_front; // distance between rep and end of last rep
    int count;    // number of replacements

    // sanity checks and initialization
    if (!orig || !rep)
        return NULL;
    len_rep = strlen(rep);
    if (len_rep == 0)
        return NULL; // empty rep causes infinite loop during count
    if (!with)
        with = "";
    len_with = strlen(with);

    // count the number of replacements needed
    ins = orig;
    for (count = 0; (tmp = strstr(ins, rep)); ++count) {
        ins = tmp + len_rep;
    }

    tmp = result = malloc(strlen(orig) + (len_with - len_rep) * count + 1);

    if (!result)
        return NULL;

    // first time through the loop, all the variable are set correctly
    // from here on,
    //    tmp points to the end of the result string
    //    ins points to the next occurrence of rep in orig
    //    orig points to the remainder of orig after "end of rep"
    while (count--) {
        ins = strstr(orig, rep);
        len_front = ins - orig;
        tmp = strncpy(tmp, orig, len_front) + len_front;
        tmp = strcpy(tmp, with) + len_with;
        orig += len_front + len_rep; // move to next "end of rep"
    }
    strcpy(tmp, orig);
    return result;
}

//read polynomial to csv
//idea: use regex to extract monomials
//returns poly_len
//Not working: bus error
int read_poly_better(char *polystr,char *csvstr){
	//int deg=0;
	int term_count=0;
	int status;
	regex_t re;
	regmatch_t term[5];
	char *p=polystr;
	char *pattern="([-])*\\s*(\\d*)(x?)([\\^\\d]*)";

	//regularize terms
	p=str_replace(p,"-","+-");

	//count number of terms
	while(*p){
		if(*p=='+'){
			term_count++;	
		}
		p++;
	}
	term_count++;

	//use regex to find terms
	p = polystr;
	if (regcomp(&re,pattern,REG_EXTENDED) != 0 ){
		//error in compiling regex
		return(0);
	}
	printf("term_count=%d\n",term_count);	
	int start,stop;
	status=regexec(&re,p,term_count+2,term,0);
	int counter=0;
	while(status==0 && counter<5){
		start=(int)(term[0].rm_so);
		stop=(int)(term[0].rm_eo);
		printf("start=%d,stop=%d\n",start,stop);
		//matched term as string
		char termbuff[stop-start+1];
		memcpy(termbuff,&p[start],stop-start);
		termbuff[stop-start]='\0';

		printf("%s\n",termbuff);

		//find next match
		p=p+stop;
		printf("%s\n",p);
		status=regexec(&re,p,term_count+2,term,0);
		counter++;
	}
	
	regfree(&re);

	return(1);
}

//read polynomial from polystr; return length (=degree+1), length 0 if failed
//can only process coefficients of type long for now
//input format requires that the coefficients be two characters before the power of x (e.g. 10x^2)
//made poorly, not working
//int read_poly(char *polystr,char *csvstr,int poly_len){
int read_poly(char *polystr,mpz_t *poly,int poly_len){
	char *p=polystr;
	int deg=0;
	int is_neg;
	int max_deg=0;

	//polystr=replace_char(polystr)

	signed long coef;
	while(*p){
		if(!strncmp(p,"-",1))
			is_neg=-1;
		if(isdigit(*p)){
			coef=strtol(p,&p,10);
			if(isdigit(*(p+2))&&!strncmp(p+1,"^",1)){//terms of degree 2 or more
				p=p+2;
				deg=strtol(p,&p,10); //degree is two characters after coefficient
				if(deg>max_deg)
					max_deg=deg;
				if(deg<poly_len){
					mpz_set_si(poly[deg],is_neg*coef);
					is_neg=1;
				}
				else{
					printf("Error: polynomial not allocated enough memory\n");
					return 0;
				}
			}
			else if(!strncmp(p," ",1)||!strncmp(p,"+",1)||!strncmp(p,"-",1)||!strncmp(p,"\0",1)){//constant term
				mpz_set_si(poly[0],is_neg*coef);
				is_neg=1;
			}
			else if(!strncmp(p,"x",1)&&(!strncmp(p+1," ",1)||!strncmp(p+1,"+",1)||!strncmp(p+1,"-",1)||!strncmp(p+1,"\0",1))){//degree 1 term
				p++;
				mpz_set_si(poly[1],is_neg*coef);
				is_neg=1;
			}
			else
				return 0; //improper format
			
		}//next come coefficient \pm 1 terms
		else if(!strncmp(p,"x",1)&&(!strncmp(p+1," ",1)||!strncmp(p+1,"+",1)||!strncmp(p+1,"-",1)||!strncmp(p+1,"\0",1))){//degree 1 term
			mpz_set_si(poly[1],is_neg);
			is_neg=1;
			p++;
		}
		else if(!strncmp(p,"x",1)&&(!strncmp(p+1,"^",1))){//degree 2 or more terms
			if(isdigit(*(p+2))){
				p=p+2;
				deg=strtol(p,&p,10);
				if(deg>max_deg)
					max_deg=deg;
				if(deg<poly_len){
					mpz_set_si(poly[deg],is_neg*coef);
					is_neg=1;
				}
				else{
					printf("Error: polynomial not allocated enough memory\n");
					return 0;
				}
			}
			else
				return 0;//improper format
		}
		else
			p++;
	}
	return max_deg+1;
}

//read degree of polynomial string, return 0 if failed
int read_degree(char *polystr){
	char *p=polystr;
	int max_deg=0;
	int deg;
	while(*p){
		if(!strncmp(p,"^",1)){
			if(isdigit(*(p+1))){
				p++;
				deg=strtol(p,&p,10);
				if(deg>max_deg)
					max_deg=deg;
			}
			else
				return 0;
		}
		else if(!strncmp(p,"x",1)){
			if(1>max_deg)
				max_deg=1;
			p++;
		}
		else
			p++;
	}
	return max_deg;
}

//read csv to polynomial
//return zero if failed, nonzero value if succeeded. 
int read_csv(char *polystr,mpz_t *poly,int poly_len){
	char *p=polystr;
	int is_neg=1; //keep track of negative 
	int current_deg=0;
	signed long coef;
	while(*p){
		if(!strncmp(p,",",1))
			p++;
		else if(!strncmp(p,"-",1)){
			is_neg=-1;
			p++;
		}
		else if(isdigit(*p)){
			coef=strtol(p,&p,10);
			if(current_deg>=poly_len){
				printf("Error: polynomial not allocated enough memory\n");
				return 0;
			}
			else{
				mpz_set_si(poly[current_deg],is_neg*coef);
				is_neg=1;
				current_deg++;
			}
		}
		else
			return 0; //not a number or sign or comma (improper format)
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
int parameter_set(int argc, char *argv[],int *PRECISION, int *verbosity,int *timer,double *delta,int *poly_len){
	int i;

	//no arguments passed
	if(argc==1){
		printf("Input is a monic polynomial in Z[x], written without spaces (e.g. x^2-x+2)\nFormat: <polynomial> <OPTS>\n        OPTS: -v: verbosity\n              -t: timer\n              -p: precision in bits (e.g. -p 150). Default is 64, minimum of 16.\n              -d: LLL parameter (0.25<d<1). Default is 0.3.\n");
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
				*PRECISION=MAX(*PRECISION,16); //Assume at least 16 bits of precision
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
				printf("Input error. Unrecognized options.\n");
				return 0;
			}
		}
	}
	return 1;
}
















