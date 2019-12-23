//REGEX PROTOTYPE FOR CLI INPUT

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

