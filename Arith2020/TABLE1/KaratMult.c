/*********************************************

	Fonctions :

multiplication multiprécision vectorisée :

- zmm
- schoolbook

ligne de compilation :

gcc -std=c99 -mavx2 Stern700-350.c -o Stern700-350


https://software.intel.com/sites/landingpage/IntrinsicsGuide/#!=undefined


*********************************************/
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <immintrin.h>
//#include <gmp.h>
#include <time.h>

#include "ccount.h"


#include "fonctions.h"

//pour force inline
#undef force_inline
#define force_inline __attribute__((always_inline))

#define SIZE_MAX_TEST 131072 // maximum size for the tests



#define NSAMPLES 50

#define NTEST 1000

unsigned long long int START, STOP, START1,STOP1;


/********************************************************************************
*
* MAIN
*
*********************************************************************************/



int main(int argc, char* argv[]){


	int flag=0, counter=0;
	int PARAM_N, SIZE_N, t;
	 	
	static unsigned long int nA[SIZE_MAX_TEST/WORD], nB[SIZE_MAX_TEST/WORD], res[2*SIZE_MAX_TEST/WORD], resMul[2*SIZE_MAX_TEST/WORD], tmp[2*SIZE_MAX_TEST/WORD];
	uint64_t mini = (uint64_t)-1L, mini1 = (uint64_t)-1L;
 
	unsigned long long int timer=0, timer1=0;
	
	
	
	
	PARAM_N = (argc == 1 ? 16384 : atoi(argv[1]));
	
	if (argc > 1)
	{
	  if ((PARAM_N != 512) && (PARAM_N != 1024) && (PARAM_N != 2048) && (PARAM_N != 4096) && (PARAM_N != 8192)
	  && (PARAM_N != 16384) && (PARAM_N != 32768) && (PARAM_N != 65536) && (PARAM_N != 131072))
	  {
	    printf("Usage : KaratMult size (where size is 2^x with 9 <= x <= 17)\n");
	    exit(0);
	  }
	}
	
	SIZE_N = 256*CEIL_DIVIDE(PARAM_N,256);
	t = SIZE_N/WORD;	

	printf("PARAM_N = %d, t = %d, t/(256/WORD) =%d\n",PARAM_N,t,t/(256/WORD));
	printf("%d\n",((PARAM_N/WORD) + (PARAM_N%WORD == 0 ? 0 : 1)));
	srand(time(NULL));


	for(int j=0; j<t;j++){
		nA[j] = (((unsigned long int)(rand()+rand())<<32)^(rand()+rand()));
		nB[j] = (((unsigned long int)(rand()+rand())<<32)^(rand()+rand()));
	}
		

	
	afficheVect(nA,"nA",t);
	afficheVect(nB,"nB",t);


	/***********************************************/
	printf("\ngf2x_mul :\n----------\n");
	
	
	KaratRecPclmul(nA,nB,resMul,t);
	gf2x_mul(res,nA,t,nB,t);
		
		
	afficheVect(res,"res",t<<1);
	printf("STOP-START = %lld\n\n",(STOP-START));
	
	printf("Comparaison avec KaratRecPclmul :\n---------------------------------\n");
	afficheVect(resMul,"resMul",t<<1);
	
	for(int i=0; i<t<<1;i++) tmp[i] = res[i]^resMul[i];
	
	printf("\n");
	afficheVect(tmp,"cmp",t<<1);
	
	printf("\n");
	
	for(int i=0; i<t<<1;i++)
		if(res[i]^resMul[i]) flag++;

	printf("flag = %d ; ",flag);
	
	if(!flag) printf("Victory !!!!!!!!!!!!!!!!!!!\n\n");
	else printf("Too bad !!!!!!!!!!!!!!!!!!!\n\n");
	
	flag=0;
	printf("\n\n");
	

	
	
	
	
	/********************
	//Comparaison gf2x KaratRecMult
	*********************/
	
	printf("\t  /***********************************/\n");
	printf("\t /   Test sur 1000 jeux de données   /\n");
	printf("\t/***********************************/\n\n");

	
	for(int i=0;i<NTEST;i++)
	{


		for(int j=0; j<t;j++){
			nA[j] = (((unsigned long int)(rand()+rand())<<32)^(rand()+rand()));
			nB[j] = (((unsigned long int)(rand()+rand())<<32)^(rand()+rand()));
		}
		

		gf2x_mul(res,nA,t,nB,t);
		KaratRecPclmul(nA,nB,resMul,t);
		
		for(int i=0; i<t<<1;i++)
			if(res[i]^resMul[i]) flag++;
		flag?counter++,flag=0:counter,flag=0;
	
	}
	if(counter) printf("%d erreurs !\nToo bad !!!!!!!!!!!!!!!!!!!\n\n",counter),counter=0;
	else printf("Victory !!!!!!!!!!!!!!!!!!!\n\n");
	counter=0;


	
	//goto fin;
	

	
	/********************
	//Chronométrage !!!!!!!!!!!!!
	*********************/
	
	//goto fin;
	printf("\t  /***************************/\n");
	printf("\t / Chronométrage !!!!!!!!!!!!/\n");
	printf("\t/***************************/\n\n");
	
	printf("\t\tPARAM_N = %d\n\t\t Size =%d bits\n\n",PARAM_N,SIZE_N);
	

	printf("t = %d, t/(256/WORD) =%d, \n",t,t/(256/WORD));

	
	printf("\ngf2x_mul vs KaratRec\n");
	printf("---------------------------------\n");

	
	for(int k=0; k<NSAMPLES;k++){
	
		mini = (uint64_t)-1L, mini1 = (uint64_t)-1L;

		for(int j=0; j<t;j++){
			nA[j] = (((unsigned long int)(rand()+rand())<<32)^(rand()+rand()));
			nB[j] = (((unsigned long int)(rand()+rand())<<32)^(rand()+rand()));
		}
		
		for(int i=0;i<NTEST;i++)
		{
			gf2x_mul(res,nA,t,nB,t);
			KaratRecPclmul(nA,nB,resMul,t);
		}
		
		for(int i=0;i<NTEST;i++)
		{
			
			STAMP(START)
			gf2x_mul(res,nA,t,nB,t);
			STAMP(STOP)


			STAMP(START1)
			KaratRecPclmul(nA,nB,resMul,t);
			STAMP(STOP1)
			
			if(mini>STOP-START) mini = STOP-START;

			if(mini1>STOP1-START1) mini1 = STOP1-START1;
			
		}

		timer += mini;

		timer1 += mini1;
	}
	
	printf("timer gf2x_mul        : %llu\n",timer/NSAMPLES);	
	
	printf("timer KaratRecPclmul  : %llu\n",timer1/NSAMPLES);



fin:
	printf("\n");//*/
	printf("Au revoir et merci !\n\n");

}
