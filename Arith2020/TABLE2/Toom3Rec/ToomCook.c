/*********************************************

	Fonctions :

multiplication multiprécision vectorisée :

https://software.intel.com/sites/landingpage/IntrinsicsGuide/#!=undefined


*********************************************/
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <immintrin.h>
#include <time.h>

#include "ccount.h"


#include "fonctions.h"

//pour force inline
#undef force_inline
#define force_inline __attribute__((always_inline))


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

	static unsigned long int nA[tTM3R], nB[tTM3R], res[tTM3R<<1], resMul[tTM3R<<1], tmp[tTM3R<<1];
	uint64_t mini = (uint64_t)-1L, mini1 = (uint64_t)-1L;
 
	unsigned long long int timer=0, timer1=0;
	
	
	
	printf("PARAM_N = %d, tTM3R = %d, tTM3R/(256/WORD) =%d\n",PARAM_N,tTM3R,tTM3R/(256/WORD));
	printf("%d\n",((PARAM_N/WORD) + (PARAM_N%WORD == 0 ? 0 : 1)));
	srand(time(NULL));


	for(int j=0; j<t;j++){
		nA[j] = (((unsigned long int)(rand()+rand())<<32)^(rand()+rand()));
		nB[j] = (((unsigned long int)(rand()+rand())<<32)^(rand()+rand()));
	}
		

	
	afficheVect(nA,"nA",tTM3R);
	afficheVect(nB,"nB",tTM3R);


	/***********************************************/
	printf("\ngf2x_mul :\n----------\n");
	
	
	Toom3RecMult(nA,nB,resMul);
	gf2x_mul(res,nA,t,nB,t);
		
	afficheVect(res,"res",tTM3R<<1);
		
	printf("\nComparaison avec ToomCookMult :\n-----------------------------\n");
	
	afficheVect(resMul,"resMul",tTM3R<<1);
	
	for(int i=0; i<tTM3R<<1;i++) tmp[i] = res[i]^resMul[i];
	
	printf("\n");
	afficheVect(tmp,"cmp",tTM3R<<1);
	
	printf("\n");
	
	for(int i=0; i<tTM3R<<1;i++)
		if(res[i]^resMul[i]) flag++;

	printf("flag = %d ; ",flag);
	
	if(!flag) printf("Victory !!!!!!!!!!!!!!!!!!!\n\n");
	else printf("Too bad !!!!!!!!!!!!!!!!!!!\n\n");
	
	flag=0;
	printf("\n\n");
	

	
	
	
	
	/********************
	//Comparaison gf2x ToomCookMult
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
		Toom3RecMult(nA,nB,resMul);
		
		for(int i=0; i<tTM3R<<1;i++)
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
	printf("\t  /***************************/\n");
	printf("\t / Chronométrage !!!!!!!!!!!!/\n");
	printf("\t/***************************/\n\n");
	
	printf("\t\tPARAM_N = %d\n\t\t Size =%d bits\n\n",PARAM_N,SIZE_N);
	
	printf("\ngf2x_mul vs ToomCookMult\n");
	printf("-------------------------\n");

	
	for(int k=0; k<NSAMPLES;k++){
	
		mini = (uint64_t)-1L, mini1 = (uint64_t)-1L;

		for(int j=0; j<t;j++){
			nA[j] = (((unsigned long int)(rand()+rand())<<32)^(rand()+rand()));
			nB[j] = (((unsigned long int)(rand()+rand())<<32)^(rand()+rand()));
		}
		
		for(int i=0;i<NTEST;i++)
		{
			gf2x_mul(res,nA,t,nB,t);
			Toom3Mult(nA,nB,resMul);
		}
		
		for(int i=0;i<NTEST;i++)
		{
			
			STAMP(START)
			gf2x_mul(res,nA,t,nB,t);
			STAMP(STOP)


			STAMP(START1)
			Toom3RecMult(nA,nB,resMul);
			STAMP(STOP1)
			
			if(mini>STOP-START) mini = STOP-START;

			if(mini1>STOP1-START1) mini1 = STOP1-START1;
			
		}

		timer += mini;

		timer1 += mini1;
	}
	
	printf("timer gf2x_mul      : %llu\n",timer/NSAMPLES);	
	
	printf("timer ToomCookMult  : %llu\n",timer1/NSAMPLES);



fin:
	printf("\n");
	printf("Au revoir et merci !\n\n");

}
