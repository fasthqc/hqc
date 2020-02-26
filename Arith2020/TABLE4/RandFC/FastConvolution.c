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

/*//pour force inline
#undef force_inline
#define force_inline __attribute__((always_inline))*/


#define NSAMPLES 50

#define NTEST 1000

unsigned long long int START, STOP, START1,STOP1;


void initFWVector(uint64_t * nB,uint32_t * vB, int w, int size)
{

	int dec;//, tab[w];
	bool flag;
	
	for(int i=0;i<size;i++) nB[i]=0UL;
	
	for(int i=0;i<w;i++)
	{
		do
		{
			flag=false;
			dec = rand()%(PARAM_N);//16;//(i<<4)+i;//i;//7807;//20;//
			for(int j=0;j<i;j++)
			{
				if(dec==vB[j]) flag = true;
			}
			vB[i]=dec;
		}
		while(flag==true);
		
		nB[dec>>6] ^= 1UL<< (dec%WORD);

	
	
	}


}

/********************************************************************************
*
* MAIN
*
*********************************************************************************/



int main(int argc, char* argv[]){


	int flag=0, counter=0, t=SIZE64;

	#ifdef P67699
	unsigned long int nA[SIZE64], nB[SIZE64], res[SIZE64<<1], resMul[SIZE64<<1], tmp[SIZE64<<1];
	#elif P70853
	unsigned long int nA[SIZE64], nB[SIZE64], res[SIZE64<<1], resMul[SIZE64<<1], tmp[SIZE64<<1];
	
	#else
	
	static __m256i nA256[SIZE256], nB256[SIZE256], res256[SIZE256<<1], resMul256[SIZE256<<1], tmp256[SIZE256<<1];
	uint64_t* nA = (uint64_t *) nA256, *nB = (uint64_t *) nB256, *res = (uint64_t *) res256,
			 *resMul = (uint64_t *) resMul256, *tmp = (uint64_t *) tmp256;
 	#endif
	
	uint64_t mini = (uint64_t)-1L, mini1 = (uint64_t)-1L,mask=0UL;
	
	#ifdef P24677
	static uint32_t vB[PARAM_OMEGA];
	#else
	uint32_t vB[PARAM_OMEGA];
 	#endif
 
	unsigned long long int timer=0, timer1=0;
	
	
	
	printf("PARAM_N = %d, t = %d, t/(256/WORD) =%d\n",PARAM_N,t,t/(256/WORD));
	printf("%d\n",((PARAM_N/WORD) + (PARAM_N%WORD == 0 ? 0 : 1)));
	printf("PARAM_OMEGA = %d\n",PARAM_OMEGA);
	
	
	srand(time(NULL));

	for(int j=0; j<t;j++){
		nA[j] = 0x0UL;
		//nB[j] = (((unsigned long int)(rand()+rand())<<32)^(rand()+rand()));
	}

	for (int i=0;i<(PARAM_N>>6)+1;i++)
	{
		// pour avoir A[i] et B[i] entre 0x00000000 et 0xffffffff
		nA[i] = (((unsigned long int)(rand()+rand())<<32)^(rand()+rand()));
	}
	nA[(PARAM_N>>6)] &=mask;
	
		

	initFWVector(nB,vB,PARAM_OMEGA,t);
	

	
	afficheVect(nA,"nA",t);
	afficheVect(nB,"nB",t);


	/***********************************************/
	printf("\nsparse_dense_mul :\n----------\n");
	
	
	//printf("\nsortie fastConvolution :\n-----------------------------\n");
	
	
	sparse_dense_mul((uint8_t *)res,vB,(uint8_t *)nA,PARAM_OMEGA);
		
	afficheVect(res,"res",SIZE64);
		
	printf("\nComparaison avec fastConvolution :\n-----------------------------\n");
	
	fastConvolutionMult(nA,vB,resMul,PARAM_OMEGA);
	afficheVect(resMul,"resMul",SIZE64);
	
	for(int i=0; i<SIZE64;i++) tmp[i] = res[i]^resMul[i];
	
	printf("\n");
	afficheVect(tmp,"cmp",SIZE64);
	
	printf("\n");
	
	for(int i=0; i<SIZE64;i++)
		if(res[i]^resMul[i]) flag++;

	printf("flag = %d ; ",flag);
	
	if(!flag) printf("Victory !!!!!!!!!!!!!!!!!!!\n\n");
	else printf("Too bad !!!!!!!!!!!!!!!!!!!\n\n");
	
	flag=0;
	printf("\n\n");
	
	//goto fin;
	

	
	/********************
	//Comparaison sparse_dense_mul fastConvolution
	*********************/
	
	
	printf("\t  /***********************************/\n");
	printf("\t /   Test sur 1000 jeux de données   /\n");
	printf("\t/***********************************/\n\n");

	

	
	for(int i=0;i<NTEST;i++)
	{

		for (int i=0;i<(PARAM_N>>6);i++)
		{
			// pour avoir A[i] et B[i] entre 0x00000000 et 0xffffffff
			nA[i] = 0x0UL;
		}

		for (int i=0;i<(PARAM_N>>6);i++)
		{
			// pour avoir A[i] et B[i] entre 0x00000000 et 0xffffffff
			nA[i] = (((unsigned long int)(rand()+rand())<<32)^(rand()+rand()));
		}
		initFWVector(nB,vB,PARAM_OMEGA,t);
		
		fastConvolutionMult(nA,vB,resMul,PARAM_OMEGA);
		sparse_dense_mul((uint8_t *)res,vB,(uint8_t *)nA,PARAM_OMEGA);

		
		for(int i=0; i<SIZE64;i++)
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
	
	printf("\t\tPARAM_N = %d\n\t\t Size = %d 64 bit words\n\n",PARAM_N,SIZE64);
	printf("PARAM_OMEGA = %d\n",PARAM_OMEGA);
	

	printf("t = %d, t/(256/WORD) =%d, \n",t,t/(256/WORD));

	
	printf("\nsparse_dense_mul vs fastConvolution\n");
	printf("-----------------------------------\n");

	
	for(int k=0; k<NSAMPLES;k++){
	
		mini = (uint64_t)-1L, mini1 = (uint64_t)-1L;

		for (int i=0;i<(PARAM_N>>6);i++)
		{
			// pour avoir A[i] et B[i] entre 0x00000000 et 0xffffffff
			nA[i] = (((unsigned long int)(rand()+rand())<<32)^(rand()+rand()));
		}
		initFWVector(nB,vB,PARAM_OMEGA,t);
		
		for(int i=0;i<NTEST;i++)
		{
			sparse_dense_mul((uint8_t *)res,vB,(uint8_t *)nA,PARAM_OMEGA);

		}
		
		for(int i=0;i<NTEST;i++)
		{
			
			STAMP(START)
			sparse_dense_mul((uint8_t *)res,vB,(uint8_t *)nA,PARAM_OMEGA);
			STAMP(STOP)
			
			if(mini>STOP-START) mini = STOP-START;

		}

		for(int i=0;i<NTEST;i++)
		{
			fastConvolutionMult(nA,vB,resMul,PARAM_OMEGA);

		}
		
		for(int i=0;i<NTEST;i++)
		{

			STAMP(START1)
			fastConvolutionMult(nA,vB,resMul,PARAM_OMEGA);
			STAMP(STOP1)

			if(mini1>STOP1-START1) mini1 = STOP1-START1;
			
		}
		timer += mini;

		timer1 += mini1;
	}
	
	printf("timer sparse_dense_mul      : %llu\n",timer/NSAMPLES);	
	
	printf("timer fastConvolution  : %llu\n",timer1/NSAMPLES);



fin:
	printf("\n");
	printf("Au revoir et merci !\n\n");

}
