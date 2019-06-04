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


#define NTEST 1000

unsigned long long int START, STOP, START1,STOP1;

void initFWVector(uint64_t * nB,int * vB, int w, int size)
{

	int dec;//, tab[w];
	bool flag;
	
	for(int i=0;i<size;i++) nB[i]=0UL;
	
	for(int i=0;i<w;i++)
	{
		do
		{
			flag=false;
			dec = rand()%(size*WORD);
			//dec=63;
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
//#define ALIGNED


int main(){

#ifndef ALIGNED
	unsigned long int /*a, b, c, d, t1,tmp64[t],*///,t2_512;
	  nA[t], nB[t], res[t*2], resMul[t*2], tmp[t*2];
#else
		__m256i /*a, b, c, d, t1,tmp64[t],*///,t2_512;
	  nA256[t>>2], nB256[t>>2], res256[t>>1], resMul256[t>>1], tmp256[t>>1];
	unsigned long int *nA = (unsigned long int *) nA256;
	  
	unsigned long int *nB = (unsigned long int *) nB256;
	unsigned long int *res = (unsigned long int *) res256;
	unsigned long int *resMul = (unsigned long int *) resMul256;
	unsigned long int *tmp = (unsigned long int *) tmp256;
	
	printf("ALIGNED !!!!!!!\n\n");
	//getchar();
	  
#endif

	 int flag=0, counter=0, vB[W];
	 unsigned long long int timer, timer1;
	
	printf("FastConvolution vs polynomial multiplication\n\t\t Size =%d bits\n\n",SIZE_N);
	

	printf("t = %d, t/(256/WORD) =%d\n",t,t/(256/WORD));
	srand(time(NULL));

	for (int i=0;i<t;i++)
	{
		// pour avoir A[i] et B[i] entre 0x00000000 et 0xffffffff
		nA[i] = /*0xf0000000fUL;//*/(((unsigned long int)(rand()+rand()+rand()/RAND_MAX))<<32)^(rand()+rand()+rand()/RAND_MAX);//*/
		//nB[i] = (((unsigned long int)(rand()+rand()+rand()/RAND_MAX))<<32)^(rand()+rand()+rand()/RAND_MAX);
		//dec = rand()%64;
		//nB[i] = 1UL<< dec;
	}
		//nA[0] =(((unsigned long int)(rand()+rand()+rand()/RAND_MAX))<<32)^(rand()+rand()+rand()/RAND_MAX);

	initFWVector(nB,vB,W,t);
	
	//vB[W-1]=0;

	
	afficheVect(nA,"nA",t);
	//matVectProd700_350(H,nA,res);
	afficheVect(nB,"nB",t);
	
	printf("vB = ");
	
	for(int i=0;i<W;i++) printf("%d ",vB[i]);
	printf("\n\n");

	/***********************************************/
	printf("multPclmulSsRed\n-----------------\n");
	STAMP(START)
	{
		multPclmulSsRed(nA,nB,resMul,t);
		/*mul1n(&t2_512,&b,&t1_512);
		mul1n(&d512,&b,&t2_512);*/
	}
	
	STAMP(STOP)
	afficheVect(resMul,"resMul",t<<1);
	printf("\n\n");
	
	printf("STOP-START = %lld\n\n",(STOP-START));
	
	
	printf("\n\n");
	//res[t] = mpn_mul_1(res, nA,t,b);
	//KaratRecPclmul128((__m128i*)nA,(__m128i*)nB,(__m128i*)res,t>>2);
	
	/***********************************************/
	if(!(SIZE_N%10240))
	{
		printf("FiveWayMult\n--------------\n");
		STAMP(START)
		fiveWayMult(nA,nB,res);
		STAMP(STOP)
		//comparaison des deux résultats
	}
	else
	if(!(SIZE_N%6144))
	{
		printf("threeWayMult\n--------------\n");
		STAMP(START)
		threeWayMult(nA,nB,res);
		STAMP(STOP)
		//comparaison des deux résultats
	}
	else
	{
		printf("KaratRecPclmul\n--------------\n");
		STAMP(START)
		KaratRecPclmul(nA,nB,res,t);
		STAMP(STOP)
		//comparaison des deux résultats
	}
	afficheVect(res,"res",t<<1);
	printf("STOP-START = %lld\n\n",(STOP-START));
	
	printf("Comparaison :\n");
	
	for(int i=0; i<t<<1;i++) printf("%16.16lX ",res[(t<<1)-1-i]^resMul[(t<<1)-1-i]);
	
	printf("\n\n");
	
	for(int i=0; i<t<<1;i++)
		if(res[i]^resMul[i]) flag++;

	printf("flag = %d ; ",flag);
	
	if(!flag) printf("Victory !!!!!!!!!!!!!!!!!!!\n\n");
	else printf("Too bad !!!!!!!!!!!!!!!!!!!\n\n");
	
	flag=0;
	printf("\n\n");
	//for(int i=0; i<t<<1;i++) res[i]=0;
	
	
	/***********************************************/
	printf("fastConvolutionMult\n-------------------\n");
	STAMP(START)
	fastConvolutionMult(nA,vB,res,t);
	STAMP(STOP)
	//comparaison des deux résultats
	afficheVect(res,"res",t<<1);
	printf("STOP-START = %lld\n\n",(STOP-START));
	
	printf("Comparaison :\n");
	
	for(int i=0; i<t<<1;i++) tmp[i] = res[i]^resMul[i];
	
	printf("\n");
	afficheVect(tmp,"cmp",t<<1);
	
	/*for(int i=0; i<t<<1;i++) printf("%16.16lX ",res[(t<<1)-1-i]^resMul[(t<<1)-1-i]);
	printf("\n");
	for(int i=0; i<t<<1;i++) printf("%16.16lX ",res[i]^resMul[i]);*/
	
	
	printf("\n");
	
	for(int i=0; i<t<<1;i++)
		if(res[i]^resMul[i]) flag++;

	printf("flag = %d ; ",flag);
	
	if(!flag) printf("Victory !!!!!!!!!!!!!!!!!!!\n\n");
	else printf("Too bad !!!!!!!!!!!!!!!!!!!\n\n");
	
	flag=0;
	printf("\n\n");
	
	//goto fin;
	
	/********************
	//Comparaison fastConvolution vs KaratRecPclmul
	*********************/
	
	printf("Comparaison sur 1000 jeux de données\n");
	printf("  fastConvolution vs KaratRecPclmul \n");
	printf("------------------------------------\n");
	
	
	for(int i=0;i<1000;i++)
	{
		for(int j=0; j<t;j++) nA[j] = (((unsigned long int)(rand()+rand()+rand()/RAND_MAX))<<32)^(rand()+rand()+rand()/RAND_MAX);
		initFWVector(nB,vB,W,t);
		
		if(!(SIZE_N%10240))
			fiveWayMult(nA,nB,resMul);
		else
		if(!(SIZE_N%6144))
			threeWayMult(nA,nB,resMul);
		else KaratRecPclmul(nA,nB,resMul,t);
	
		fastConvolutionMult(nA,vB,res,t);
		
		for(int i=0; i<t<<1;i++)
			if(res[i]^resMul[i]) flag++;
		//printf("flag = %d ; ",flag);
		//getchar();*/
		flag?counter++,flag=0:counter,flag=0;
	
	}
	if(counter) printf("%d erreurs !\nToo bad !!!!!!!!!!!!!!!!!!!\n\n",counter),counter=0;
	else printf("Victory !!!!!!!!!!!!!!!!!!!\n\n");
	
	
	/********************
	//Chronométrage !!!!!!!!!!!!!
	*********************/
	
	//goto fin;
	printf("\t  /***************************/\n");
	printf("\t / Chronométrage !!!!!!!!!!!!/\n");
	printf("\t/***************************/\n\n");
	printf("FastConvolution vs polynomial multiplication\n\t\t Size =%d bits\n\t\t w_e =%d\n\n",SIZE_N,W);
	

	printf("t = %d, t/(256/WORD) =%d\n",t,t/(256/WORD));
	
	uint64_t mini = (uint64_t)-1L, mini1 = (uint64_t)-1L;
	
	if(!(SIZE_N%10240)) printf("multPclmulSsRed vs fiveWayMult\n");
		else
		if(!(SIZE_N%6144)) printf("multPclmulSsRed vs threeWayMult\n");
		else  printf("multPclmulSsRed vs KaratRecPclmul\n");
	printf("---------------------------------\n");
	
	for(int i=0;i<NTEST;i++)
	{
		for(int j=0; j<t;j++) nA[j] = (((unsigned long int)(rand()+rand()+rand()/RAND_MAX))<<32)^(rand()+rand()+rand()/RAND_MAX);
		initFWVector(nB,vB,W,t);
		//initPCMat(H);
		STAMP(START)
		multPclmulSsRed(nA,nB,res,t);
		STAMP(STOP)
		//tmp64 = d;
		//for(int i = 0;i<t;i++) resMul[i] = tmp64[i];
		if(!(SIZE_N%10240)){
			STAMP(START1)
			fiveWayMult(nA,nB,resMul);
			STAMP(STOP1)
		}	else
		if(!(SIZE_N%6144))
		{
			STAMP(START1)
			threeWayMult(nA,nB,resMul);
			STAMP(STOP1)
			//comparaison des deux résultats
		}

		else{
			STAMP(START1)
			KaratRecPclmul(nA,nB,resMul,t);
			STAMP(STOP1)
		}
			for(int i=0; i<t<<1;i++)
				if(res[i]^resMul[i]) flag++;
			flag?counter++,flag=0:counter,flag=0;
			
			if(mini>STOP-START) mini = STOP-START;

			if(mini1>STOP1-START1) mini1 = STOP1-START1;
		}
	timer = mini;

	timer1 = mini1;
	
	printf("Nombre de tests = %d\n",NTEST);
	if(counter){
		printf("%d erreurs !\n",counter);
		printf("Too bad !!!!!!!!!!!!!!!!!!!\n\n");
	}
	else printf("Victory !!!!!!!!!!!!!!!!!!!\n\n");

	
	printf("timer multPclmulSsRed : %llu\n",timer);
	if(!(SIZE_N%10240)) 	printf("timer fiveWayMult     :  %llu\n\n",timer1);
		else
		if(!(SIZE_N%6144)) printf("timer threeWayMult    :  %llu\n\n",timer1);
	else 				printf("timer KaratRecPclmul  :  %llu\n\n",timer1);	
	

	mini = (uint64_t)-1L, mini1 = (uint64_t)-1L;
	
	if(!(SIZE_N%10240)) printf("fiveWayMult vs fastConvolution\n");
		else
		if(!(SIZE_N%6144)) printf("threeWayMult vs fastConvolution\n");
	else printf("KaratRecPclmul vs fastConvolution\n");
	printf("---------------------------------\n");
	
	for(int i=0;i<NTEST;i++)
	{
		for(int j=0; j<t;j++) nA[j] = (((unsigned long int)(rand()+rand()+rand()/RAND_MAX))<<32)^(rand()+rand()+rand()/RAND_MAX);
		initFWVector(nB,vB,W,t);
		if(!(SIZE_N%10240)){
			STAMP(START)
			fiveWayMult(nA,nB,resMul);
			STAMP(STOP)
		}	else
		if(!(SIZE_N%6144))
		{
			STAMP(START)
			threeWayMult(nA,nB,resMul);
			STAMP(STOP)
			//comparaison des deux résultats
		}

		else{
			STAMP(START)
			KaratRecPclmul(nA,nB,resMul,t);
			STAMP(STOP)
		}
		STAMP(START1)
		fastConvolutionMult(nA,vB,res,t);
		STAMP(STOP1)
		for(int i=0; i<t<<1;i++)
			if(res[i]^resMul[i]) flag++;
		flag?counter++,flag=0:counter,flag=0;
		
		if(mini>STOP-START) mini = STOP-START;

		if(mini1>STOP1-START1) mini1 = STOP1-START1;
	}
	timer = mini;

	timer1 = mini1;
	
	printf("Nombre de tests = %d\n",NTEST);
	if(counter){
		printf("%d erreurs !\n",counter);
		printf("Too bad !!!!!!!!!!!!!!!!!!!\n\n");
	}
	else printf("Victory !!!!!!!!!!!!!!!!!!!\n\n");

	
	if(!(SIZE_N%10240)) 	printf("timer fiveWayMult     : %llu\n",timer);
		else
		if(!(SIZE_N%6144)) printf("timer threeWayMult    : %llu\n",timer);
	else printf("timer KaratRecPclmul  : %llu\n",timer);	
	printf("timer fastConvolution : %llu\n",timer1);


fin:
	printf("\n");//*/
	printf("Au revoir et merci !\n\n");

}
