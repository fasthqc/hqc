/*****************************************************************



*****************************************************************/



#include "fonctions.h"

/***************************************************************

	Fonctions d'affichage

*/
  

void afficheVect(unsigned long int *A, char *var, int size)
{
	int i;
	unsigned long int tmp;
	printf("%s := ",var);
	
	for(i=0;i<size;i++){
		tmp=0;
		for(int j=0;j<WORD;j++) tmp^= ((A[i]>>j)&1UL)<<(WORD-1-j);
		printf("%16.16lX ",tmp);
	}
	printf("\n");
}


/***************************************************************

	Reverses the bits: MSB -> LSB and vice-versa

*/


void reverse_u32(uint32_t * A, uint64_t * Out, int s)
{
	int i;
	uint32_t * out = (uint32_t *) Out;
	//printf("%s := ",var);
	
	for(i=0;i<s;i++){
		out[i]=0;
		for(int j=0;j<32;j++) out[i]^= ((A[i]>>j)&1)<<(WORD-1-j);
		//printf("%16.16lX ",tmp);
	}
	
	for(;i<t;i++){
		out[i]=0;
	}
	
	//printf("\n");
}



/*************************************************************************************
//
//                       MULTIPLICATION
//
//	Schoolbook 128 bits avec PCLMULQDQ intrinsic
//  multiplication élémentaire 256 bits
//                   SANS REDUCTION
//
*************************************************************************************/

int multPclmulSsRed(uint64_t * A, uint64_t * B, uint64_t * C, int size)//size est en nombre de mots de 64 bits et on doit avoir size%4=0 !!!
{
	__m128i tmp128[4];//,A1,B1;
	
	__m128i * A128 = (__m128i *) A;
	__m128i * B128 = (__m128i *) B;
	//__m128i * C128 = (__m128i *) C;
	
	__m128i C128[4*size];
	uint64_t * C64 = (uint64_t *) C128;
	
	int s = size/4;
	
	__m128i *C128prime = (__m128i *)(C64+1);
	for(int i=0;i<size;i++) C128[i]=(__m128i){0x0UL,0x0UL};

	for(int i=0;i<2*s;i++){
		for(int j=0;j<2*s;j++)
		{
			tmp128[0] = _mm_clmulepi64_si128 (A128[i], B128[j], 0x00);//a0b0 et a2b2
			tmp128[1] = _mm_clmulepi64_si128 (A128[i], B128[j], 0x01);//a1b0 et a3b2
			tmp128[2] = _mm_clmulepi64_si128 (A128[i], B128[j], 0x10);//a0b1 et a2b3
			tmp128[3] = _mm_clmulepi64_si128 (A128[i], B128[j], 0x11);//a1b1 et a3b3
			
			C128[i+j]^=tmp128[0];
			C128[i+j+1]^=tmp128[3];
			C128prime[i+j]^=tmp128[2]^tmp128[1];
			
		}
	}
	

	
	for(int i=0;i<size*2;i++) C[i] = C64[i];

	return 0;

}


//--------------------------------------------------------------------------------------------------
//
//                       MULTIPLICATION
//
//	Karatsuba PCLMULQDQ intrinsic (256 bits)
//                   SANS REDUCTION
//
//--------------------------------------------------------------------------------------------------


inline static int Pclmul256(__m128i * A, __m128i * B, __m128i * C)// Karatsuba_256_256_512
{


	__m128i DD0 , DD1, DD2,
			D0l, D1l, D2l, D0h, D1h, D2h,
	  AlpAh, BlpBh, AAlpAAh, BBlpBBh;


	//--------------------------------------------------------------------------------------------------
	//	Calcul de Al.Bl=D0 on le fait récursivement !
	//

	//
	// Boucle de calcul de DD0 = AAAl*BBBl
	//
	
	DD0 = _mm_clmulepi64_si128 (A[0], B[0], 0x00);
		
	// Calcul de AAAh.BBBh soit DD2

	DD2 = _mm_clmulepi64_si128 (A[0], B[0], 0x11);

	// Calcul de AAAlpAAAh.BBBlpBBBh soit DD1

	AAlpAAh = _mm_xor_si128(A[0],_mm_shuffle_epi32(A[0],0x4e));
	BBlpBBh = _mm_xor_si128(B[0],_mm_shuffle_epi32(B[0],0x4e));

	DD1 = _mm_clmulepi64_si128 (AAlpAAh,BBlpBBh, 0);

	//Additions finales

	D0l = _mm_xor_si128(DD0,_mm_xor_si128(_mm_slli_si128(DD0,8),
									_mm_xor_si128(_mm_slli_si128(DD2,8),
									_mm_slli_si128(DD1,8))));

	D0h = _mm_xor_si128(DD2,_mm_xor_si128(_mm_srli_si128(DD0,8),
									_mm_xor_si128(_mm_srli_si128(DD2,8),
									_mm_srli_si128(DD1,8))));
	



	//--------------------------------------------------------------------------------------------------
	// Calcul de Ah.Bh = D2 soit la taille d'un mot DD1 toujours par Karatsuba
	//

	//
	// Boucle de calcul de DD0 = AAAl*BBBl
	//
	
	DD0 = _mm_clmulepi64_si128 (A[1], B[1], 0x00);
		
	// Calcul de AAAh.BBBh soit DD2

	DD2 = _mm_clmulepi64_si128 (A[1], B[1], 0x11);

	// Calcul de AAAlpAAAh.BBBlpBBBh soit DD1

	AAlpAAh = _mm_xor_si128(A[1],_mm_shuffle_epi32(A[1],0x4e));
	BBlpBBh = _mm_xor_si128(B[1],_mm_shuffle_epi32(B[1],0x4e));

	DD1 = _mm_clmulepi64_si128 (AAlpAAh,BBlpBBh, 0);

	//Additions finales

	D2l = _mm_xor_si128(DD0,_mm_xor_si128(_mm_slli_si128(DD0,8),
									_mm_xor_si128(_mm_slli_si128(DD2,8),
									_mm_slli_si128(DD1,8))));

	D2h = _mm_xor_si128(DD2,_mm_xor_si128(_mm_srli_si128(DD0,8),
									_mm_xor_si128(_mm_srli_si128(DD2,8),
									_mm_srli_si128(DD1,8))));
	


	//--------------------------------------------------------------------------------------------------
	// Calcul de D1
	//
	//


	// initialisation de AlpAh et BlpBh

	AlpAh = _mm_xor_si128(A[0],A[1]);
	BlpBh = _mm_xor_si128(B[0],B[1]);



	// Calcul de AAl.BBl par Karatsuba

	DD0 = _mm_clmulepi64_si128(AlpAh,BlpBh,0);

	// Calcul de AAAh.BBBh soit DD2

	DD2 = _mm_clmulepi64_si128(AlpAh,BlpBh,0x11);

	// Calcul de AAAlpAAAh.BBBlpBBBh soit DD1

	AAlpAAh = _mm_xor_si128(AlpAh,_mm_shuffle_epi32(AlpAh,0x4e));
	BBlpBBh = _mm_xor_si128(BlpBh,_mm_shuffle_epi32(BlpBh,0x4e));

	DD1 = _mm_clmulepi64_si128 (AAlpAAh,BBlpBBh, 0);


	//Additions finales

	D1l = _mm_xor_si128(DD0,_mm_xor_si128(_mm_slli_si128(DD0,8),
									_mm_xor_si128(_mm_slli_si128(DD2,8),
									_mm_slli_si128(DD1,8))));

	D1h = _mm_xor_si128(DD2,_mm_xor_si128(_mm_srli_si128(DD0,8),
									_mm_xor_si128(_mm_srli_si128(DD2,8),
									_mm_srli_si128(DD1,8))));


	//--------------------------------------------------------------------------------------------------
	//		CALCUL FINAL DE C !!!!!!
	//


		C[0] = D0l;
		C[1] = _mm_xor_si128(D0h,_mm_xor_si128(D0l,_mm_xor_si128(D2l,D1l)));
		C[2] = _mm_xor_si128(D2l,_mm_xor_si128(D0h,_mm_xor_si128(D2h,D1h)));
		C[3] = D2h;


	return 0;

}



/*************************************************************************************
//
//                       MULTIPLICATION
//
//	Karatsuba récursif 256 bits avec PCLMULQDQ intrinsic
//  multiplication élémentaire 256 bits
//                   SANS REDUCTION
//
*************************************************************************************/

inline static int KaratRecPclmul256(__m256i * A, __m256i * B, __m256i * C, int size)//size est en nombre de mots de 256 bits !!!
{
	//printf("entrée dans KaratRecPclmul, size = %d\n",size);
	if(size==1){
		//KaratPclmul256(A,B,C);
		//afficheVect512((unsigned long int *) C,"C");
		Pclmul256((__m128i*)A,(__m128i*)B,(__m128i*)C);
		//multPclmul(A, B, C, size);
		
		//afficheVect512((unsigned long int *) C,"C");
		//afficheVect512((unsigned long int *) C,"C");
		
	}
	else{
		__m256i D0[size],D1[size],D2[size],
				SAA[size>>1],SBB[size>>1];

		size/=2;
		
		/*printf("size = %d ;\n",size);
		
		afficheVect((unsigned long int *) A,"A",4*size);
		afficheVect((unsigned long int *) B,"B",4*size);*/

		KaratRecPclmul256(A,B,D0,size);
		
		/*afficheVect((unsigned long int *) D0,"D0",8*size);
		getchar();
		printf("(unsigned long int *) &A[size] := %16.16lX\n",((unsigned long int *) &A[size])[0]);
		printf("(unsigned long int *) &B[size] := %16.16lX\n",((unsigned long int *) &B[size])[0]);
		afficheVect((unsigned long int*) &A[size],"&A[size]",4*size);
		afficheVect((unsigned long int*) &B[size],"&B[size]",4*size);//*/
		
		KaratRecPclmul256(&A[size],&B[size],D2,size);
		
		/*afficheVect((unsigned long int *) D2,"D2",8*size);
		getchar();*/
		
		for(int i=0;i<size;i++) {SAA[i]=A[i]^A[i+size];SBB[i]=B[i]^B[i+size];}
		
		/*printf("size = %d ;\n",size);
		afficheVect((unsigned long int *) SAA,"SAA",4*size);
		afficheVect((unsigned long int *) SBB,"SBB",4*size);//*/
		
		KaratRecPclmul256(SAA,SBB,D1,size);
		
		/*afficheVect((unsigned long int *) D0,"D0",8*size);
		afficheVect((unsigned long int *) D1,"D1",8*size);
		afficheVect((unsigned long int *) D2,"D2",8*size);*/
		
		for(int i=0;i<size;i++)
		{
			C[i] = D0[i];
			C[i+size] = D0[i+size]^D0[i]^D1[i]^D2[i];
			C[i+2*size] = D2[i]^D0[i+size]^D1[i+size]^D2[i+size];
			C[i+3*size] = D2[i+size];
		}
	}
	

	return 0;

}




/*************************************************************************************
//
//                       MULTIPLICATION
//
//	Wrapper 64 bits pour
//  Karatsuba récursif 128 bits avec PCLMULQDQ intrinsic
//                   SANS REDUCTION
//
*************************************************************************************/

int KaratRecPclmul(uint64_t * A64, uint64_t * B64, uint64_t * C64, int size)//size est en nombre de mots de 64 bits !!!
{

	__m256i * A = (__m256i *) A64;
	__m256i * B = (__m256i *) B64;

	__m256i * C = (__m256i *) C64;

	KaratRecPclmul256(A,B,C,size>>2);
		

	return 0;

}



/*************************************************************************************
//
//                       MULTIPLICATION
//
//	Wrapper 128 bits pour
//  Karatsuba récursif 128 bits avec PCLMULQDQ intrinsic
//                   SANS REDUCTION
//
*************************************************************************************/

int KaratRecPclmul128(__m128i *A128, __m128i *B128,  __m128i *C128,int size)//size est en nombre de mots de 256 bits !!!
{


	__m256i * A = (__m256i *) A128;
	__m256i * B = (__m256i *) B128;

	__m256i * C = (__m256i *) C128;

	KaratRecPclmul256(A,B,C,size);
		

	return 0;

}


/*************************************************************************************
//
//                       MULTIPLICATION
//
//	Multiplication 3 way, 43669, opérandes vector_u32 (voir hqc) bits
//  puis multiplications KARATSUBA récursive 16384 (*6)
//  multiplication élémentaire 256 bits avec PCLMULQDQ intrinsic
//                   AVEC REDUCTION
//
*************************************************************************************/
#define T_3W_256 (SIZE_N/(768))

#define T2_3W_256  (2*T_3W_256)
#define t3 t

inline int threeWayMult(uint64_t* A, uint64_t* B, uint64_t* Out)//size est en nombre de mots de 256 bits !!!
{
	
	__m256i a0[T_3W_256], b0[T_3W_256], a1[T_3W_256], b1[T_3W_256], a2[T_3W_256], b2[T_3W_256];
	
	__m256i aa01[T_3W_256], bb01[T_3W_256], aa02[T_3W_256], bb02[T_3W_256], aa12[T_3W_256], bb12[T_3W_256];
	
	__m256i D0[T2_3W_256], D1[T2_3W_256], D2[T2_3W_256], D3[T2_3W_256], D4[T2_3W_256], D5[T2_3W_256];
	
	__m256i ro256[t/2];
	

	
	
	for(int i=0;i<T_3W_256;i++)
	{
		a0[i]= (__m256i){A[4*i],A[4*i+1],A[4*i+2],A[4*i+3]};//a0[i]= A256[i]
		b0[i]= (__m256i){B[4*i],B[4*i+1],B[4*i+2],B[4*i+3]};//b0[i]= B256[i];
		a1[i]= (__m256i){A[4*(i+T_3W_256)],A[4*(i+T_3W_256)+1],A[4*(i+T_3W_256)+2],A[4*(i+T_3W_256)+3]};//a1[i]= A256[i+T_3W_256];
		b1[i]= (__m256i){B[4*(i+T_3W_256)],B[4*(i+T_3W_256)+1],B[4*(i+T_3W_256)+2],B[4*(i+T_3W_256)+3]};//b1[i]= B256[i+T_3W_256];
		a2[i]= (__m256i){A[4*(i+2*T_3W_256)],A[4*(i+2*T_3W_256)+1],A[4*(i+2*T_3W_256)+2],A[4*(i+2*T_3W_256)+3]};//a2[i]= A256[i+2*T_3W_256];
		b2[i]= (__m256i){B[4*(i+2*T_3W_256)],B[4*(i+2*T_3W_256)+1],B[4*(i+2*T_3W_256)+2],B[4*(i+2*T_3W_256)+3]};//b2[i]= B256[i+2*T_3W_256];
	}
	
	
	
	for(int i=0;i<T_3W_256;i++)
	{
		//printf("i = %d\n",i);
		aa01[i]=a0[i]^a1[i];
		bb01[i]=b0[i]^b1[i];
		
		//printf("i = %d\n",i);
		aa12[i]=a2[i]^a1[i];
		bb12[i]=b2[i]^b1[i];
		
		//printf("i = %d\n",i);
		aa02[i]=a0[i]^a2[i];
		bb02[i]=b0[i]^b2[i];
	}
	
	
	
	KaratRecPclmul256(a0, b0, D0,T_3W_256);
	KaratRecPclmul256( a1, b1, D1,T_3W_256);
	KaratRecPclmul256( a2, b2, D2,T_3W_256);
	
	KaratRecPclmul256( aa01, bb01, D3,T_3W_256);
	KaratRecPclmul256( aa02, bb02, D4,T_3W_256);
	KaratRecPclmul256( aa12, bb12, D5,T_3W_256);
	
	
	
	for(int i=0;i<T_3W_256;i++)
	{
		ro256[i]            = D0[i];
		ro256[i+T_3W_256]   = D3[i]^D0[i]^D1[i]^D0[i+T_3W_256];
		ro256[i+2*T_3W_256] = D4[i]^D0[i]^D1[i]^D2[i]^D3[i+T_3W_256]^D0[i+T_3W_256]^D1[i+T_3W_256];
		ro256[i+3*T_3W_256] = D5[i]^D1[i]^D2[i]^D4[i+T_3W_256]^D0[i+T_3W_256]^D1[i+T_3W_256]^D2[i+T_3W_256];
		ro256[i+4*T_3W_256] = D2[i]^D5[i+T_3W_256]^D1[i+T_3W_256]^D2[i+T_3W_256];
		ro256[i+5*T_3W_256] = D2[i+T_3W_256];
	}


	

	uint64_t * ro64 = (uint64_t *) ro256;
	
	for(int i=0;i<2*t;i++) Out[i]=ro64[i];

	//reduction40597((uint64_t *) row256, Out);


	return 0;
}



/*************************************************************************************
//
//                       MULTIPLICATION
//
//	Multiplication 5 way, 40597, opérandes vector_u32 (voir hqc) bits
//  puis multiplications KARATSUBA récursive 8192 (*15)
//  multiplication élémentaire 256 bits avec PCLMULQDQ intrinsic
//                   AVEC REDUCTION
//
*************************************************************************************/
#define T_5W_256 (SIZE_N/1280)

#define T2_5W_256  (2*T_5W_256)
#define t5 t
int fiveWayMult(uint64_t* A, uint64_t* B, uint64_t* Out)//size est en nombre de mots de 256 bits !!!
{
	/*__m256i *A256 = (__m256i *)A, *B256 = (__m256i *)B;
	
	__m256i * a0, * b0, *a1, *b1, *a2, *b2;
	
	__m256i * aa01, * bb01, *aa12, *bb12, *aa02, *bb02;
	
	__m256i * D0, * D1, *D2, *D3, *D4, *D5;*/
	
	__m256i a0[T_5W_256], b0[T_5W_256], a1[T_5W_256], b1[T_5W_256], a2[T_5W_256], b2[T_5W_256],
			a3[T_5W_256], b3[T_5W_256], a4[T_5W_256], b4[T_5W_256];
	
	__m256i aa01[T_5W_256], bb01[T_5W_256], aa02[T_5W_256], bb02[T_5W_256], aa03[T_5W_256], bb03[T_5W_256], aa04[T_5W_256], bb04[T_5W_256], 
			aa12[T_5W_256], bb12[T_5W_256], aa13[T_5W_256], bb13[T_5W_256], aa14[T_5W_256], bb14[T_5W_256],
			aa23[T_5W_256], bb23[T_5W_256], aa24[T_5W_256], bb24[T_5W_256],
			aa34[T_5W_256], bb34[T_5W_256];
	
	__m256i D0[T2_5W_256], D1[T2_5W_256], D2[T2_5W_256], D3[T2_5W_256], D4[T2_5W_256], 
			D01[T2_5W_256], D02[T2_5W_256], D03[T2_5W_256], D04[T2_5W_256],
			D12[T2_5W_256], D13[T2_5W_256], D14[T2_5W_256],
			D23[T2_5W_256], D24[T2_5W_256],
			D34[T2_5W_256];
	
	__m256i ro256[t5/2];
	
	// Ne marche pas : erreur de segmentation !!!!!!!
	/*__m256i *ro256 = (__m256i*) Out;
	a0 = A256;
	a1 = A256+T_5W_256;
	a2 = A256+(T_5W_256<<1);
	
	b0 = B256;
	b1 = B256+T_5W_256;
	b2 = B256+(T_5W_256<<1);*/
	
	/*printf("affectations a0,b0, etc.\n");
	printf("a0 = %lX\n",(long)a0);
	printf("b0 = %lX\n",(long)b0);
	printf("a1 = %lX\n",(long)a1);
	printf("b1 = %lX\n",(long)b1);
	printf("a2 = %lX\n",(long)a2);
	printf("b2 = %lX\n",(long)b2);*/
	
	
	for(int i=0;i<T_5W_256;i++)
	{
		a0[i]= (__m256i){A[4*i],A[4*i+1],A[4*i+2],A[4*i+3]};
		b0[i]= (__m256i){B[4*i],B[4*i+1],B[4*i+2],B[4*i+3]};
		a1[i]= (__m256i){A[4*(i+T_5W_256)],A[4*(i+T_5W_256)+1],A[4*(i+T_5W_256)+2],A[4*(i+T_5W_256)+3]};
		b1[i]= (__m256i){B[4*(i+T_5W_256)],B[4*(i+T_5W_256)+1],B[4*(i+T_5W_256)+2],B[4*(i+T_5W_256)+3]};
		a2[i]= (__m256i){A[4*(i+2*T_5W_256)],A[4*(i+2*T_5W_256)+1],A[4*(i+2*T_5W_256)+2],A[4*(i+2*T_5W_256)+3]};
		b2[i]= (__m256i){B[4*(i+2*T_5W_256)],B[4*(i+2*T_5W_256)+1],B[4*(i+2*T_5W_256)+2],B[4*(i+2*T_5W_256)+3]};
		a3[i]= (__m256i){A[4*(i+3*T_5W_256)],A[4*(i+3*T_5W_256)+1],A[4*(i+3*T_5W_256)+2],A[4*(i+3*T_5W_256)+3]};
		b3[i]= (__m256i){B[4*(i+3*T_5W_256)],B[4*(i+3*T_5W_256)+1],B[4*(i+3*T_5W_256)+2],B[4*(i+3*T_5W_256)+3]};
		a4[i]= (__m256i){A[4*(i+4*T_5W_256)],A[4*(i+4*T_5W_256)+1],A[4*(i+4*T_5W_256)+2],A[4*(i+4*T_5W_256)+3]};
		b4[i]= (__m256i){B[4*(i+4*T_5W_256)],B[4*(i+4*T_5W_256)+1],B[4*(i+4*T_5W_256)+2],B[4*(i+4*T_5W_256)+3]};
	}
	
	
	
	
		
	/*printf("calcul des aaxx,bbxx\n");
	fflush(stdout);*/
	
	for(int i=0;i<T_5W_256;i++)
	{
		//printf("i = %d\n",i);
		aa01[i]=a0[i]^a1[i];
		bb01[i]=b0[i]^b1[i];
		
		//printf("i = %d\n",i);
		aa02[i]=a0[i]^a2[i];
		bb02[i]=b0[i]^b2[i];
		
		//printf("i = %d\n",i);
		aa03[i]=a0[i]^a3[i];
		bb03[i]=b0[i]^b3[i];
		
		//printf("i = %d\n",i);
		aa04[i]=a0[i]^a4[i];
		bb04[i]=b0[i]^b4[i];
		
		//printf("i = %d\n",i);
		aa12[i]=a2[i]^a1[i];
		bb12[i]=b2[i]^b1[i];
		
		//printf("i = %d\n",i);
		aa13[i]=a3[i]^a1[i];
		bb13[i]=b3[i]^b1[i];
		
		//printf("i = %d\n",i);
		aa14[i]=a4[i]^a1[i];
		bb14[i]=b4[i]^b1[i];
		
		//printf("i = %d\n",i);
		aa23[i]=a2[i]^a3[i];
		bb23[i]=b2[i]^b3[i];
		
		//printf("i = %d\n",i);
		aa24[i]=a2[i]^a4[i];
		bb24[i]=b2[i]^b4[i];
		
		//printf("i = %d\n",i);
		aa34[i]=a3[i]^a4[i];
		bb34[i]=b3[i]^b4[i];	
		
	}
	
	
	/*printf("calcul des Dx\n");
	fflush(stdout);*/
	
	KaratRecPclmul256( a0, b0, D0,T_5W_256);
	KaratRecPclmul256( a1, b1, D1,T_5W_256);
	KaratRecPclmul256( a2, b2, D2,T_5W_256);
	KaratRecPclmul256( a3, b3, D3,T_5W_256);
	KaratRecPclmul256( a4, b4, D4,T_5W_256);
	
	KaratRecPclmul256( aa01, bb01, D01,T_5W_256);
	KaratRecPclmul256( aa02, bb02, D02,T_5W_256);
	KaratRecPclmul256( aa03, bb03, D03,T_5W_256);
	KaratRecPclmul256( aa04, bb04, D04,T_5W_256);
	
	KaratRecPclmul256( aa12, bb12, D12,T_5W_256);
	KaratRecPclmul256( aa13, bb13, D13,T_5W_256);
	KaratRecPclmul256( aa14, bb14, D14,T_5W_256);
	
	KaratRecPclmul256( aa23, bb23, D23,T_5W_256);
	KaratRecPclmul256( aa24, bb24, D24,T_5W_256);
	
	KaratRecPclmul256( aa34, bb34, D34,T_5W_256);
	
	
	/*printf("reconstruction finale\n");
	fflush(stdout);*/
	
	
	for(int i=0;i<T_5W_256;i++)
	{
		ro256[i]            = D0[i];
		ro256[i+T_5W_256]   = D0[i+T_5W_256]^D01[i]^D0[i]^D1[i];
		ro256[i+2*T_5W_256] = D1[i]^D02[i]^D0[i]^D2[i]^D01[i+T_5W_256]^D0[i+T_5W_256]^D1[i+T_5W_256];
		ro256[i+3*T_5W_256] = D1[i+T_5W_256]^D03[i]^D0[i]^D3[i]^D12[i]^D1[i]^D2[i]^D02[i+T_5W_256]^D0[i+T_5W_256]^D2[i+T_5W_256];
		ro256[i+4*T_5W_256] = D2[i]^D04[i]^D0[i]^D4[i]^D13[i]^D1[i]^D3[i]^D03[i+T_5W_256]^D0[i+T_5W_256]^D3[i+T_5W_256]^D12[i+T_5W_256]^D1[i+T_5W_256]^D2[i+T_5W_256];
		ro256[i+5*T_5W_256] = D2[i+T_5W_256]^D14[i]^D1[i]^D4[i]^D23[i]^D2[i]^D3[i]^D04[i+T_5W_256]^D0[i+T_5W_256]^D4[i+T_5W_256]^D13[i+T_5W_256]^D1[i+T_5W_256]^D3[i+T_5W_256];
		ro256[i+6*T_5W_256] = D3[i]^D24[i]^D2[i]^D4[i]^D14[i+T_5W_256]^D1[i+T_5W_256]^D4[i+T_5W_256]^D23[i+T_5W_256]^D2[i+T_5W_256]^D3[i+T_5W_256];
		ro256[i+7*T_5W_256] = D3[i+T_5W_256]^D34[i]^D3[i]^D4[i]^D24[i+T_5W_256]^D2[i+T_5W_256]^D4[i+T_5W_256];
		ro256[i+8*T_5W_256] = D4[i]^D34[i+T_5W_256]^D3[i+T_5W_256]^D4[i+T_5W_256];
		ro256[i+9*T_5W_256] = D4[i+T_5W_256];
	}


	uint64_t * ro64 = (uint64_t *) ro256;
	
	for(int i=0;i<2*t5;i++) Out[i]=ro64[i];

	//reduction43669((uint64_t *) row256, Out);


	return 0;
}




/*************************************************************************************
//
//                       MULTIPLICATION
//
//	Fast convolution 64 bits avec left shift
//                   SANS REDUCTION
//
*************************************************************************************/


/*

On multiplie sur F_2[x] un polynôme creux avec un polynôme de densité quelconque de la façon suivante :

vB = \sum_0^(W-1) x^(vB[i]) (polynôme creux de poids de Hamming W)
A (polynôme représenté par un tableau de mots de 64 bits) A[0] contient les 64 premiers monômes de degré 0 à 63.

C = \sum_0^(W-1) x^(vB[i])*A

c'est à dire qu'on décale à gauche de vB[i] le polynôme A et on additionne modulo 2 (xor bit à bit) chacun de ces décalages.

*/


// VERSION REGULIERE, mais vectorisée et plus rapide !!!

/******************
version entièrement vectorisée avec utilisation d'une petite tabulation pour les décalages des "extrémités"
empreinte mémoire : 320 octets (10*8)
*/

int fastConvolutionMult(uint64_t * A, int * vB, uint64_t * C, int size)//size est en nombre de mots de 64 bits !!!
{
	int dec,d0,dec64;
	
	__m256i A256;
	__m256i * C256 = (__m256i *) C;
	__m256i carry256;
	
	int size256 = size>>2;
	int rsize = (size<<1)%4,i64;
	
	for(int j=0;j<size256<<1;j++) C256[j] = (__m256i){0x0UL,0x0UL,0x0UL,0x0UL};
	
	for(int j=(size<<1)-rsize;j<size<<1;j++) C[j] = 0x0UL;
	
	//pour la troisième version sans switch : tableau !
	
	__m256i tab[5], tabprime[5];
	tab[0] = (__m256i){A[0],A[1],A[2],A[3]};
	tab[1] = (__m256i){0x0UL,A[0],A[1],A[2]};
	tab[2] = (__m256i){0x0UL,0x0UL,A[0],A[1]};
	tab[3] =(__m256i){0x0UL,0x0UL,0x0UL,A[0]};
	tab[4] = (__m256i){0x0UL,0x0UL,0x0UL,0x0UL};
	
	tabprime[0] = (__m256i){0x0UL,0x0UL,0x0UL,0x0UL};
	tabprime[1] = (__m256i){A[size-1],0x0UL,0x0UL,0x0UL};
	tabprime[2] = (__m256i){A[size-2],A[size-1],0x0UL,0x0UL};
	tabprime[3] =(__m256i){A[size-3],A[size-2],A[size-1],0x0UL};
	tabprime[4] = (__m256i){A[size-4],A[size-3],A[size-2],A[size-1]};
	

	for(int i=0;i<W;i++){
		dec = vB[i]&0xff;
		dec64 = dec&0x3f;
		i64=(dec>>6);
		d0=WORD-dec64;
		
		// prologue de la boucle de décalage
		// pour les shifts vectorisés ( _mm256_slli_epi64() ou _mm256_srli_epi64()), il n'y a pas nécessité de vérifier la taille du décalage < 64.
		
		A256 = tab[i64];
		carry256 = tab[i64+1];
		A256 = _mm256_slli_epi64(A256,dec64);
		carry256 = _mm256_srli_epi64(carry256,d0);
		A256^=carry256;
		C256[vB[i]>>8] ^= A256;
		
		// boucle vectorisée de décalage
		for(int j=1;j<size256;j++)
		{
			A256 = (__m256i){A[4*j-i64],A[4*j-i64+1],A[4*j-i64+2],A[4*j-i64+3]};
			A256 = _mm256_slli_epi64(A256,dec64);
			carry256 = (__m256i){A[4*j-i64-1],A[4*j-i64],A[4*j-i64+1],A[4*j-i64+2]};
			carry256 = _mm256_srli_epi64(carry256,d0);
			A256^=carry256;
			C256[j+(vB[i]>>8)] ^= A256;
		}	
	
		
		// épilogue de la boucle de décalage
		// pour les shifts vectorisés ( _mm256_slli_epi64() ou _mm256_srli_epi64()), il n'y a pas nécessité de vérifier la taille du décalage < 64.
		
		A256 = tabprime[i64];
		carry256 = tabprime[i64+1];
		A256 = _mm256_slli_epi64(A256,dec64);
		carry256 = _mm256_srli_epi64(carry256,d0);
		A256^=carry256;
		C256[size256+(vB[i]>>8)] ^= A256;
		
	}
	return 0;

}




