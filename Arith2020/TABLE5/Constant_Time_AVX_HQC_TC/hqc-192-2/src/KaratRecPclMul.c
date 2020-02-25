

#define SIZE_N 48768

#define tTM3 (SIZE_N/WORD)//(T_TM3/WORD)
#define T_TM3_3W (SIZE_N/3)
#define T_TM3_3W_256 ((T_TM3_3W+128)/(4*WORD))// on rajoute 128 bits pour avoir un nombre entier de mots de 256 bits.
#define T_TM3_3W_64 (T_TM3_3W_256<<2)











//--------------------------------------------------------------------------------------------------
//
//                       Fixed weight vector converion
//
//	Convert a fixed weight vector stored by positions to a coordinate vector
//                   of size PARAM_N bits (plus the padding...)
//
//--------------------------------------------------------------------------------------------------


void initFWVector(uint64_t * nB,uint32_t * vB, int w)
{
	
	//for(int i=0;i<SIZE64;i++) nB[i]=0UL;
	
	for(int i=0;i<w;i++)
		nB[vB[i]>>6] ^= 1UL<< (vB[i]%WORD);
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

	if(size==1){

		Pclmul256((__m128i*)A,(__m128i*)B,(__m128i*)C);

		
	}
	else{
		__m256i D0[size],D1[size],D2[size],
				SAA[size>>1],SBB[size>>1];

		size/=2;


		KaratRecPclmul256(A,B,D0,size);
		

		
		KaratRecPclmul256(&A[size],&B[size],D2,size);
		

		
		for(int i=0;i<size;i++) {int is = i+size; SAA[i]=A[i]^A[is];SBB[i]=B[i]^B[is];}
		

		KaratRecPclmul256(SAA,SBB,D1,size);
		

		
		for(int i=0;i<size;i++)
		{
			int is = i+size;
			C[i] = D0[i];
			C[is] = D0[is]^D0[i]^D1[i]^D2[i];
			C[i+2*size] = D2[i]^D0[is]^D1[is]^D2[is];
			C[i+3*size] = D2[is];
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

inline static int KaratRecPclmul(uint64_t * A64, uint64_t * B64, uint64_t * C64, int size)//size est en nombre de mots de 64 bits !!!
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
//	Multiplication TOOM3_5, opérandes vector_u32 (voir hqc) bits
//  evaluation en 0,1,x,x+1,\infty (inspiré de Bodrato et Zimmermann (ntl))
//	interpolation utilisant :
//  	division par x par shift à droite d'un mot (256 bits, x = X^256)
//  	division par x+1 par algo type Quercia
//  puis multiplications KARATSUBA récursive 16384 (*5, ou 8192 ou 4096)
//  multiplication élémentaire 256 bits avec PCLMULQDQ intrinsic
//                   SANS REDUCTION
//
*************************************************************************************/




// avec décalage d'un mot de 64 bits...

static inline int divByXplus1(__m256i* in,__m256i* out,int size){//mod X^N !!! N = 


	uint64_t * A = (uint64_t *) in;

	uint64_t * B = (uint64_t *) out;
	

	
	B[0] = A[0];
	
	for(int i=1;i<2*(size<<2);i++)
		B[i]= B[i-1]^A[i];
	


	return 0;
}

int TOOM3Mult(uint64_t* A, uint64_t* B, uint64_t* Out)
{


	
	static __m256i U0[T_TM3_3W_256], V0[T_TM3_3W_256], U1[T_TM3_3W_256], V1[T_TM3_3W_256], U2[T_TM3_3W_256], V2[T_TM3_3W_256];
	
	static __m256i W0[2*(T_TM3_3W_256)], W1[2*(T_TM3_3W_256)], W2[2*(T_TM3_3W_256)], W3[2*(T_TM3_3W_256)], W4[2*(T_TM3_3W_256)];
	static __m256i tmp[2*(T_TM3_3W_256)];// 
	
	static __m256i ro256[6*(T_TM3_3W_256)];
 
	const __m256i zero = (__m256i){0ul,0ul,0ul,0ul};
	
	int T2 = T_TM3_3W_64<<1;
	for(int i=0;i<T_TM3_3W_256-1;i++)
	{
		
		int i4 = i<<2;
		int i42 = i4-2;
		U0[i]= _mm256_lddqu_si256((__m256i const *)(& A[i4]));
		V0[i]= _mm256_lddqu_si256((__m256i const *)(& B[i4]));
		U1[i]= _mm256_lddqu_si256((__m256i const *)(& A[i42+T_TM3_3W_64]));
		V1[i]= _mm256_lddqu_si256((__m256i const *)(& B[i42+T_TM3_3W_64]));
		U2[i]= _mm256_lddqu_si256((__m256i const *)(& A[i4+T2-4]));
		V2[i]= _mm256_lddqu_si256((__m256i const *)(& B[i4+T2-4]));
	}
	
	for(int i=T_TM3_3W_256-1;i<T_TM3_3W_256;i++)
	{
		int i4 = i<<2;
		int i41 = i4+1;
		
		U0[i]= (__m256i){A[i4],A[i41],0x0ul,0x0ul};
		V0[i]= (__m256i){B[i4],B[i41],0x0ul,0x0ul};
		U1[i]= (__m256i){A[i4+T_TM3_3W_64-2],A[i41+T_TM3_3W_64-2],0x0ul,0x0ul};
		V1[i]= (__m256i){B[i4+T_TM3_3W_64-2],B[i41+T_TM3_3W_64-2],0x0ul,0x0ul};
		U2[i]= (__m256i){A[i4-4+T2],A[i4-3+T2],0x0ul,0x0ul};
		V2[i]= (__m256i){B[i4-4+T2],B[i4-3+T2],0x0ul,0x0ul};

	}
	
	
	// EVALUATION PHASE : x= X^64
	// P(X): P0=(0); P1=(1); P2=(x); P3=(1+x); P4=(\infty)
	// Evaluation: 5*2 add, 2*2 shift; 5 mul (n)
	
	
	//W3 = U2 + U1 + U0 ; W2 = V2 + V1 + V0

	for(int i=0;i<T_TM3_3W_256;i++)
	{
		//printf("i = %d\n",i);
		W3[i]=U0[i]^U1[i]^U2[i];
		W2[i]=V0[i]^V1[i]^V2[i];
	}
	
	//W1 = W2 * W3
	KaratRecPclmul256( W2, W3, W1,T_TM3_3W_256);
	
	
	//W0 =(U1 + U2*x)*x ; W4 =(V1 + V2*x)*x (SIZE = T_TM3_3W_256 !)
	//printf("W0 =(U1 + U2*x)*x ; W4 =(V1 + V2*x)*x !!!!!!!!\n");
	// décaler de 64 bits, x = X^64 !!!!!!!!!!!!!!!!!!!!!!!!!
	
	uint64_t * U1_64 = ((uint64_t *) U1);
	uint64_t * U2_64 = ((uint64_t *) U2);
	
	uint64_t * V1_64 = ((uint64_t *) V1);
	uint64_t * V2_64 = ((uint64_t *) V2);
	
	W0[0] = (__m256i){0ul,U1_64[0],U1_64[1]^U2_64[0],U1_64[2]^U2_64[1]};
	W4[0] = (__m256i){0ul,V1_64[0],V1_64[1]^V2_64[0],V1_64[2]^V2_64[1]};
	
	U1_64 = ((uint64_t *) U1)-1;
	U2_64 = ((uint64_t *) U2)-2;
	
	V1_64 = ((uint64_t *) V1)-1;
	V2_64 = ((uint64_t *) V2)-2;
	
	for(int i=1;i<T_TM3_3W_256;i++)
	{
		int i4 = i<<2;
		W0[i] = _mm256_lddqu_si256((__m256i const *)(& U1_64[i4]));
		W0[i] ^= _mm256_lddqu_si256((__m256i const *)(& U2_64[i4]));
		
		W4[i] = _mm256_lddqu_si256((__m256i const *)(& V1_64[i4]));
		W4[i] ^= _mm256_lddqu_si256((__m256i const *)(& V2_64[i4]));
	}
	
	

	
	//W3 = W3 + W0      ; W2 = W2 + W4
	for(int i=0;i<T_TM3_3W_256;i++)
	{
		W3[i] ^= W0[i];
		W2[i] ^= W4[i];
	}
	

	//W0 = W0 + U0      ; W4 = W4 + V0
	for(int i=0;i<T_TM3_3W_256;i++)
	{
		W0[i] ^= U0[i];
		W4[i] ^= V0[i];
	}

	//W3 = W3 * W2      ; W2 = W0 * W4
	KaratRecPclmul256( W3, W2, tmp,T_TM3_3W_256);
	for(int i=0;i<2*(T_TM3_3W_256);i++) W3[i] = tmp[i];
	KaratRecPclmul256( W0, W4, W2,T_TM3_3W_256);
	

	//W4 = U2 * V2      ; W0 = U0 * V0
	KaratRecPclmul256( U2, V2, W4,T_TM3_3W_256);
	KaratRecPclmul256( U0, V0, W0,T_TM3_3W_256);

	
	
	
	//INTERPOLATION PHASE
	//9 add, 1 shift, 1 Smul, 2 Sdiv (2n)
	
	//W3 = W3 + W2
	for(int i=0;i<2*(T_TM3_3W_256);i++)
		W3[i] ^= W2[i];
		
	
	//W1 = W1 + W0
	for(int i=0;i<2*(T_TM3_3W_256);i++)
		W1[i] ^= W0[i];
		
	
	//W2 =(W2 + W0)/x -> x = X^64
	
	U1_64 = ((uint64_t *) W2)+1;
	U2_64 = ((uint64_t *) W0)+1;
	for(int i=0;i<(T_TM3_3W_256<<1);i++)
		{
			int i4 = i<<2;
			W2[i] = _mm256_lddqu_si256((__m256i const *)(& U1_64[i4]));
			W2[i] ^= _mm256_lddqu_si256((__m256i const *)(& U2_64[i4]));
		}
	
	

	//W2 =(W2 + W3 + W4*(x^3+1))/(x+1)
	
	U1_64 = ((uint64_t *) W4);
	__m256i * U1_256 = (__m256i *) (U1_64+1);
	
	tmp[0] = W2[0]^W3[0]^W4[0]^(__m256i){0x0ul,0x0ul,0x0ul,U1_64[0]};
	for(int i=1;i<(T_TM3_3W_256<<1)-1;i++)
		tmp[i] = W2[i]^W3[i]^W4[i]^U1_256[i-1];
		
	
	divByXplus1(tmp,W2,T_TM3_3W_256);
	W2[2*(T_TM3_3W_256)-1] = zero;
	
	//W3 =(W3 + W1)/(x*(x+1))
	
	U1_64 = (uint64_t *) W3;
	U1_256 = (__m256i *) (U1_64+1);
	
	U2_64 = (uint64_t *) W1;
	__m256i * U2_256 = (__m256i *) (U2_64+1);
	
	for(int i=0;i<2*(T_TM3_3W_256)-1;i++)
		{tmp[i] = U1_256[i]^U2_256[i];}
	
	divByXplus1(tmp,W3,T_TM3_3W_256);
	W3[2*(T_TM3_3W_256)-1] = zero;
	
	
	//W1 = W1 + W4 + W2
	for(int i=0;i<2*(T_TM3_3W_256);i++)
		W1[i] ^= W2[i]^W4[i];
	
	//W2 = W2 + W3
	for(int i=0;i<2*(T_TM3_3W_256);i++)
		W2[i] ^= W3[i];
	
	
	
	// Recomposition
	//W  = W0+ W1*x+ W2*x^2+ W3*x^3 + W4*x^4
	//Attention : W0, W1, W4 of size 2*T_TM3_3W_256, W2 and W3 of size 2*(T_TM3_3W_256)


	for(int i=0;i<(T_TM3_3W_256<<1)-1;i++)
	{
		ro256[i]=W0[i];
		ro256[i+2*T_TM3_3W_256-1] = W2[i];
		ro256[i+4*T_TM3_3W_256-2] = W4[i];//*/
	}
	
	ro256[(T_TM3_3W_256<<1)-1]=W0[(T_TM3_3W_256<<1)-1]^W2[0];
	ro256[(T_TM3_3W_256<<2)-2]=W2[(T_TM3_3W_256<<1)-1]^W4[0];
	ro256[(T_TM3_3W_256*6) -3]=W4[(T_TM3_3W_256<<1)-1];
	
	U1_64 = ((uint64_t *) &ro256[T_TM3_3W_256]);
	U1_256 = (__m256i *) (U1_64-2);
	
	U2_64 = ((uint64_t *) &ro256[3*T_TM3_3W_256-1]);
	U2_256 = (__m256i *) (U2_64-2);
	
	for(int i=0;i<T_TM3_3W_256<<1;i++){
		U1_256[i]^=W1[i];
		U2_256[i]^=W3[i];
	}
	
	
	reduction256((uint64_t *) ro256, Out);
	
	return 0;
}






