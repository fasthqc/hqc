

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

inline int KaratRecPclmul(uint64_t * A64, uint64_t * B64, uint64_t * C64, int size)//size est en nombre de mots de 64 bits !!!
{

	__m256i * A = (__m256i *) A64;
	__m256i * B = (__m256i *) B64;

	__m256i * C = (__m256i *) C64;

	KaratRecPclmul256(A,B,C,size>>2);
		

	return 0;

}


