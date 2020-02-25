/**
 * \file gf2x_avx2.c
 * \brief AVX2 implementation of multiplication of two polynomials
 */

#include "gf2x_avx2.h"


#define MASK 0x7ffffffffffffUL

	//PARAM_N 24677 : mask = 0x1fffffffffUL;
	//PARAM_N 43669 : mask = 0x1fffffUL;
	//PARAM_N 46747 : mask = 0x7ffffffUL;
	//PARAM_N 63587 : mask = 0x7ffffffffUL;
	//PARAM_N 67699 : mask = 0x7ffffffffffffUL
	//PARAM_N 70853 : mask = 0x1fUL;

uint64_t C2[SIZE64<<1];
uint64_t TT[SIZE64];
	//__m256i * B256 = (__m256i *) TT;


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
	
	//for(int i=0;i<size;i++) nB[i]=0UL;
	
	for(int i=0;i<w;i++)
		nB[vB[i]>>6] ^= 1UL<< (vB[i]%WORD);
}





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



/*************************************************************************************
//
//                       REDUCTION
//
//	Reduction modulo (X^N-1)
//  256 bits
//
*************************************************************************************/

static inline int reduction256(uint64_t * A, uint64_t * B)
{
	__m256i r256, carry256;//, mask256;
	__m256i * A256 = (__m256i *) A;

	__m256i * B256 = (__m256i *) TT;

	
	const int dec64 = PARAM_N&0x3f, i64=PARAM_N>>6, d0=WORD-dec64;

		int i,i2;
	
	for(i=i64;i<(PARAM_N>>5)-4;i+=4)
	{
		
		r256 = _mm256_lddqu_si256((__m256i const *)(& A[i]));
		
		r256 = _mm256_srli_epi64(r256,dec64);
		carry256 = _mm256_lddqu_si256((__m256i const *)(& A[i+1]));
		carry256 = _mm256_slli_epi64(carry256,d0);
		r256 ^= carry256;
		i2 =(i-i64)>>2;
		B256[i2] = A256[i2]^r256;
		
	}
	
	
	//PARAM_N
	/*r256 = (__m256i){A[i],0x0UL,0x0UL,0x0UL};
	r256 = _mm256_srli_epi64(r256,dec64);
	i2 =(i-i64)>>2;
	B256[i2] = (A256[i2]^r256);*/
	
	//PARAM_N  24677 - 63587 - 67699
	r256 = (__m256i){A[i],A[i+1],0x0UL,0x0UL};
	carry256 = _mm256_lddqu_si256((__m256i const *)(& A[i+1]));
	r256 = _mm256_srli_epi64(r256,dec64);
	carry256 = _mm256_slli_epi64(carry256,d0);
	r256 ^= carry256;
	i2 = (i-i64)>>2;
	B256[i2] = (A256[i2]^r256);
	TT[LAST64] &= MASK;
	
	
	
	//PARAM_N 43699 - 47647
	/*r256 = (__m256i){A[i],A[i+1],A[i+2],0x0UL};
	r256 = _mm256_srli_epi64(r256,dec64);
	carry256 = _mm256_lddqu_si256((__m256i const *)(& A[i+1]));
	carry256 = _mm256_slli_epi64(carry256,d0);
	r256 ^= carry256;
	i2 = (i-i64)>>2;
	B256[i2] = (A256[i2]^r256);
	TT[LAST64] &= MASK;*/
	
	
	//PARAM_N 70853
	/*r256 = _mm256_lddqu_si256((__m256i const *)(& A[i]));
	r256 = _mm256_srli_epi64(r256,dec64);
	carry256 = _mm256_lddqu_si256((__m256i const *)(& A[i+1]));
	carry256 = _mm256_slli_epi64(carry256,d0);
	r256 ^= carry256;
	i2 = (i-i64)>>2;
	B256[i2] = (A256[i2]^r256);
	TT[LAST64] &= MASK;*/
	
	//printf("TT[LAST64] = %16.16x\n",TT[LAST64]);

	memcpy(B, TT, VEC_N_SIZE_BYTES);


	return 0;
}



/*************************************************************************************
//
//                       MULTIPLICATION
//
//	avec Reduction modulo (X^N-1)
//  256 bits
//
*************************************************************************************/


int gf2x_Mult(uint64_t * Out, uint64_t * A, uint64_t * B)//size est en nombre de mots de 64 bits !!!
{
	
	gf2x_mul( C2,A,SIZE64,B,SIZE64);

	reduction256(C2,Out);
	//afficheVect(Out,"Out   ",SIZE64/*tTM3R 24*(T_TM3R_3W_256+2)*/);
	
	
	return 0;

}


