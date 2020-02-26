/*****************************************************************



*****************************************************************/



#include "fonctions.h"


uint64_t C2[SIZE64<<1];

uint64_t TT[SIZE64];
__m256i * B256 = (__m256i *) TT;





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


/*
	procédure d'échange dans un tableau entre deux éléments, pour P_omega
*/

static inline void swap(int * tab, int elt1, int elt2)
{
	int sw = tab[elt1];
	
	tab[elt1] = tab[elt2];
	tab[elt2] = sw;
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
	//__m256i * B256 = (__m256i *) B;
	
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
	
	
	
	//PARAM_N  24677 - 63587 - 67699
	
	#ifdef P24677
	r256 = (__m256i){A[i],A[i+1],0x0UL,0x0UL};
	carry256 = _mm256_lddqu_si256((__m256i const *)(& A[i+1]));// JMR improved Oleg's line
	r256 = _mm256_srli_epi64(r256,dec64);
	carry256 = _mm256_slli_epi64(carry256,d0);
	r256 ^= carry256;
	i2 = (i-i64)>>2;
	B256[i2] = (A256[i2]^r256);//&mask256;
	TT[LAST64] &= 0x1fffffffffUL;//63587 = 0x7ffffffffUL ;67699 = 0x7ffffffffffffUL ; 24677 = 0x1fffffffffUL*/
	#endif

	#ifdef P63587
	r256 = (__m256i){A[i],A[i+1],0x0UL,0x0UL};
	carry256 = _mm256_lddqu_si256((__m256i const *)(& A[i+1]));// JMR improved Oleg's line
	r256 = _mm256_srli_epi64(r256,dec64);
	carry256 = _mm256_slli_epi64(carry256,d0);
	r256 ^= carry256;
	i2 = (i-i64)>>2;
	B256[i2] = (A256[i2]^r256);//&mask256;
	TT[LAST64] &= 0x7ffffffffUL;//63587 = 0x7ffffffffUL ;67699 = 0x7ffffffffffffUL ; 24677 = 0x1fffffffffUL*/
	#endif
	
	#ifdef P67699
	r256 = (__m256i){A[i],A[i+1],0x0UL,0x0UL};
	carry256 = _mm256_lddqu_si256((__m256i const *)(& A[i+1]));// JMR improved Oleg's line
	r256 = _mm256_srli_epi64(r256,dec64);
	carry256 = _mm256_slli_epi64(carry256,d0);
	r256 ^= carry256;
	i2 = (i-i64)>>2;
	B256[i2] = (A256[i2]^r256);//&mask256;
	TT[LAST64] &= 0x7ffffffffffffUL;//63587 = 0x7ffffffffUL ;67699 = 0x7ffffffffffffUL ; 24677 = 0x1fffffffffUL*/
	#endif
	
	
	//PARAM_N 43669 - 46747
	#ifdef P43669
	r256 = (__m256i){A[i],A[i+1],A[i+2],0x0UL};
	r256 = _mm256_srli_epi64(r256,dec64);
	carry256 = _mm256_lddqu_si256((__m256i const *)(& A[i+1]));
	carry256 = _mm256_slli_epi64(carry256,d0);
	r256 ^= carry256;
	i2 = (i-i64)>>2;
	B256[i2] = (A256[i2]^r256);//&mask256;
	TT[LAST64] &= 0x1fffffUL;//47647 - 0x7ffffffUL;//43669 - 0x1fffffUL;//*/
	#endif
	
	
	#ifdef P46747
	r256 = (__m256i){A[i],A[i+1],A[i+2],0x0UL};
	r256 = _mm256_srli_epi64(r256,dec64);
	carry256 = _mm256_lddqu_si256((__m256i const *)(& A[i+1]));
	carry256 = _mm256_slli_epi64(carry256,d0);
	r256 ^= carry256;
	i2 = (i-i64)>>2;
	B256[i2] = (A256[i2]^r256);//&mask256;
	TT[LAST64] &= 0x7ffffffUL;//47647 - 0x7ffffffUL;//43669 - 0x1fffffUL;//*/
	#endif
	
	
	//PARAM_N 70853
	#ifdef P70853
	r256 = _mm256_lddqu_si256((__m256i const *)(& A[i]));
	r256 = _mm256_srli_epi64(r256,dec64);
	carry256 = _mm256_lddqu_si256((__m256i const *)(& A[i+1]));
	carry256 = _mm256_slli_epi64(carry256,d0);
	r256 ^= carry256;
	i2 = (i-i64)>>2;
	B256[i2] = (A256[i2]^r256);//&mask256;
	TT[LAST64] &= 0x1fUL;//*/
	#endif
	

	memcpy(B, TT, VEC_N_SIZE_BYTES);




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


static inline int fastConvolutionMultWithoutRed(uint64_t * A, int * vB, uint64_t * C, int w)//size est en nombre de mots de 64 bits !!!Origin
{
	int dec,d0,dec64,vbi,i64;
	
	__m256i A256;
	__m256i * C256 = (__m256i *) C;
	__m256i carry256;
	__m256i zero256 = (__m256i){0x0UL,0x0UL,0x0UL,0x0UL};
	
	int P_omega[PARAM_OMEGA];
	
	for(int i=0;i<w;i++) P_omega[i]=i;
	
	for(int i=0;i<w-1;i++) swap(P_omega+i,0,rand()%(w-i));
	
	for(int j=0;j<SIZE256<<1;j++) C256[j] = zero256;
	
	
	//pour la troisième version sans switch : tableau !
	
	__m256i tab[5], tabprime[5];
	tab[0] = (__m256i){A[0],A[1],A[2],A[3]};
	tab[1] = (__m256i){0x0UL,A[0],A[1],A[2]};
	tab[2] = (__m256i){0x0UL,0x0UL,A[0],A[1]};
	tab[3] = (__m256i){0x0UL,0x0UL,0x0UL,A[0]};
	tab[4] = (__m256i){0x0UL,0x0UL,0x0UL,0x0UL};
	
	tabprime[0] = (__m256i){0x0UL,0x0UL,0x0UL,0x0UL};
	tabprime[1] = (__m256i){A[SIZE64-1],0x0UL,0x0UL,0x0UL};
	tabprime[2] = (__m256i){A[SIZE64-2],A[SIZE64-1],0x0UL,0x0UL};
	tabprime[3] = (__m256i){A[SIZE64-3],A[SIZE64-2],A[SIZE64-1],0x0UL};
	tabprime[4] = (__m256i){A[SIZE64-4],A[SIZE64-3],A[SIZE64-2],A[SIZE64-1]};
	

	for(int i=0;i<w;i++){
		dec = vB[P_omega[i]]&0xff;
		dec64 = dec&0x3f;
		i64=(dec>>6);
		d0=WORD-dec64;
		vbi = (vB[P_omega[i]]>>8);
		// prologue de la boucle de décalage
		// pour les shifts vectorisés ( _mm256_slli_epi64() ou _mm256_srli_epi64()), il n'y a pas nécessité de vérifier la taille du décalage < 64.
		
		A256 = tab[i64];
		carry256 = tab[i64+1];
		A256 = _mm256_slli_epi64(A256,dec64);
		carry256 = _mm256_srli_epi64(carry256,d0);
		A256^=carry256;
		C256[vbi] ^= A256;
		
		// boucle vectorisée de décalage
		
		
		for(int j=1;j<SIZE256;j++)
		{
			int k = (j<<2)-i64;
			A256 = _mm256_lddqu_si256((__m256i const *)(& A[k]));
			A256 = _mm256_slli_epi64(A256,dec64);
			carry256 = _mm256_lddqu_si256((__m256i const *)(& A[--k]));
			carry256 = _mm256_srli_epi64(carry256,d0);
			A256^=carry256;
			C256[j+vbi] ^= A256;
		}	
	
		
		// épilogue de la boucle de décalage
		// pour les shifts vectorisés ( _mm256_slli_epi64() ou _mm256_srli_epi64()), il n'y a pas nécessité de vérifier la taille du décalage < 64.
		
		A256 = tabprime[i64];
		carry256 = tabprime[i64+1];
		A256 = _mm256_slli_epi64(A256,dec64);
		carry256 = _mm256_srli_epi64(carry256,d0);
		A256^=carry256;
		C256[SIZE256+vbi] ^= A256;
		
	}
	return 0;

}


// Wrapper fastConvolution et réduction
int fastConvolutionMult(uint64_t * A, uint32_t * vB, uint64_t * C, int w)
{
	fastConvolutionMultWithoutRed(A, vB, C2, w);
	
	reduction256(C2,C);
	
	
	return 0;

}


#include "gf2x_avx2.c"

