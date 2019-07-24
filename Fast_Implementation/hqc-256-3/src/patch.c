#include <immintrin.h>

#define WORD 64
#define SIZE256 ((PARAM_N>>8)+1)
#define SIZE64 (SIZE256<<2)
#define SIZE8 (SIZE256<<5)

#define LAST64 (PARAM_N>>6)


//const int size1 = (PARAM_N>>6)+1;
uint8_t tmp8[VEC_N_SIZE_BYTES];
uint64_t C2[SIZE64<<1];

uint64_t TT[SIZE64];
__m256i * B256 = (__m256i *) TT;


/**
 * \fn void vect_fixed_weight_by_coord(uint8_t* v, const uint16_t weight, AES_XOF_struct* ctx)
 * \brief Generates a vector of a given Hamming weight
 *
 * This function generates uniformly at random a binary vector of a Hamming weight equal to the parameter <b>weight</b>. The vector 
 * is stored by position. 
 * To generate the vector we have to sample uniformly at random values in the interval [0, PARAM_N -1]. Suppose the PARAM_N is equal to \f$ 70853 \f$, to select a position \f$ r\f$ the function works as follow:
 *  1. It makes a call to the seedexpander function to obtain a random number \f$ x\f$ in \f$ [0, 2^{24}[ \f$.
 *  2. Let \f$ t = \lfloor {2^{24} \over 70853} \rfloor \times  70853\f$
 *  3. If \f$ x \geq t\f$, go to 1
 *  4. It return \f$ r = x \mod 70853\f$ 
 *
 * The parameter \f$ t \f$ is precomputed and it's denoted by UTILS_REJECTION_THRESHOLD (see the file parameters.h).
 *
 * \param[in] v Pointer to an array
 * \param[in] weight Integer that is the Hamming weight
 * \param[in] ctx Pointer to the context of the seed expander
 */


void vect_fixed_weight_by_position(uint32_t* v, const uint16_t weight, AES_XOF_struct* ctx) {
  unsigned long random_bytes_size = 3 * weight;
  unsigned char* rand_bytes = (unsigned char*) calloc(random_bytes_size, sizeof(unsigned char));

  seedexpander(ctx, rand_bytes, random_bytes_size);
  
  unsigned long j =0;
  uint32_t random_data;
  uint8_t exist;

  for(uint32_t i = 0; i < weight ; ++i) {
    exist = 0;
    do {
      if(j == random_bytes_size) {
        seedexpander(ctx, rand_bytes, random_bytes_size);
        j = 0;
      }

      random_data  = ((uint32_t) rand_bytes[j++]) << 16;
      random_data |= ((uint32_t) rand_bytes[j++]) << 8;
      random_data |= rand_bytes[j++]; 

    } while(random_data >= UTILS_REJECTION_THRESHOLD);

    random_data = random_data % PARAM_N;

    for(uint32_t k = 0 ; k < i ; k++) {
      if(v[k] == random_data) exist = 1;
    }
    if(exist == 1) {
      i--;
    } else {
      v[i] = random_data;
    }
  }


  free(rand_bytes);

}

void vector_add_by_position_and_coordinate(uint8_t * o, uint32_t * v1, uint8_t * v2, uint32_t W) {
	uint32_t * tmp = (uint32_t *) v2;
	for(uint16_t i = 0 ; i < W ; ++i)	{
		int index = v1[i] >> 5;
		tmp[index] ^= 1U << (v1[i] & 0x1fU);
	}
	
	memcpy(o, tmp, VEC_N_SIZE_BYTES);
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
	__m256i r256, carry256, mask256;
	__m256i * A256 = (__m256i *) A;
	
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
	
	
	//PARAM_N 70853
	r256 = _mm256_lddqu_si256((__m256i const *)(& A[i]));
	r256 = _mm256_srli_epi64(r256,dec64);
	carry256 = _mm256_lddqu_si256((__m256i const *)(& A[i+1]));
	carry256 = _mm256_slli_epi64(carry256,d0);
	r256 ^= carry256;
	i2 = (i-i64)>>2;
	B256[i2] = (A256[i2]^r256);
	TT[LAST64] &= 0x1fUL;

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


int fastConvolutionMult(uint64_t * A, uint32_t * vB, uint64_t * C, int W)//size est en nombre de mots de 64 bits !!!
{
	
	int dec,d0,dec64, i64,i,j,vbi;

	__m256i A256;
	__m256i * C256 = (__m256i *) C2;
	__m256i carry256;
	__m256i zero256 = (__m256i){0x0UL,0x0UL,0x0UL,0x0UL};
	
	
	
	
	for(int j=0;j<SIZE256<<1;j++) C256[j] = zero256;
	
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
	

	for(i=0;i<W;i++){
		dec = vB[i]&0xff;
		dec64 = dec&0x3f;
		i64=(dec>>6);
		d0=WORD-dec64;
		vbi = (vB[i]>>8);
		// prologue de la boucle de décalage
		// pour les shifts vectorisés ( _mm256_slli_epi64() ou _mm256_srli_epi64()), il n'y a pas nécessité de vérifier la taille du décalage < 64.
		
		A256 = tab[i64];
		carry256 = tab[i64+1];
		A256 = _mm256_slli_epi64(A256,dec64);
		carry256 = _mm256_srli_epi64(carry256,d0);
		A256^=carry256;
		C256[vbi] ^= A256;
		
		// boucle vectorisée de décalage
		
		
		for(j=1;j<SIZE256;j++)
		{
			A256 = _mm256_lddqu_si256((__m256i const *)(& A[4*j-i64]));
			A256 = _mm256_slli_epi64(A256,dec64);
			carry256 = _mm256_lddqu_si256((__m256i const *)(& A[4*j-i64-1]));
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
	
	reduction256(C2,C);
	
	
	
	return 0;

}


