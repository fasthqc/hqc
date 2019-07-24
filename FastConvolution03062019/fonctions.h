/*****************************************************************



*****************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>



#include <immintrin.h>
#include <gmp.h>
#include "ccount.h"


#ifndef _FONCTIONS_H
#define _FONCTIONS_H

#define WORD 64

#define SIZE_N 49152
#define t (SIZE_N/WORD)

#define W 117


/***************************************************************

	Fonctions d'affichage

*/
  

void afficheVect(uint64_t *A, char *var, int size);


/***************************************************************

	Reverses the bits: MSB -> LSB and vice-versa

*/


void reverse_u32(uint32_t * A, uint64_t * Out, int s);


/***************************************************************

	Polynomial Multiplication (in GF2[X])

*/



int multPclmulSsRed(uint64_t * A, uint64_t * B, uint64_t * C, int size);//size est en nombre de mots de 256 bits !!!

int KaratRecPclmul(uint64_t * A, uint64_t * B, uint64_t * C, int size);//size est en nombre de mots de 256 bits !!!

int fiveWayMult(uint64_t * A, uint64_t * B, uint64_t * C);//size est en nombre de mots de 256 bits !!!
// inline static int KaratRecPclmul128(__m128i * A, __m128i * B, __m128i * C, int size);//size est en nombre de mots de 256 bits !!!
int threeWayMult(uint64_t* A, uint64_t* B, uint64_t* Out);//size est en nombre de mots de 256 bits !!!

/***************************************************************

	Polynomial Multiplication with fast convolution algorithm (in GF2[X])
	
	(with one sparse polynomial as operand)

*/



int fastConvolutionMult(uint64_t * A, int * B, uint64_t * C, int size);//size est en nombre de mots de 256 bits !!!


#endif
