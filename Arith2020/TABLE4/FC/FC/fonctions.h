/*****************************************************************



*****************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <immintrin.h>


#include "ccount.h"
//include HQC :

#include "gf2x_avx2.h"

#ifndef _FONCTIONS_H
#define _FONCTIONS_H

#define CEIL_DIVIDE(a, b)  ((a/b) + (a % b == 0 ? 0 : 1)) /*!< Divide a by b and ceil the result*/

#define WORD 64

	#define SIZE256 CEIL_DIVIDE(PARAM_N,256)
//#endif



#define SIZE64 (SIZE256<<2)
#define SIZE32 (SIZE256<<3)
#define SIZE8 (SIZE256<<5)


#define LAST64 (PARAM_N>>6)




/***************************************************************

	Fonctions d'affichage

*/
  

void afficheVect(uint64_t *A, char *var, int size);


/***************************************************************

	Polynomial Multiplication (in GF2[X])

*/

int fastConvolutionMult(uint64_t * A, uint32_t * B, uint64_t * C, int w);//size est en nombre de mots de 256 bits !!!


#endif
