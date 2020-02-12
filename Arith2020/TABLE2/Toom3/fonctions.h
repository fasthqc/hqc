/*****************************************************************



*****************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>


//include HQC :

#include <gf2x.h>



#include <immintrin.h>
#include <gmp.h>
#include "ccount.h"

#ifndef _FONCTIONS_H
#define _FONCTIONS_H

#define CEIL_DIVIDE(a, b)  ((a/b) + (a % b == 0 ? 0 : 1)) /*!< Divide a by b and ceil the result*/

#define WORD 64

#define SIZE_N (256*CEIL_DIVIDE(PARAM_N, 256))
#define t (SIZE_N/WORD)



#define T_TM3 SIZE_N
#define tTM3 (T_TM3/WORD)
#define T_TM3_3W (T_TM3/3)
#define T_TM3_3W_256 ((T_TM3_3W+128)/(4*WORD))// on rajoute 128 bits pour avoir un nombre entier de mots de 256 bits.
#define T_TM3_3W_64 (T_TM3_3W_256<<2)




/***************************************************************

	Fonctions d'affichage

*/
  

void afficheVect(uint64_t *A, char *var, int size);


/***************************************************************

	Reverses the bits: MSB -> LSB and vice-versa

*/


void reverse_u32(uint32_t * A, uint64_t * Out, int s);


void reverse_u64(uint64_t * A, uint64_t * Out, int s);

/***************************************************************

	Polynomial Multiplication (in GF2[X])

*/




int Toom3Mult(uint64_t * A, uint64_t * B, uint64_t * C);//size est en nombre de mots de 256 bits !!!


#endif
