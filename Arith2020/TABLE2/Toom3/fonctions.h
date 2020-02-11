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




#define t (SIZE_N/WORD)
#define T_3W (SIZE_N/3)
#define T_3W_256 (T_3W/256)
#define T2_3W_256 (2*T_3W_256)

#define TOOM64
//constantes TOOM3 avec décalage de 64 bits -> SIZE_N doit être divisible par 3 et par 64, soit 192 !!!!!!!
// et T_TM3_3W_256 correspond à une puissance de 2 (on appelle KaratRecMult !) :
// d'où SIZE_N = 2688, 5760, 11904, 24192 ou 48768

#define T_TM3 SIZE_N
#define tTM3 (T_TM3/WORD)
#define T_TM3_3W (T_TM3/3)
#define T_TM3_3W_256 ((T_TM3_3W+128)/(4*WORD))// on rajoute 128 bits pour avoir un nombre entier de mots de 256 bits.
#define T_TM3_3W_64 (T_TM3_3W_256<<2)


//constantes TOOM3REC

/*#define T_TM3R 67584
#define tTM3R (T_TM3R/WORD)
#define T_TM3R_3W (T_TM3R/3)
#define T_TM3R_3W_256 (T_TM3R_3W/(4*WORD))
#define T_TM3R_3W_64 (T_TM3R_3W/WORD)*/

#define T_TM3R_3W (PARAM_N-512)
#define T_TM3R (3*(T_TM3R_3W+128))
#define tTM3R (T_TM3R/WORD)
#define T_TM3R_3W_256 ((T_TM3R_3W+128)/(4*WORD))
#define T_TM3R_3W_64 (T_TM3R_3W_256<<2)




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
