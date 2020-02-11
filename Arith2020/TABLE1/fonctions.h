/*****************************************************************



*****************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>


//include HQC :

#include <gf2x.h>



#include <immintrin.h>
#include "ccount.h"

#ifndef _FONCTIONS_H
#define _FONCTIONS_H

#define CEIL_DIVIDE(a, b)  ((a/b) + (a % b == 0 ? 0 : 1)) /*!< Divide a by b and ceil the result*/

#define WORD 64




/***************************************************************

	Fonctions d'affichage

*/
  

void afficheVect(uint64_t *A, char *var, int size);


/***************************************************************

	Polynomial Multiplication (in GF2[X])

*/


int KaratRecPclmul(uint64_t * A, uint64_t * B, uint64_t * C, int size);//size est en nombre de mots de 256 bits !!!


#endif
