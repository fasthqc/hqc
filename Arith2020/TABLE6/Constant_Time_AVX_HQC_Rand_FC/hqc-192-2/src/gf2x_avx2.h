/**
 * @file gf2x_avx2.h
 * @brief Header file for gf2x_avx2.c
 */

#ifndef GF2X_AVX2_H
#define GF2X_AVX2_H

#include <stdlib.h>
#include <string.h>

#include "parameters.h"
#include "vector.h"

#include <immintrin.h>

#define WORD 64
#define SIZE256 ((PARAM_N>>8)+1)
#define SIZE64 (SIZE256<<2)
#define SIZE8 (SIZE256<<5)

#define LAST64 (PARAM_N>>6)

/**
 * fast convolution avx2
 */

int fastConvolutionMult(uint64_t * A, uint32_t * vB, uint64_t * C, int W);//size est en nombre de mots de 64 bits !!!

#endif
