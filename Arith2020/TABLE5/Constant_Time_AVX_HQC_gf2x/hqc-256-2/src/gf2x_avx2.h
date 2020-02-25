/**
 * @file gf2x_avx2.h
 * @brief Header file for gf2x_avx2.c
 */

#ifndef GF2X_AVX2_H
#define GF2X_AVX2_H

#include <stdlib.h>
#include <string.h>
#include <gf2x.h>

#include "parameters.h"

#include <immintrin.h>

#define WORD 64

#define SIZE64 CEIL_DIVIDE(PARAM_N, 64)//pour avoir un nombre entier de mots de 64 bits ((PARAM_N>>6)+1)
#define SIZE8 ((SIZE64<<3))

#define LAST64 (PARAM_N>>6)

/**
 * multiplication gf2x mod x^n-1
 */



void initFWVector(uint64_t * nB,uint32_t * vB, int w);

int gf2x_Mult(uint64_t* Out, uint64_t* A, uint64_t* B);




#endif
