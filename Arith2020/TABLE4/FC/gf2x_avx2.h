/**
 * \file gf2x_avx2.h
 * \brief Header file for gf2x_avx2.c
 */

#ifndef GF2X_AVX2_H
#define GF2X_AVX2_H

#include <stdlib.h>
#include <string.h>



//#include "vector.h"
#define CEIL_DIVIDE(a, b)  ((a/b) + (a % b == 0 ? 0 : 1)) /*!< Divide a by b and ceil the result*/
#define BITMASK(a, size) ((1UL << (a % size)) - 1) /*!< Create a mask*/

#define VECTORSIZE 256 /*!< The size in bits of a vector*/
#define BYTESPERVECTOR (VECTORSIZE/8) /*!< The number of bytes in a vector*/
#define BITSPERINT 32 /*!< The number of bits per int*/
#define INTSPERVECTOR (VECTORSIZE/BITSPERINT) /*!< The number of ints per vector*/

//#include "parameters.h"

#define VEC_N_SIZE_BYTES                    CEIL_DIVIDE(PARAM_N, 8)
#define VEC_N1_SIZE_BYTES                   CEIL_DIVIDE(PARAM_N1, 8)
#define VEC_N1N2_SIZE_BYTES                 CEIL_DIVIDE(PARAM_N1N2, 8)
#define VEC_K_SIZE_BYTES                    CEIL_DIVIDE(PARAM_K, 8)

#define VEC_N_ARRAY_SIZE_VEC CEIL_DIVIDE(PARAM_N,256) /*!< The number of needed vectors to store PARAM_N bits*/



/**
 * v8uint32: a vector of 8 unsigned ints
 */
typedef uint32_t v8uint32 __attribute__ ((vector_size (BYTESPERVECTOR)))
                        __attribute__ ((aligned (BYTESPERVECTOR)));

void shiftXor(v8uint32 dest[], v8uint32 source[], unsigned int length, unsigned int distance);

void vect_sparse_from_bytes(uint32_t* o, const uint8_t* a);

void sparse_dense_mul(uint8_t* o, const uint32_t* v1, const uint8_t* v2, const uint32_t weight);

#endif
