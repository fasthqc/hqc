/**
 * @file fft.h
 * Header file of fft.c
 */
#ifndef FFT_H
#define FFT_H

#include <stddef.h>
#include <stdint.h>

/**
 * The additive FFT takes a 2^PARAM_FFT polynomial as input.
 * We use the FFT to compute the roots of sigma, whose degree if PARAM_DELTA=60.
 * The smallest power of 2 greater than 60+1 is 64=2^6.
 */
#define PARAM_FFT     6

/**
 * The additive FFT transpose computes a (2^PARAM_FFT_T)-sized syndrome vector.
 * We want to compute 2*PARAM_DELTA=120 syndromes.
 * The smallest power of 2 greater than 120 is 2^7.
 */
#define PARAM_FFT_T   7

void compute_fft_lut(uint16_t* gammas_sums, uint16_t* fft_final_betas,
		     uint16_t* fft_t_final_betas_sums, uint16_t* betas_pows);

void compute_fft_betas(uint16_t* betas);

void compute_subset_sums(uint16_t* subset_sums, uint16_t* set, size_t set_size);

void radix_t_1024_4(uint16_t* f, const uint16_t* f0, const uint16_t* f1);

void radix_t_1024_8(uint16_t* f, const uint16_t* f0, const uint16_t* f1);

void radix_t_1024_16(uint16_t* f, const uint16_t* f0, const uint16_t* f1);

void radix_t_1024_32(uint16_t* f, const uint16_t* f0, const uint16_t* f1);

void radix_t_1024_64(uint16_t* f, const uint16_t* f0, const uint16_t* f1);

void radix_t_1024_128(uint16_t* f, const uint16_t* f0, const uint16_t* f1);

void fft_t_1024_2(uint16_t* f, const uint16_t* w);

void fft_t_1024_4(uint16_t* f, const uint16_t* w, size_t f_coeffs);

void fft_t_1024_8(uint16_t* f, const uint16_t* w, size_t f_coeffs);

void fft_t_1024_16(uint16_t* f, const uint16_t* w, size_t f_coeffs);

void fft_t_1024_32(uint16_t* f, const uint16_t* w, size_t f_coeffs);

void fft_t_1024_64(uint16_t* f, const uint16_t* w, size_t f_coeffs);

void fft_t_1024_128(uint16_t* f, const uint16_t* w, size_t f_coeffs);

void radix_1024_4(uint16_t* f0, uint16_t* f1, const uint16_t* f);

void radix_1024_8(uint16_t* f0, uint16_t* f1, const uint16_t* f);

void radix_1024_16(uint16_t* f0, uint16_t* f1, const uint16_t* f);

void radix_1024_32(uint16_t* f0, uint16_t* f1, const uint16_t* f);

void radix_1024_64(uint16_t* f0, uint16_t* f1, const uint16_t* f);

void fft_1024_2 (uint16_t* w, const uint16_t* f);

void fft_1024_4(uint16_t* w, uint16_t* f, size_t f_coeffs);

void fft_1024_8(uint16_t* w, uint16_t* f, size_t f_coeffs);

void fft_1024_16(uint16_t* w, uint16_t* f, size_t f_coeffs);

void fft_1024_32(uint16_t* w, uint16_t* f, size_t f_coeffs);

void fft_1024_64(uint16_t* w, uint16_t* f, size_t f_coeffs);

void fft_t_preprocess_bch_codeword(uint16_t* w, const uint8_t* rcv);

void fft_retrieve_bch_error_poly(uint8_t* err, const uint16_t* w);

#endif // FFT_H
