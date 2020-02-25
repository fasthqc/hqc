/**
 * @file fft.c
 * Implementation of the additive FFT and its transpose.
 * This implementation is based on the paper from Gao and Mateer: <br>
 * Shuhong Gao and Todd Mateer, Additive Fast Fourier Transforms over Finite Fields,
 * IEEE Transactions on Information Theory 56 (2010), 6265--6272.
 * http://www.math.clemson.edu/~sgao/papers/GM10.pdf <br>
 * and includes improvements proposed by Bernstein, Chou and Schwabe here:
 * https://binary.cr.yp.to/mcbits-20130616.pdf
 */
#include "fft.h"

#include "gf.h"
#include "parameters.h"

#include <string.h>

/**
 * Constants for the additive FFT and its transpose.
 */
const uint16_t gammas_sums[512+256+128+64+32+16] = {
  0,512,256,768,128,640,384,896,64,576,320,832,192,704,448,960,32,544,288,800,160,672,416,928,96,608,352,864,224,736,480,992,16,528,272,784,144,656,400,912,80,592,336,848,208,720,464,976,48,560,304,816,176,688,432,944,112,624,368,880,240,752,496,1008,8,520,264,776,136,648,392,904,72,584,328,840,200,712,456,968,40,552,296,808,168,680,424,936,104,616,360,872,232,744,488,1000,24,536,280,792,152,664,408,920,88,600,344,856,216,728,472,984,56,568,312,824,184,696,440,952,120,632,376,888,248,760,504,1016,2,514,258,770,130,642,386,898,66,578,322,834,194,706,450,962,34,546,290,802,162,674,418,930,98,610,354,866,226,738,482,994,18,530,274,786,146,658,402,914,82,594,338,850,210,722,466,978,50,562,306,818,178,690,434,946,114,626,370,882,242,754,498,1010,10,522,266,778,138,650,394,906,74,586,330,842,202,714,458,970,42,554,298,810,170,682,426,938,106,618,362,874,234,746,490,1002,26,538,282,794,154,666,410,922,90,602,346,858,218,730,474,986,58,570,314,826,186,698,442,954,122,634,378,890,250,762,506,1018,237,749,493,1005,109,621,365,877,173,685,429,941,45,557,301,813,205,717,461,973,77,589,333,845,141,653,397,909,13,525,269,781,253,765,509,1021,125,637,381,893,189,701,445,957,61,573,317,829,221,733,477,989,93,605,349,861,157,669,413,925,29,541,285,797,229,741,485,997,101,613,357,869,165,677,421,933,37,549,293,805,197,709,453,965,69,581,325,837,133,645,389,901,5,517,261,773,245,757,501,1013,117,629,373,885,181,693,437,949,53,565,309,821,213,725,469,981,85,597,341,853,149,661,405,917,21,533,277,789,239,751,495,1007,111,623,367,879,175,687,431,943,47,559,303,815,207,719,463,975,79,591,335,847,143,655,399,911,15,527,271,783,255,767,511,1023,127,639,383,895,191,703,447,959,63,575,319,831,223,735,479,991,95,607,351,863,159,671,415,927,31,543,287,799,231,743,487,999,103,615,359,871,167,679,423,935,39,551,295,807,199,711,455,967,71,583,327,839,135,647,391,903,7,519,263,775,247,759,503,1015,119,631,375,887,183,695,439,951,55,567,311,823,215,727,471,983,87,599,343,855,151,663,407,919,23,535,279,791,
  0,786,832,82,16,770,848,66,100,886,804,54,116,870,820,38,41,827,873,123,57,811,889,107,77,863,781,31,93,847,797,15,272,514,592,322,256,530,576,338,372,614,564,294,356,630,548,310,313,555,633,363,297,571,617,379,349,591,541,271,333,607,525,287,72,858,776,26,88,842,792,10,44,830,876,126,60,814,892,110,97,883,801,51,113,867,817,35,5,791,837,87,21,775,853,71,344,586,536,266,328,602,520,282,316,558,636,366,300,574,620,382,369,611,561,291,353,627,545,307,277,519,597,327,261,535,581,343,6,788,838,84,22,772,854,68,98,880,802,48,114,864,818,32,47,829,879,125,63,813,895,109,75,857,779,25,91,841,795,9,278,516,598,324,262,532,582,340,370,608,562,288,354,624,546,304,319,557,639,365,303,573,623,381,347,585,539,265,331,601,523,281,78,860,782,28,94,844,798,12,42,824,874,120,58,808,890,104,103,885,807,53,119,869,823,37,3,785,835,81,19,769,851,65,350,588,542,268,334,604,526,284,314,552,634,360,298,568,618,376,375,613,567,293,359,629,551,309,275,513,595,321,259,529,579,337,
  0,304,3,307,717,1021,718,1022,362,90,361,89,935,663,932,660,657,929,658,930,92,364,95,367,1019,715,1016,712,310,6,309,5,620,860,623,863,161,401,162,402,774,566,773,565,459,251,456,248,253,461,254,462,560,768,563,771,407,167,404,164,858,618,857,617,292,20,295,23,1001,729,1002,730,78,382,77,381,643,947,640,944,949,645,950,646,376,72,379,75,735,1007,732,1004,18,290,17,289,840,632,843,635,389,181,390,182,546,786,545,785,239,479,236,476,473,233,474,234,788,548,791,551,179,387,176,384,638,846,637,845,
  0,807,1010,213,318,537,716,491,746,461,280,575,980,243,38,769,831,24,205,1002,513,294,499,724,469,754,551,256,235,972,793,62,887,80,133,930,585,366,443,668,413,698,623,328,163,900,849,118,72,879,954,157,374,593,644,419,674,389,336,631,924,187,110,841,
  0,355,74,297,684,975,742,901,49,338,123,280,669,1022,727,948,644,999,718,941,40,331,98,257,693,982,767,924,25,378,83,304,
  0,819,997,214,804,23,193,1010,246,965,787,32,978,225,55,772 };

/**
 * Constants for the additive FFT and its transpose.
 */
const uint16_t betas_pows[64+32+16+8+4] = {
  1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
  1,18,260,620,592,488,367,3,54,780,701,761,568,945,5,90,285,974,770,577,218,15,238,807,91,271,714,366,17,306,352,237,
  1,893,558,500,335,606,1005,12,561,510,98,863,958,248,80,378,
  1,285,785,353,595,1021,620,577,
  1,790,582,919 };

/**
 * Constants for the additive FFT transpose.
 */
const uint16_t fft_t_final_betas_sums[16] = {
  0,365,27,374,111,258,116,281,351,50,324,41,304,93,299,70 };

/**
 * Constants for the additive FFT.
 */
const uint16_t fft_final_betas[5] = {
  779,42,887,313,790 };

/**
 * Computes all constants involved in the additive FFT and its transpose.
 */
void compute_fft_lut(uint16_t* gammas_sums, uint16_t* fft_final_betas,
		     uint16_t* fft_t_final_betas_sums, uint16_t* betas_pows) {
  uint16_t gammas[PARAM_M-1];
  compute_fft_betas(gammas);
  
  // Compute the first gammas subset sums.
  // Since the PARAM_M-th beta is 1, there's no twisting on the first step,
  // i.e. gammas are betas.
  compute_subset_sums(gammas_sums, gammas, PARAM_M-1);

  gammas_sums += 1 << (PARAM_M-1);
  
  // First call is done. Compute constants for the recursive calls
  for(size_t j = 1 ; j < PARAM_FFT ; ++j) {
    
    // Compute deltas (which are the next betas)
    for(size_t k = 0 ; k < PARAM_M-j ; ++k)
      gammas[k] = gf_square(gammas[k]) ^ gammas[k];

    // Compute last beta powers
    uint16_t beta_m = gammas[PARAM_M-j-1];
    betas_pows[0] = 1;
    for(size_t i = 1 ; i < (1U << (PARAM_FFT_T-j)) ; ++i)
      betas_pows[i] = gf_mul(betas_pows[i-1], beta_m);

    if(j == PARAM_FFT-1) // Compute betas from the last fft call
      memcpy(fft_final_betas, gammas, 2 * (PARAM_M - (PARAM_FFT-1)));
    
    // Compute gammas
    for(size_t k = 0 ; k < PARAM_M-j-1 ; ++k)
      gammas[k] = gf_mul(gammas[k], gf_inverse(beta_m));

    // Compute gammas subset sums
    compute_subset_sums(gammas_sums, gammas, PARAM_M-j-1);

    betas_pows += 1 << (PARAM_FFT_T-j);
    gammas_sums += 1 << (PARAM_M-1-j);
  }

  // Compute betas sums from the last fft transpose call
  fft_t_final_betas_sums[0] = 0;
  for(size_t k = 0 ; k < PARAM_M-PARAM_FFT ; ++k) {
    gammas[k] = gf_square(gammas[k]) ^ gammas[k];
    for(size_t l = 0 ; l < (1U << k) ; ++l)
      fft_t_final_betas_sums[(1 << k) + l] = fft_t_final_betas_sums[l] ^ gammas[k];
  }
}

/**
 * Computes the basis of betas (omitting 1) used in the additive FFT and its transpose.
 * @param[out] betas Array of size PARAM_M-1
 */
void compute_fft_betas(uint16_t* betas) {
  for(size_t i = 0 ; i < PARAM_M-1 ; ++i)
    betas[i] = 1 << (PARAM_M-1-i);

  // Optimisation suggested by Bernstein et al. to avoid twisting.
  // Implemented for GF(2^10) only.
  betas[7] = 2;
  betas[8] = 237;
}

/**
 * Computes the subset sums of the given set.
 * The array subset_sums is such that its ith element is
 * the subset sum of the set elements given by the binary form of i.
 * @param[out] subset_sums Array of size 2^set_size receiving the subset sums
 * @param[in] set Array of set_size elements
 * @param[in] set_size Size of the array set
 */
void compute_subset_sums(uint16_t* subset_sums, uint16_t* set, size_t set_size) {
  subset_sums[0] = 0;
  for(size_t i = 0 ; i < set_size ; ++i)
    for(size_t j = 0 ; j < 1U<<i ; ++j)
      subset_sums[(1 << i) + j] = set[i] ^ subset_sums[j];        
}

/**
 * Last recursive call of the radix conversion transpose on GF(2^10)
 * with polynomials f0 and f1 having 64 coefficients.
 * @param[out] f Array of 4 elements
 * @param[in] f0 Array of 2 elements
 * @param[in] f1 Array of 2 elements
 */
void radix_t_1024_4(uint16_t* f, const uint16_t* f0, const uint16_t* f1) {
  f[0] = f0[0];
  f[1] = f1[0];
  f[2] = f0[1] ^ f1[0];
  f[3] = f[2] ^ f1[1];
}

/**
 * Fourth recursive call of the radix conversion transpose on GF(2^10)
 * with polynomials f0 and f1 having 64 coefficients.
 * @param[out] f Array of 8 elements
 * @param[in] f0 Array of 4 elements
 * @param[in] f1 Array of 4 elements
 */
void radix_t_1024_8(uint16_t* f, const uint16_t* f0, const uint16_t* f1) {
  f[0] = f0[0];
  f[1] = f1[0];
  f[2] = f0[1] ^ f1[0];
  f[3] = f[2] ^ f1[1];
  f[4] = f[2] ^ f0[2];
  f[5] = f[3] ^ f1[2];
  f[6] = f[4] ^ f0[3] ^ f1[2];
  f[7] = f[3] ^ f0[3] ^ f1[3];
}

/**
 * Third recursive call of the radix conversion transpose on GF(2^10)
 * with polynomials f0 and f1 having 64 coefficients.
 * @param[out] f Array of 16 elements
 * @param[in] f0 Array of 8 elements
 * @param[in] f1 Array of 8 elements
 */
void radix_t_1024_16(uint16_t* f, const uint16_t* f0, const uint16_t* f1) {
  f[0] = f0[0];
  f[1] = f1[0];
  f[2] = f0[1] ^ f1[0];
  f[3] = f[2] ^ f1[1];
  f[4] = f[2] ^ f0[2];
  f[5] = f[3] ^ f1[2];
  f[6] = f[4] ^ f0[3] ^ f1[2];
  f[7] = f[3] ^ f0[3] ^ f1[3];
  f[8] = f[4] ^ f0[4];
  f[9] = f[5] ^ f1[4];
  f[10] = f[6] ^ f0[5] ^ f1[4];
  f[11] = f[7] ^ f0[5] ^ f1[4] ^ f1[5];
  f[12] = f[8] ^ f0[5] ^ f0[6] ^ f1[4];
  f[13] = f[7] ^ f[9] ^ f[11] ^ f1[6];
  f[14] = f[6] ^ f0[6] ^ f0[7] ^ f1[6];
  f[15] = f[7] ^ f0[7] ^ f1[7];
}

/**
 * Second recursive call of the radix conversion transpose on GF(2^10)
 * with polynomials f0 and f1 having 64 coefficients.
 * @param[out] f Array of 32 elements
 * @param[in] f0 Array of 16 elements
 * @param[in] f1 Array of 16 elements
 */
void radix_t_1024_32(uint16_t* f, const uint16_t* f0, const uint16_t* f1) {
  uint16_t Q0[8];
  uint16_t Q1[8];
  uint16_t R0[8];
  uint16_t R1[8];
        
  memcpy(Q0, f0 + 8, 2*8);
  memcpy(Q1, f1 + 8, 2*8);
  memcpy(R0, f0, 2*8);
  memcpy(R1, f1, 2*8);
    
  uint16_t Q[2*8];
  uint16_t R[2*8];
        
  radix_t_1024_16(Q, Q0, Q1);
  radix_t_1024_16(R, R0, R1);

  memcpy(f, R, 4*8);
  memcpy(f + 3*8, Q + 8, 2*8);

  for(size_t i = 0 ; i < 8 ; ++i) {
    f[2*8 + i] = Q[i] ^ R[8 + i];
    f[3*8 + i] ^= f[2*8 + i];
  }
}

/**
 * First recursive call of the radix conversion transpose on GF(2^10)
 * with polynomials f0 and f1 having 64 coefficients.
 * @param[out] f Array of 64 elements
 * @param[in] f0 Array of 32 elements
 * @param[in] f1 Array of 32 elements
 */
void radix_t_1024_64(uint16_t* f, const uint16_t* f0, const uint16_t* f1) {
  uint16_t Q0[16];
  uint16_t Q1[16];
  uint16_t R0[16];
  uint16_t R1[16];
        
  memcpy(Q0, f0 + 16, 2*16);
  memcpy(Q1, f1 + 16, 2*16);
  memcpy(R0, f0, 2*16);
  memcpy(R1, f1, 2*16);
    
  uint16_t Q[2*16];
  uint16_t R[2*16];
        
  radix_t_1024_32(Q, Q0, Q1);
  radix_t_1024_32(R, R0, R1);

  memcpy(f, R, 4*16);
  memcpy(f + 3*16, Q + 16, 2*16);

  for(size_t i = 0 ; i < 16 ; ++i) {
    f[2*16 + i] = Q[i] ^ R[16 + i];
    f[3*16 + i] ^= f[2*16 + i];
  }
}

/**
 * Transpose of the radix conversion on GF(2^10)
 * with polynomials f0 and f1 having 64 coefficients.
 * @param[out] f Array of 128 elements
 * @param[in] f0 Array of 64 elements
 * @param[in] f1 Array of 64 elements
 */
void radix_t_1024_128(uint16_t* f, const uint16_t* f0, const uint16_t* f1) {
  uint16_t Q0[32];
  uint16_t Q1[32];
  uint16_t R0[32];
  uint16_t R1[32];
        
  memcpy(Q0, f0 + 32, 2*32);
  memcpy(Q1, f1 + 32, 2*32);
  memcpy(R0, f0, 2*32);
  memcpy(R1, f1, 2*32);
    
  uint16_t Q[2*32];
  uint16_t R[2*32];
        
  radix_t_1024_64(Q, Q0, Q1);
  radix_t_1024_64(R, R0, R1);

  memcpy(f, R, 4*32);
  memcpy(f + 3*32, Q + 32, 2*32);

  for(size_t i = 0 ; i < 32 ; ++i) {
    f[2*32 + i] = Q[i] ^ R[32 + i];
    f[3*32 + i] ^= f[2*32 + i];
  }
}

/**
 * Last recursive call of the additive FFT transpose on GF(2^10)
 * for an array of 1024 elements of GF(2^10).
 * @param[out] f Array of 2 elements
 * @param[in] w Array of 16 elements
 */
void fft_t_1024_2(uint16_t* f, const uint16_t* w) {
  f[0] = 0;
  for(size_t i = 0 ; i < 16 ; ++i)
    f[0] ^= w[i];
  f[1] = 0;
  for(size_t j = 0 ; j < 4 ; ++j) {
    for(size_t k = 0 ; k < (1U << j) ; ++k) {
      size_t index = (1 << j) + k;
      f[1] ^= gf_mul(fft_t_final_betas_sums[index], w[index]);
    }
  }
}

/**
 * Fifth recursive call of the additive FFT transpose on GF(2^10)
 * for an array of 1024 elements of GF(2^10).
 * @param[out] f Array of 4 elements
 * @param[in] w Array of 32 elements
 * @param[in] f_coeffs Length of f
 */
void fft_t_1024_4(uint16_t* f, const uint16_t* w, size_t f_coeffs) {
  uint16_t u[16];
  uint16_t f0[2];
  uint16_t f1[2];

  if(f_coeffs <= 3) {
    f1[1] = 0;
    u[0] = w[0] ^ w[16];
    f1[0] = w[16];
    for(size_t i = 1 ; i < 16 ; ++i) {
      u[i] = w[i] ^ w[16+i];
      f1[0] ^= gf_mul(gammas_sums[512+256+128+64+32+i], u[i]) ^ w[16+i];
    }
    fft_t_1024_2(f0, u);
  }
  else {
    uint16_t v[16];
  
    u[0] = w[0] ^ w[16];
    v[0] = w[16];

    for(size_t i = 1 ; i < 16 ; ++i) {
      u[i] = w[i] ^ w[16+i];
      v[i] = gf_mul(gammas_sums[512+256+128+64+32+i], u[i]) ^ w[16+i];
    }

    fft_t_1024_2(f0, u);
    fft_t_1024_2(f1, v);
  }

  radix_t_1024_4(f, f0, f1);

  for(size_t i = 1 ; i < 4 ; ++i)
    f[i] = gf_mul(betas_pows[64+32+16+8+i], f[i]);
}

/**
 * Fourth recursive call of the additive FFT transpose on GF(2^10)
 * for an array of 1024 elements of GF(2^10).
 * @param[out] f Array of 8 elements
 * @param[in] w Array of 64 elements
 * @param[in] f_coeffs Length of f
 */
void fft_t_1024_8(uint16_t* f, const uint16_t* w, size_t f_coeffs) {
  uint16_t u[32];
  uint16_t v[32];
  
  u[0] = w[0] ^ w[32];
  v[0] = w[32];

  for(size_t i = 1 ; i < 32 ; ++i) {
    u[i] = w[i] ^ w[32+i];
    v[i] = gf_mul(gammas_sums[512+256+128+64+i], u[i]) ^ w[32+i];
  }

  uint16_t f0[4];
  uint16_t f1[4];

  fft_t_1024_4(f0, u, (f_coeffs+1)/2);
  fft_t_1024_4(f1, v, f_coeffs/2);

  radix_t_1024_8(f, f0, f1);

  for(size_t i = 1 ; i < 8 ; ++i)
    f[i] = gf_mul(betas_pows[64+32+16+i], f[i]);
}

/**
 * Third recursive call of the additive FFT transpose on GF(2^10)
 * for an array of 1024 elements of GF(2^10).
 * @param[out] f Array of 16 elements
 * @param[in] w Array of 128 elements
 * @param[in] f_coeffs Length of f
 */
void fft_t_1024_16(uint16_t* f, const uint16_t* w, size_t f_coeffs) {
  uint16_t u[64];
  uint16_t v[64];
  
  u[0] = w[0] ^ w[64];
  v[0] = w[64];

  for(size_t i = 1 ; i < 64 ; ++i) {
    u[i] = w[i] ^ w[64+i];
    v[i] = gf_mul(gammas_sums[512+256+128+i], u[i]) ^ w[64+i];
  }

   uint16_t f0[8];
  uint16_t f1[8];

  fft_t_1024_8(f0, u, (f_coeffs+1)/2);
  fft_t_1024_8(f1, v, f_coeffs/2);

  radix_t_1024_16(f, f0, f1);

  for(size_t i = 1 ; i < 16 ; ++i)
    f[i] = gf_mul(betas_pows[64+32+i], f[i]);
}

/**
 * Second recursive call of the additive FFT transpose on GF(2^10)
 * for an array of 1024 elements of GF(2^10).
 * @param[out] f Array of 32 elements
 * @param[in] w Array of 256 elements
 * @param[in] f_coeffs Length of f
 */
void fft_t_1024_32(uint16_t* f, const uint16_t* w, size_t f_coeffs) {
  uint16_t u[128];
  uint16_t v[128];
  
  u[0] = w[0] ^ w[128];
  v[0] = w[128];

  for(size_t i = 1 ; i < 128 ; ++i) {
    u[i] = w[i] ^ w[128+i];
    v[i] = gf_mul(gammas_sums[512+256+i], u[i]) ^ w[128+i];
  }

  uint16_t f0[16];
  uint16_t f1[16];

  fft_t_1024_16(f0, u, (f_coeffs+1)/2);
  fft_t_1024_16(f1, v, f_coeffs/2);

  radix_t_1024_32(f, f0, f1);

  for(size_t i = 1 ; i < 32 ; ++i)
    f[i] = gf_mul(betas_pows[64+i], f[i]);
}

/**
 * First recursive call of the additive FFT transpose on GF(2^10)
 * for an array of 1024 elements of GF(2^10).
 * @param[out] f Array of 64 elements
 * @param[in] w Array of 512 elements
 * @param[in] f_coeffs Length of f
 */
void fft_t_1024_64(uint16_t* f, const uint16_t* w, size_t f_coeffs) {
  uint16_t u[256];
  uint16_t v[256];
  
  u[0] = w[0] ^ w[256];
  v[0] = w[256];

  for(size_t i = 1 ; i < 256 ; ++i) {
    u[i] = w[i] ^ w[256+i];
    v[i] = gf_mul(gammas_sums[512+i], u[i]) ^ w[256+i];
  }

  /* Step 5: Compute f0 from u and f1 from v */
  uint16_t f0[32];
  uint16_t f1[32];
        
  fft_t_1024_32(f0, u, (f_coeffs+1)/2);
  fft_t_1024_32(f1, v, f_coeffs/2);

  /* Step 3: Compute g from g0 and g1 */
  radix_t_1024_64(f, f0, f1);
}

/**
 * Additive FFT Transpose for a family of 1024 elements of GF(2^10).
 * Computes the first f_coeffs syndromes of family w.
 * @param[out] f Array of 128 elements receiving the syndromes
 * @param[in] w Array of 1024 elements of GF(2^10) storing the family
 * @param[in] f_coeffs Length of f
 */
void fft_t_1024_128(uint16_t* f, const uint16_t* w, size_t f_coeffs) {
  uint16_t u[512];
  uint16_t v[512];

  u[0] = w[0] ^ w[512];
  v[0] = w[512];

  for(size_t i = 1 ; i < 512 ; ++i) {
    u[i] = w[i] ^ w[512+i];
    v[i] = gf_mul(gammas_sums[i], u[i]) ^ w[512+i];
  }

  /* Step 5: Compute f0 from u and f1 from v */
  uint16_t f0[64];
  uint16_t f1[64];
        
  fft_t_1024_64(f0, u, (f_coeffs+1)/ 2);
  fft_t_1024_64(f1, v, f_coeffs/2);

  /* Step 3: Compute g from g0 and g1 */
  radix_t_1024_128(f, f0, f1);
}

/**
 * Last recursive call of the radix conversion on GF(2^10)[X]
 * for a target polynomial of 64 coefficients.
 * @param[out] f0 Array of 2 elements
 * @param[out] f1 Array of 2 elements
 * @param[in] f Array of 4 elements
 */
void radix_1024_4(uint16_t* f0, uint16_t* f1, const uint16_t* f) {
  f0[0] = f[0];
  f0[1] = f[2]^f[3];
  f1[0] = f[1]^f0[1];
  f1[1] = f[3];
}

/**
 * Third recursive call of the radix conversion on GF(2^10)[X]
 * for a target polynomial of 64 coefficients.
 * @param[out] f0 Array of 4 elements
 * @param[out] f1 Array of 4 elements
 * @param[in] f Array of 8 elements
 */
void radix_1024_8(uint16_t* f0, uint16_t* f1, const uint16_t* f) {
  f0[0] = f[0];
  f0[2] = f[4]^f[6];
  f0[3] = f[6]^f[7];
  f1[1] = f[3]^f[5]^f[7];
  f1[2] = f[5]^f[6];
  f1[3] = f[7];
  f0[1] = f[2]^f0[2]^f1[1];
  f1[0] = f[1]^f0[1];
}

/**
 * Second recursive call of the radix conversion on GF(2^10)[X]
 * for a target polynomial of 64 coefficients.
 * @param[out] f0 Array of 8 elements
 * @param[out] f1 Array of 8 elements
 * @param[in] f Array of 16 elements
 */
void radix_1024_16(uint16_t* f0, uint16_t* f1, const uint16_t* f) {
  f0[4] = f[8]^f[12];
  f0[6] = f[12]^f[14];
  f0[7] = f[14]^f[15];
  f1[5] = f[11]^f[13];
  f1[6] = f[13]^f[14];
  f1[7] = f[15];
  f0[5] = f[10]^f[12]^f1[5];
  f1[4] = f[9]^f[13]^f0[5];

  f0[0] = f[0];
  f1[3] = f[7]^f[11]^f[15];
  f0[3] = f[6]^f[10]^f[14]^f1[3];
  f0[2] = f[4]^f0[4]^f0[3]^f1[3];
  f1[1] = f[3]^f[5]^f[9]^f[13]^f1[3];
  f1[2] = f[3]^f1[1]^f0[3];
  f0[1] = f[2]^f0[2]^f1[1];
  f1[0] = f[1]^f0[1];
}

/**
 * First recursive call of the radix conversion on GF(2^10)[X]
 * for a target polynomial of 64 coefficients.
 * @param[out] f0 Array of 16 elements
 * @param[out] f1 Array of 16 elements
 * @param[in] f Array of 32 elements
 */
void radix_1024_32(uint16_t* f0, uint16_t* f1, const uint16_t* f) {
  uint16_t Q[2*8];
  uint16_t R[2*8];

  memcpy(Q + 8, f + 3*8, 2*8);
  memcpy(R, f, 4*8);

  for(size_t i = 0 ; i < 8 ; ++i) {
    Q[i] = f[2*8 + i] ^ f[3*8 + i];
    R[8 + i] ^= Q[i];
  }

  uint16_t Q0[8];
  uint16_t Q1[8];
  uint16_t R0[8];
  uint16_t R1[8];

  radix_1024_16 (Q0, Q1, Q);
  radix_1024_16 (R0, R1, R);

  memcpy(f0, R0, 2*8);
  memcpy(f0 + 8, Q0, 2*8);
  memcpy(f1, R1, 2*8);
  memcpy(f1 + 8, Q1, 2*8);
}

/**
 * Computes the radix conversion of f,
 * that is f0 and f1 such that f(x) = f0(x^2-x) + x.f1(x^2-x)
 * @param[out] f0 Array of 32 elements
 * @param[out] f1 Array of 32 elements
 * @param[in] f Array of 64 elements
 */
void radix_1024_64(uint16_t* f0, uint16_t* f1, const uint16_t* f) {
  uint16_t Q[2*16];
  uint16_t R[2*16];

  memcpy(Q + 16, f + 3*16, 2*16);
  memcpy(R, f, 4*16);

  for(size_t i = 0 ; i < 16 ; ++i) {
    Q[i] = f[2*16 + i] ^ f[3*16 + i];
    R[16 + i] ^= Q[i];
  }

  uint16_t Q0[16];
  uint16_t Q1[16];
  uint16_t R0[16];
  uint16_t R1[16];

  radix_1024_32 (Q0, Q1, Q);
  radix_1024_32 (R0, R1, R);
     
  memcpy(f0, R0, 2*16);
  memcpy(f0 + 16, Q0, 2*16);
  memcpy(f1, R1, 2*16);
  memcpy(f1 + 16, Q1, 2*16);
}

/**
 * Fifth and last recursive call of the additive FFT on GF(2^10) for a 64-coefficient polynomial.
 * @param[out] w Array of 32 elements
 * @param[in] f Array of 2 elements
 */
void fft_1024_2 (uint16_t* w, const uint16_t* f) {
  uint16_t tmp[5];
  for(size_t i = 0 ; i < 5 ; ++i)
    tmp[i] = gf_mul(fft_final_betas[i], f[1]);
  w[0] = f[0];
  for(size_t j = 0 ; j < 5 ; ++j) {
    for(size_t k = 0 ; k < 1U << j ; ++k) {
      w[(1 << j) + k] = w[k] ^ tmp[j];
    }
  }
}

/**
 * Fourth recursive call of the additive FFT on GF(2^10) for a 64-coefficient polynomial.
 * @param[out] w Array of 64 elements
 * @param[in] f Array of 4 elements
 * @param[in] f_coeffs Number of coefficients of f (<= 4)
 */
void fft_1024_4(uint16_t* w, uint16_t* f, size_t f_coeffs) {
  for(size_t i = 1 ; i < (1 << 2) ; ++i)
    f[i] = gf_mul(betas_pows[64+32+16+i], f[i]);
  
  uint16_t f0[1 << (2-1)];
  uint16_t f1[1 << (2-1)];

  radix_1024_4 (f0, f1, f);

  uint16_t u[1 << (6-1)];
  uint16_t v[1 << (6-1)];
        
  fft_1024_2 (u, f0);
        
  if(f_coeffs <= 3) {
    size_t k = 1 << (6-1);
    w[0] = u[0];
    w[k] = u[0] ^ f1[0];
    for(size_t i = 1 ; i < k ; ++i) {
      w[i] = u[i] ^ gf_mul(gammas_sums[512+256+128+64+i], f1[0]);
      w[k+i] = w[i] ^ f1[0];
    }
  }
  else {
    fft_1024_2 (v, f1);

    size_t k = 1 << (6-1);
    memcpy(w + k, v, 2*k);
    w[0] = u[0];
    w[k] ^= u[0];
    for(size_t i = 1 ; i < k ; ++i) {
      w[i] = u[i] ^ gf_mul(gammas_sums[512+256+128+64+i], v[i]);
      w[k+i] ^= w[i];
    }
  }
}

/**
 * Third recursive call of the additive FFT on GF(2^10) for a 64-coefficient polynomial.
 * @param[out] w Array of 128 elements
 * @param[in] f Array of 8 elements
 * @param[in] f_coeffs Number of coefficients of f (<= 8)
 */
void fft_1024_8(uint16_t* w, uint16_t* f, size_t f_coeffs) {
  for(size_t i = 1 ; i < (1 << 3) ; ++i)
    f[i] = gf_mul(betas_pows[64+32+i], f[i]);

  uint16_t f0[1 << (3-1)];
  uint16_t f1[1 << (3-1)];

  radix_1024_8 (f0, f1, f);

  uint16_t u[1 << (7-1)];
  uint16_t v[1 << (7-1)];

  fft_1024_4 (u, f0, (f_coeffs+1)/2);
  fft_1024_4 (v, f1, f_coeffs/2);

  size_t k = 1 << (7-1);
  memcpy(w + k, v, 2*k);
  w[0] = u[0];
  w[k] ^= u[0];
  for(size_t i = 1 ; i < k ; ++i) {
    w[i] = u[i] ^ gf_mul(gammas_sums[512+256+128+i], v[i]);
    w[k+i] ^= w[i];
  }
}

/**
 * Second recursive call of the additive FFT on GF(2^10) for a 64-coefficient polynomial.
 * @param[out] w Array of 256 elements
 * @param[in] f Array of 16 elements
 * @param[in] f_coeffs Number of coefficients of f (<= 16)
 */
void fft_1024_16(uint16_t* w, uint16_t* f, size_t f_coeffs) {
  for(size_t i = 1 ; i < (1 << 4) ; ++i)
    f[i] = gf_mul(betas_pows[64+i], f[i]);

  uint16_t f0[1 << (4-1)];
  uint16_t f1[1 << (4-1)];

  radix_1024_16 (f0, f1, f);

  uint16_t u[1 << (8-1)];
  uint16_t v[1 << (8-1)];

  fft_1024_8 (u, f0, (f_coeffs+1)/2);
  fft_1024_8 (v, f1, f_coeffs/2);

  size_t k = 1 << (8-1);
  memcpy(w + k, v, 2*k);
  w[0] = u[0];
  w[k] ^= u[0];
  for(size_t i = 1 ; i < k ; ++i) {
    w[i] = u[i] ^ gf_mul(gammas_sums[512+256+i], v[i]);
    w[k+i] ^= w[i];
  }
}

/**
 * First recursive call of the additive FFT on GF(2^10) for a 64-coefficient polynomial.
 * @param[out] w Array of 512 elements
 * @param[in] f Array of 32 elements
 * @param[in] f_coeffs Number of coefficients of f (<= 32)
 */
void fft_1024_32(uint16_t* w, uint16_t* f, size_t f_coeffs) {
  uint16_t f0[1 << (5-1)];
  uint16_t f1[1 << (5-1)];

  radix_1024_32 (f0, f1, f);

  uint16_t u[1 << (9-1)];
  uint16_t v[1 << (9-1)];

  fft_1024_16 (u, f0, (f_coeffs+1)/2);
  fft_1024_16 (v, f1, f_coeffs/2);

  size_t k = 1 << (9-1);
  memcpy(w + k, v, 2*k);
  w[0] = u[0];
  w[k] ^= u[0];
  for(size_t i = 1 ; i < k ; ++i) {
    w[i] = u[i] ^ gf_mul(gammas_sums[512+i], v[i]);
    w[k+i] ^= w[i];
  }
}

/**
 * Additive FFT on GF(2^10) for a 64-coefficient polynomial.
 * Evaluates polynomial f of degree f_coeffs-1 (less than 64) on all elements of GF(2^10).
 * @param[out] w Array of 1024 elements
 * @param[in] f Array of 64 elements
 * @param[in] f_coeffs Number of coefficients of f (less than 64)
 */
void fft_1024_64(uint16_t* w, uint16_t* f, size_t f_coeffs) {
  uint16_t f0[32];
  uint16_t f1[32];

  radix_1024_64 (f0, f1, f);

  uint16_t u[512];
  uint16_t v[512];

  fft_1024_32 (u, f0, (f_coeffs+1)/2);
  fft_1024_32 (v, f1, f_coeffs/2);

  size_t k = 1 << (10-1);
  memcpy(w + k, v, 2*k);
  w[0] = u[0];
  w[k] ^= u[0];

  for(size_t i = 1 ; i < k ; ++i) {
    w[i] = u[i] ^ gf_mul(gammas_sums[i], v[i]);
    w[k+i] ^= w[i];
  }
}

/**
 * Arranges the received word rcv in a form w such that applying the additive FFT transpose to w
 * yields the BCH syndromes of the received word rcv.
 * Since the received word rcv gives coefficients of the primitive element alpha, we twist accordingly. <br>
 * Furthermore, the additive FFT transpose needs elements indexed by their decomposition on the chosen basis,
 * so we apply the adequate permutation.
 * @param[out] w Array of size 2^PARAM_M
 * @param[in] rcv Array of size VEC_N1_SIZE_BYTES
 */
void fft_t_preprocess_bch_codeword(uint16_t* w, const uint8_t* rcv) {
  // Unpack the received word rcv into array r
  uint16_t r[1 << PARAM_M];
  size_t i;
  for(i = 0 ; i < VEC_N1_SIZE_BYTES-(PARAM_N1%8!=0) ; ++i)
    for(size_t j = 0 ; j < 8 ; ++j)
      r[8*i+j] = (rcv[i] >> j) & 1;

  for(size_t j = 0 ; j < PARAM_N1%8 ; ++j) // last byte
    r[8*i+j] = (rcv[i] >> j) & 1;

  memset(r+PARAM_N1, 0, 2*((1<<PARAM_M)-PARAM_N1)); // complete r with zeros

  // Twist and permute r adequately to obtain w
  size_t k = 1 << (PARAM_M-1);
  w[0] = 0;
  w[k] = -r[0] & 1;
  for(size_t i = 1 ; i < k ; ++i) {
    w[i] = -r[gf_log(gammas_sums[i])] & gammas_sums[i];
    w[k+i] = -r[gf_log(gammas_sums[i]^1)] & (gammas_sums[i]^1);
  }
}

/**
 * Retrieves the error polynomial err from the evaluations w of the ELP
 * (Error Locator Polynomial) on all field elements.
 * @param[out] err Array of size VEC_N1_SIZE_BYTES
 * @param[in] w Array of size 2^PARAM_M
 */
void fft_retrieve_bch_error_poly(uint8_t* err, const uint16_t* w) {
  size_t k = 1 << (PARAM_M-1);
  err[0] ^= 1 ^ ((uint16_t)-w[0] >> 15);
  size_t index = PARAM_GF_MUL_ORDER;
  uint8_t bit = 1 ^ ((uint16_t)-w[k] >> 15);
  err[index/8] ^= bit << (index%8);
  
  for(size_t i = 1 ; i < k ; ++i) {
    index = PARAM_GF_MUL_ORDER - gf_log(gammas_sums[i]);
    bit = 1 ^ ((uint16_t)-w[i] >> 15);
    err[index/8] ^= bit << (index%8);

    index = PARAM_GF_MUL_ORDER - gf_log(gammas_sums[i] ^ 1);
    bit = 1 ^ ((uint16_t)-w[k+i] >> 15);
    err[index/8] ^= bit << (index%8);
  }
}
