/**
 * \file hqc.c
 * \brief Implementation of hqc.h
 */
 
#include "hqc.h"

// in order to use or not the patch fast convolution multiplication
#include "patch.h"


#ifdef PATCH
#include "patch.c"
#endif


/**
 *\fn void hqc_pke_keygen(unsigned char* pk, unsigned char* sk)
 *\brief Keygen of the HQC_PKE IND_CPA scheme
 *
 * The public key is composed of the syndrome <b>s</b> as well as the <b>seed</b> used to generate the vector <b>h</b>.
 *
 * The secret key is composed of the <b>seed</b> used to generate vectors <b>x</b> and  <b>y</b>.
 * As a technicality, the public key is appended to the secret key in order to respect NIST API.
 *
 * \param[out] pk String containing the public key
 * \param[out] sk String containing the secret key
 */
void hqc_pke_keygen(unsigned char* pk, unsigned char* sk) {

  // Create seed_expanders for public key and secret key
  unsigned char sk_seed[SEED_BYTES];
  randombytes(sk_seed, SEED_BYTES);
  AES_XOF_struct* sk_seedexpander = (AES_XOF_struct*) malloc(sizeof(AES_XOF_struct));
  seedexpander_init(sk_seedexpander, sk_seed, sk_seed + 32, SEEDEXPANDER_MAX_LENGTH);
    
  unsigned char pk_seed[SEED_BYTES];
  randombytes(pk_seed, SEED_BYTES);
  AES_XOF_struct* pk_seedexpander = (AES_XOF_struct*) malloc(sizeof(AES_XOF_struct));
  seedexpander_init(pk_seedexpander, pk_seed, pk_seed + 32, SEEDEXPANDER_MAX_LENGTH);

  // Compute secret key
  
  
#ifdef PATCH

	uint32_t vX[PARAM_OMEGA], vY[PARAM_OMEGA];
	vect_fixed_weight_by_position(vX, PARAM_OMEGA, sk_seedexpander);
	vect_fixed_weight_by_position(vY, PARAM_OMEGA, sk_seedexpander);
	
	
#else
  uint8_t x [VEC_N_SIZE_BYTES] = {0};
  uint8_t y [VEC_N_SIZE_BYTES] = {0};

  vect_fixed_weight(x, PARAM_OMEGA, sk_seedexpander);
  vect_fixed_weight(y, PARAM_OMEGA, sk_seedexpander);
#endif

  // Compute public key
#ifdef PATCH

	uint8_t h [SIZE8] = {0};
	vect_set_random(h, pk_seedexpander);
	uint8_t s [SIZE8] = {0};
	fastConvolutionMult((uint64_t *) h, vY,(uint64_t *) s, PARAM_OMEGA);
	vector_add_by_position_and_coordinate(s,vX,s, PARAM_OMEGA);
	
	
#else
  uint8_t h [VEC_N_SIZE_BYTES] = {0};
  vect_set_random(h, pk_seedexpander);
  uint8_t s [VEC_N_SIZE_BYTES] = {0};
  sparse_dense_mul(s, y, h, PARAM_OMEGA);
  vect_add(s, s, x, VEC_N_SIZE_BYTES);
  
  
#endif
  
  // Parse keys to string 
  hqc_public_key_to_string(pk, pk_seed, s);
  hqc_secret_key_to_string(sk, sk_seed, pk);

  #ifdef VERBOSE
    printf("\n\nsk_seed: "); for(int i = 0 ; i < SEED_BYTES ; ++i) printf("%02x", sk_seed[i]);
#ifdef PATCH
#else
    printf("\n\nx: "); vect_print(x, VEC_N_SIZE_BYTES);
    printf("\n\ny: "); vect_print(y, VEC_N_SIZE_BYTES);
#endif

    printf("\n\npk_seed: "); for(int i = 0 ; i < SEED_BYTES ; ++i) printf("%02x", pk_seed[i]);
    printf("\n\nh: "); vect_print(h, VEC_N_SIZE_BYTES);
    printf("\n\ns: "); vect_print(s, VEC_N_SIZE_BYTES);

    printf("\n\nsk: "); for(int i = 0 ; i < SECRET_KEY_BYTES ; ++i) printf("%02x", sk[i]);
    printf("\n\npk: "); for(int i = 0 ; i < PUBLIC_KEY_BYTES ; ++i) printf("%02x", pk[i]);
    printf("\n\nsk: "); for(int i = 0 ; i < SECRET_KEY_BYTES ; ++i) printf("%02x", sk[i]);
    printf("\n\npk: "); for(int i = 0 ; i < PUBLIC_KEY_BYTES ; ++i) printf("%02x", pk[i]);
  #endif

  free(sk_seedexpander);
  free(pk_seedexpander);
}

/**
 *\fn void hqc_pke_encrypt(uint8_t* u, uint8_t* v, uint8_t* m, unsigned char* theta, const unsigned char* pk)
 *\brief Encryption of the HQC_PKE IND_CPA scheme
 *
 * The cihertext is composed of vectors <b>u</b> and <b>v</b>.
 *
 * \param[out] u Vector u (first part of the ciphertext)
 * \param[out] v Vector v (second part of the ciphertext)
 * \param[in] m Vector representing the message to encrypt
 * \param[in] theta Seed used to derive randomness required for encryption
 * \param[in] pk String containing the public key
 */
void hqc_pke_encrypt(uint8_t* u, uint8_t* v, uint8_t* m, unsigned char* theta, const unsigned char* pk) {

  // Create seed_expander from theta
  AES_XOF_struct* seedexpander = (AES_XOF_struct*) malloc(sizeof(AES_XOF_struct));
  seedexpander_init(seedexpander, theta, theta + 32, SEEDEXPANDER_MAX_LENGTH);

  // Retrieve h and s from public key
  
  
  #ifdef PATCH
	uint8_t h [SIZE8] = {0};
	uint8_t s [SIZE8] = {0};
  #else
  uint8_t h [VEC_N_SIZE_BYTES] = {0};
  uint8_t s [VEC_N_SIZE_BYTES] = {0};
  #endif
  hqc_public_key_from_string(h, s, pk);
 
  // Generate r1, r2 and e
  #ifdef PATCH
	uint32_t v_r1[PARAM_OMEGA_R], v_r2[PARAM_OMEGA_R], v_e[PARAM_OMEGA_E];
	vect_fixed_weight_by_position(v_r1, PARAM_OMEGA_R, seedexpander);
	vect_fixed_weight_by_position(v_r2, PARAM_OMEGA_R, seedexpander);
	vect_fixed_weight_by_position(v_e, PARAM_OMEGA_E, seedexpander);
	
  #else
  uint8_t r1 [VEC_N_SIZE_BYTES] = {0};
  uint8_t r2 [VEC_N_SIZE_BYTES] = {0};
  uint8_t e [VEC_N_SIZE_BYTES] = {0};
  
  vect_fixed_weight(r1, PARAM_OMEGA_R, seedexpander);
  vect_fixed_weight(r2, PARAM_OMEGA_R, seedexpander);
  vect_fixed_weight(e, PARAM_OMEGA_E, seedexpander); 
  #endif

  // Compute u = r1 + r2.h
  #ifdef PATCH
	uint8_t tmp1 [SIZE8] = {0};
	uint8_t tmp2 [SIZE8] = {0};
	fastConvolutionMult((uint64_t *) h, v_r2,(uint64_t *) tmp1, PARAM_OMEGA_R);
	vector_add_by_position_and_coordinate(u,v_r1,tmp1, PARAM_OMEGA_R);
  #else
  sparse_dense_mul(u, r2, h, PARAM_OMEGA_R);
  vect_add(u, u, r1, VEC_N_SIZE_BYTES);
  #endif

  // Compute v = m.G by encoding the message
  tensor_code_encode(v, m);

  // Compute v = m.G + s.r2 + e
  #ifdef PATCH
	fastConvolutionMult((uint64_t *) s, v_r2,(uint64_t *) tmp2, PARAM_OMEGA_R);
	vect_resize(tmp1, PARAM_N, v, PARAM_N1N2);
	vector_add_by_position_and_coordinate(tmp2,v_e,tmp2, PARAM_OMEGA_R);
	vect_add(tmp2,tmp2,tmp1, VEC_N_SIZE_BYTES);
    vect_resize(v, PARAM_N1N2, tmp2, PARAM_N);
  #else
  uint8_t tmp1 [VEC_N_SIZE_BYTES];
  vect_resize(tmp1, PARAM_N, v, PARAM_N1N2);
  uint8_t tmp2 [VEC_N_SIZE_BYTES];
  sparse_dense_mul(tmp2, r2, s, PARAM_OMEGA_R);

  vect_add(tmp2, tmp2, e, VEC_N_SIZE_BYTES);
  vect_add(tmp2, tmp1, tmp2, VEC_N_SIZE_BYTES);
  vect_resize(v, PARAM_N1N2, tmp2, PARAM_N);
  #endif

  #ifdef VERBOSE
    printf("\n\nh: "); vect_print(h, VEC_N_SIZE_BYTES);
    printf("\n\ns: "); vect_print(s, VEC_N_SIZE_BYTES);
    #ifndef PATCH
    printf("\n\nr1: "); vect_print(r1, VEC_N_SIZE_BYTES);
    printf("\n\nr2: "); vect_print(r2, VEC_N_SIZE_BYTES);
    printf("\n\ne: "); vect_print(e, VEC_N_SIZE_BYTES);
    #endif
    printf("\n\ntmp2: "); vect_print(tmp2, VEC_N_SIZE_BYTES);

    printf("\n\nu: "); vect_print(u, VEC_N_SIZE_BYTES);
    printf("\n\nv: "); vect_print(v, VEC_N1N2_SIZE_BYTES);
  #endif

  free(seedexpander);
}

/**
 *\fn void hqc_pke_decrypt(uint8_t* m, uint8_t* u, uint8_t* v, const unsigned char* sk)
 *\brief Decryption of the HQC_PKE IND_CPA scheme
 *
 * \param[out] m Vector representing the decrypted message
 * \param[in] u Vector u (first part of the ciphertext)
 * \param[in] v Vector v (second part of the ciphertext)
 * \param[in] sk String containing the secret key
 */
void hqc_pke_decrypt(uint8_t* m, uint8_t* u, uint8_t* v, const unsigned char* sk) {

  
  #ifdef PATCH
  // Retrieve y, pk from secret key
	uint32_t vX[PARAM_OMEGA];
	uint32_t vY[PARAM_OMEGA];
	
	unsigned char sk_seed[SEED_BYTES];
	memcpy(sk_seed, sk, SEED_BYTES);
	AES_XOF_struct* sk_seedexpander = (AES_XOF_struct*) malloc(sizeof(AES_XOF_struct));
	seedexpander_init(sk_seedexpander, sk_seed, sk_seed + 32, SEEDEXPANDER_MAX_LENGTH);

	vect_fixed_weight_by_position(vX, PARAM_OMEGA, sk_seedexpander);
	vect_fixed_weight_by_position(vY, PARAM_OMEGA, sk_seedexpander);

	free(sk_seedexpander);
	
	// Compute v - u.y
	uint8_t tmp1 [SIZE8] = {0};
		
	uint8_t tmp2 [SIZE8] = {0};
	
	memcpy(tmp1,u,VEC_N_SIZE_BYTES);
	
	fastConvolutionMult((uint64_t *) tmp1, vY,(uint64_t *) tmp2, PARAM_OMEGA);
	
	vect_resize(tmp1, PARAM_N, v, PARAM_N1N2);
	
  #else
  // Retrieve x, y, pk from secret key
  uint8_t x [VEC_N_SIZE_BYTES] = {0};
  uint8_t y [VEC_N_SIZE_BYTES] = {0};
  unsigned char pk[PUBLIC_KEY_BYTES];
  
  hqc_secret_key_from_string(x, y, pk, sk);
  // Compute v - u.y
  uint8_t tmp1 [VEC_N_SIZE_BYTES] = {0};
  vect_resize(tmp1, PARAM_N, v, PARAM_N1N2);
  uint8_t tmp2 [VEC_N_SIZE_BYTES] = {0};

  sparse_dense_mul(tmp2, y, u, PARAM_OMEGA);
  #endif
  
  vect_add(tmp2, tmp1, tmp2, VEC_N_SIZE_BYTES);

  #ifdef VERBOSE
    printf("\n\nu: "); vect_print(u, VEC_N_SIZE_BYTES);
    printf("\n\nv: "); vect_print(v, VEC_N1N2_SIZE_BYTES);
    #ifndef PATCH
    printf("\n\ny: "); vect_print(y, VEC_N_SIZE_BYTES);
	#endif
    printf("\n\nv - u.y: "); vect_print(tmp2, VEC_N_SIZE_BYTES);
  #endif
  
  // Compute m by decoding v - u.y
  tensor_code_decode(m, tmp2);
}
