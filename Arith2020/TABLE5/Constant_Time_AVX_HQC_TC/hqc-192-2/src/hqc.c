/**
 * @file hqc.c
 * @brief Implementation of hqc.h
 */
 
#include "hqc.h"

/**
 * @brief Keygen of the HQC_PKE IND_CPA scheme
 *
 * The public key is composed of the syndrome <b>s</b> as well as the <b>seed</b> used to generate the vector <b>h</b>.
 *
 * The secret key is composed of the <b>seed</b> used to generate vectors <b>x</b> and  <b>y</b>.
 * As a technicality, the public key is appended to the secret key in order to respect NIST API.
 *
 * @param[out] pk String containing the public key
 * @param[out] sk String containing the secret key
 */
void hqc_pke_keygen(unsigned char* pk, unsigned char* sk) {

  // Create seed_expanders for public key and secret key
  unsigned char sk_seed[SEED_BYTES];
  randombytes(sk_seed, SEED_BYTES);
  AES_XOF_struct sk_seedexpander;
  seedexpander_init(&sk_seedexpander, sk_seed, sk_seed + 32, SEEDEXPANDER_MAX_LENGTH);
    
  unsigned char pk_seed[SEED_BYTES];
  randombytes(pk_seed, SEED_BYTES);
  AES_XOF_struct pk_seedexpander;
  seedexpander_init(&pk_seedexpander, pk_seed, pk_seed + 32, SEEDEXPANDER_MAX_LENGTH);

  // Compute secret key
  uint8_t x [VEC_N_SIZE_BYTES] = {0};
  uint32_t y [PARAM_OMEGA] = {0};

  vect_fixed_weight(x, PARAM_OMEGA, &sk_seedexpander);
  vect_fixed_weight_by_coordinates(y, PARAM_OMEGA, &sk_seedexpander);

  // Compute public key
  uint8_t h [SIZE8] = {0};
  vect_set_random(h, &pk_seedexpander);
  uint8_t s [SIZE8] = {0};
  
 
  uint8_t y_c [SIZE8] = {0};
  initFWVector((uint64_t *) y_c, y, PARAM_OMEGA );

  TOOM3Mult((uint64_t*) h, (uint64_t*) y_c, (uint64_t*) s);
 
  vect_add(s, x, s, VEC_N_SIZE_BYTES);
  
  // Parse keys to string 
  hqc_public_key_to_string(pk, pk_seed, s);
  hqc_secret_key_to_string(sk, sk_seed, pk);

  #ifdef VERBOSE
    printf("\n\nsk_seed: "); for(int i = 0 ; i < SEED_BYTES ; ++i) printf("%02x", sk_seed[i]);
    printf("\n\nx: "); vect_print(x, VEC_N_SIZE_BYTES);
    printf("\n\ny: "); vect_print_sparse(y, PARAM_OMEGA);

    printf("\n\npk_seed: "); for(int i = 0 ; i < SEED_BYTES ; ++i) printf("%02x", pk_seed[i]);
    printf("\n\nh: "); vect_print(h, VEC_N_SIZE_BYTES);
    printf("\n\ns: "); vect_print(s, VEC_N_SIZE_BYTES);

    printf("\n\nsk: "); for(int i = 0 ; i < SECRET_KEY_BYTES ; ++i) printf("%02x", sk[i]);
    printf("\n\npk: "); for(int i = 0 ; i < PUBLIC_KEY_BYTES ; ++i) printf("%02x", pk[i]);
  #endif
}

/**
 * @brief Encryption of the HQC_PKE IND_CPA scheme
 *
 * The cihertext is composed of vectors <b>u</b> and <b>v</b>.
 *
 * @param[out] u Vector u (first part of the ciphertext)
 * @param[out] v Vector v (second part of the ciphertext)
 * @param[in] m Vector representing the message to encrypt
 * @param[in] theta Seed used to derive randomness required for encryption
 * @param[in] pk String containing the public key
 */
void hqc_pke_encrypt(uint8_t* u, uint8_t* v, uint8_t* m, unsigned char* theta, const unsigned char* pk) {

  // Create seed_expander from theta
  AES_XOF_struct seedexpander;
  seedexpander_init(&seedexpander, theta, theta + 32, SEEDEXPANDER_MAX_LENGTH);

  // Retrieve h and s from public key
  uint8_t h [SIZE8] = {0};
  uint8_t s [SIZE8] = {0};
  hqc_public_key_from_string(h, s, pk);
 
  // Generate r1, r2 and e
  uint8_t r1 [VEC_N_SIZE_BYTES] = {0};
  uint32_t r2 [PARAM_OMEGA_R] = {0};
  uint8_t e [VEC_N_SIZE_BYTES] = {0};
  
  vect_fixed_weight(r1, PARAM_OMEGA_R, &seedexpander);
  vect_fixed_weight_by_coordinates(r2, PARAM_OMEGA_R, &seedexpander);
  vect_fixed_weight(e, PARAM_OMEGA_E, &seedexpander); 

  // Compute u = r1 + r2.h
	uint8_t tmp2 [SIZE8] = {0};
	uint8_t r2_c [SIZE8] = {0};
	initFWVector((uint64_t *) r2_c, r2, PARAM_OMEGA );
	TOOM3Mult((uint64_t*) h, (uint64_t*) r2_c, (uint64_t*) tmp2);
	vect_add(u, r1, tmp2, VEC_N_SIZE_BYTES);
 
  // Compute v = m.G by encoding the message
  tensor_code_encode(v, m);
  uint8_t tmp1 [VEC_N_SIZE_BYTES];
  vect_resize(tmp1, PARAM_N, v, PARAM_N1N2);

  // Compute v = m.G + s.r2 + e
	TOOM3Mult((uint64_t*) s, (uint64_t*) r2_c, (uint64_t*) tmp2);

  vect_add(tmp2, e, tmp2, VEC_N_SIZE_BYTES);
  vect_add(tmp2, tmp1, tmp2, VEC_N_SIZE_BYTES);
  vect_resize(v, PARAM_N1N2, tmp2, PARAM_N);

  #ifdef VERBOSE
    printf("\n\nh: "); vect_print(h, VEC_N_SIZE_BYTES);
    printf("\n\ns: "); vect_print(s, VEC_N_SIZE_BYTES);
    printf("\n\nr1: "); vect_print(r1, VEC_N_SIZE_BYTES);
    printf("\n\nr2: "); vect_print_sparse(r2, PARAM_OMEGA_R);
    printf("\n\ne: "); vect_print(e, VEC_N_SIZE_BYTES);
    printf("\n\ntmp2: "); vect_print(tmp2, VEC_N_SIZE_BYTES);

    printf("\n\nu: "); vect_print(u, VEC_N_SIZE_BYTES);
    printf("\n\nv: "); vect_print(v, VEC_N1N2_SIZE_BYTES);
  #endif
}

/**
 * @brief Decryption of the HQC_PKE IND_CPA scheme
 *
 * @param[out] m Vector representing the decrypted message
 * @param[in] u Vector u (first part of the ciphertext)
 * @param[in] v Vector v (second part of the ciphertext)
 * @param[in] sk String containing the secret key
 */
void hqc_pke_decrypt(uint8_t* m, uint8_t* u, uint8_t* v, const unsigned char* sk) {

  // Retrieve x, y, pk from secret key
  uint8_t x [VEC_N_SIZE_BYTES] = {0};
  uint32_t y [PARAM_OMEGA] = {0};
  unsigned char pk[PUBLIC_KEY_BYTES];
  hqc_secret_key_from_string(x, y, pk, sk);
  
  // Compute v - u.y
	uint8_t tmp1 [SIZE8] = {0};
		
	uint8_t tmp2 [SIZE8] = {0};
	
	memcpy(tmp1,u,VEC_N_SIZE_BYTES);
	
	
	uint8_t y_c [SIZE8] = {0};
	initFWVector((uint64_t *) y_c, y, PARAM_OMEGA );
	
	TOOM3Mult((uint64_t*) u, (uint64_t*) y_c, (uint64_t*) tmp2);

	vect_resize(tmp1, PARAM_N, v, PARAM_N1N2);


  vect_add(tmp2, tmp1, tmp2, VEC_N_SIZE_BYTES);

  #ifdef VERBOSE
    printf("\n\nu: "); vect_print(u, VEC_N_SIZE_BYTES);
    printf("\n\nv: "); vect_print(v, VEC_N1N2_SIZE_BYTES);
    printf("\n\ny: "); vect_print_sparse(y, PARAM_OMEGA);
    printf("\n\nv - u.y: "); vect_print(tmp2, VEC_N_SIZE_BYTES);
  #endif
  
  // Compute m by decoding v - u.y
  tensor_code_decode(m, tmp2);
}
