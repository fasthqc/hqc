#include <immintrin.h>



#define PARAM_N1_256 ((PARAM_N1/256)+1)

// number of 256-bit words to store 2*PARAM_DELTA syndromes.
// each syndrome is 16-bit long. 16 syndromes can be encoded
// in one 256-bit word.
#define PARAM_2DELTA_SIZE_256         (((PARAM_DELTA>>4)+1)<<1)

static const unsigned char mask2[32] = {
        0x01, 0x02, 0x04, 0x08,0x10, 0x20, 0x40, 0x80,
        0x01, 0x02, 0x04, 0x08,0x10, 0x20, 0x40, 0x80,
        0x01, 0x02, 0x04, 0x08,0x10, 0x20, 0x40, 0x80,
        0x01, 0x02, 0x04, 0x08,0x10, 0x20, 0x40, 0x80
    };



static const unsigned char mask1[32] = {
0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 
0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 
0x02, 0x02, 0x02, 0x02, 0x02, 0x02, 0x02, 0x02, 
0x03, 0x03, 0x03, 0x03, 0x03, 0x03, 0x03, 0x03
};


__m256i zero256 = (__m256i){0UL,0UL,0UL,0UL};
__m256i mask_un, mask_deux ;


//tables alpha_ij globale
int16_t *table_alpha_ij;
   gf_tables *tables;

syndrome_set* syndrome_init2() {
  syndrome_set* synd_set = (syndrome_set*) malloc(sizeof(syndrome_set));
  synd_set->size = ((PARAM_DELTA>>4)+1)<<5;
  synd_set->tab = (int16_t*) calloc((synd_set->size), sizeof(int16_t));
  
  return synd_set;
}



// function to initialize a table which contains values alpha^ij for i in [0,N1[ and j in [1,2*PARAM_DELTA]
// these values are used in order to compute the syndromes of the received word v(x)=v_0+v_1x+...+v_{n1-1}x^{n1-1}
// value alpha^ij is stored in alpha_ij_table[2*PARAM_DELTA*i+j-1]
// The syndromes are equal to v(alpha^k) for k in [1,2*PARAM_DELTA]
// Size of the table is fixed to match 256 bit representation
// Useless values are filled with 0.
int16_t* alpha_ij_table_init()
{
  int16_t* alpha_ij_table = (int16_t*)calloc(PARAM_N1*(PARAM_2DELTA_SIZE_256<<4),sizeof(int16_t));
  return alpha_ij_table;
}

void gf_generation2(gf_tables* gf_tables, int16_t* alpha_ij_table) {
  const uint16_t k  = 1 << PARAM_M; // k = 2^m = 2^10
  const uint16_t poly = PARAM_POLY; // get the primitive polynomial
  uint16_t val = 1;
  uint16_t alpha = 2; // alpha the root of the primitive polynomial is the primitive element of GF(2^10)
  int tmp_value;
  int16_t* alpha_temp;

  for(int i = 0 ; i < PARAM_GF_MUL_ORDER ; ++i){
    gf_tables->antilog_tab[i] = val;
    gf_tables->log_tab[val] = i;
    val = val * alpha; // by multiplying by alpha and reducing later if needed we generate all the elements of GF(2^10)
    if(val >= k){ // if val is greater than 2^10
      val ^= poly; // replace alpha^10 by alpha^3 + 1
    }
  }

  gf_tables->antilog_tab[PARAM_GF_MUL_ORDER] = 1; 
  gf_tables->log_tab[0] = -1; // by convention 
  
  // pre-computation of alpha^ij for i in [0, N1[ and j in [1, 2*PARAM_DELTA]
  // see comment of alpha_ij_table_init() function.
  for(uint16_t i = 0; i < PARAM_N1 ; ++i) {
    	tmp_value = 0;
	alpha_temp = alpha_ij_table+i*(PARAM_DELTA<<1);
      	for(uint16_t j = 0 ; j < (PARAM_DELTA<<1) ; j++) {      
		tmp_value = gf_mod(tmp_value + i);
		alpha_temp[j] = gf_get_antilog(gf_tables, tmp_value);
      	}
    
  }
} 


// Function to precompute table_alpha_ij

void syndrome_gen_init(){
   //gf_tables *tables;
	tables = gf_tables_init();
	table_alpha_ij = alpha_ij_table_init();

	gf_generation2(tables, table_alpha_ij);
  //gf_tables_clear(tables);  




	mask_un = _mm256_loadu_si256((__m256i*)mask1);
	mask_deux = _mm256_loadu_si256((__m256i*)mask2);


	
	printf("sortie syndrom_gen_init()\n");

}







// Our new AVX2 optimization
void syndrome_gen2(syndrome_set* synd_set, int16_t* table_alpha_ij, uint8_t* v) {

	// static variable so that it is stored in the DATA segment
	// not in the STACK segment
	static uint8_t tmp_array[PARAM_N1];

	__m256i y;
	__m256i S;
	// 256 bit register used to store 16 syndromes
	__m256i * synd_set256 = (__m256i *) synd_set->tab;
	__m256i *z = (__m256i*) tmp_array;

	uint32_t * aux;
	int16_t* alpha_temp;
	uint32_t i;

	//vectorized version of the separation of the coordinates of the vector v in order to put each coordinate in an unsigned char
	// aux is ued to consider 4 elements in v at each step of the loop
	aux = (uint32_t *) v;
	for (i = 0; i < (VEC_N1_SIZE_BYTES >> 2) << 2 ; i+=4) {
		 // duplicate aux 8 times in y , i.e y= (aux aux aux .... aux)
		 y = _mm256_set1_epi32(*aux);
		 // shuffle the bytes of y so that if aux=(a0 a1 a2 a3) 
		 // then y = (a0 a0 a0 a0 a0 a0 a0 a0 a1 a1 a1 a1 a1 a1 a1 a1 .... a3) 
		 y = _mm256_shuffle_epi8(y,mask_un);
		 // apply a mask on each byte of y to determine if jth bit of a_k is 0 or 1
		 z[i>>2] = _mm256_and_si256(y,mask_deux);
		 aux++;
	  }
	// last part of v
	uint8_t rest_bytes = (VEC_N1_SIZE_BYTES & 0x03) ; //(VEC_N1_SIZE_BYTES % 4)
	for (uint32_t k = 0; k < rest_bytes; k++)
	 {
	  for (uint32_t j = 0; j < 8 ; ++j) {
		tmp_array [(k+i) * 8 + j] = (v[i+k] >> j) & 0x01;
	 }
	}
	
	// Evaluation of the polynomial corresponding to the vector v in alpha^i for i in {1, ..., 2 * PARAM_DELTA}
	  for(uint16_t j = 0 ; j < PARAM_2DELTA_SIZE_256 ; ++j) {
	  	S = zero256;
		alpha_temp = table_alpha_ij+(j<<4);
		      
		for(uint16_t i = 0; i < PARAM_N1 ; ++i) {
			if(tmp_array[i]) {
				S^=_mm256_lddqu_si256((__m256i*)(alpha_temp+i*(PARAM_DELTA<<1)));
			}
		}
		_mm256_storeu_si256(synd_set256+j, S);
	   }
}







