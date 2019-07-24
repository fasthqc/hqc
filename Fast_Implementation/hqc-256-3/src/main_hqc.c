
#include <stdio.h>

#include "parameters.h"

#include "patch.h"

#define CPU_SPEED 3400000.0
//#define PATCH2
#ifdef PATCH2    
  #include "bch.h"
#endif



#include "../ccount.h"

#define NSAMPLES 50

#define NTEST 1000



int main() {

  printf("\n");
  printf("********************\n");
  printf("**** HQC-%d-%d ****\n", PARAM_SECURITY, PARAM_DFR_EXP);
  printf("********************\n");

  printf("\n");
  printf("N: %d   ", PARAM_N);
  printf("N1: %d   ", PARAM_N1);
  printf("N2: %d   ", PARAM_N2);
  printf("OMEGA: %d   ", PARAM_OMEGA);
  printf("OMEGA_R: %d   ", PARAM_OMEGA_R);
  printf("Failure rate: 2^-%d   ", PARAM_DFR_EXP);
  printf("Sec: %d bits", PARAM_SECURITY);
  printf("\n");

  unsigned char pk[PUBLIC_KEY_BYTES];
  unsigned char sk[SECRET_KEY_BYTES];
  unsigned char ct[CIPHERTEXT_BYTES];
  unsigned char key1[SHARED_SECRET_BYTES];
  unsigned char key2[SHARED_SECRET_BYTES];



#ifdef PATCH2    
  syndrome_gen_init();
#endif

  crypto_kem_keypair(pk, sk);
  
  
  
  crypto_kem_enc(ct, key1, pk);
  
  crypto_kem_dec(key2, ct, sk);

  printf("\n\nsecret1: ");
  for(int i = 0 ; i < SHARED_SECRET_BYTES ; ++i) printf("%x", key1[i]);

  printf("\nsecret2: ");
  for(int i = 0 ; i < SHARED_SECRET_BYTES ; ++i) printf("%x", key2[i]);
  printf("\n\n");
  //goto fin;
  
	/**************** CHRONOMETRAGE ********************/

	unsigned long long int timer, meanTimer = 0,START,STOP;
	
	printf("Chronométrage crypto_kem_keypair(pk, sk) !\n\n");
	meanTimer = 0;
	// chauffer les caches
	for(int i=0;i<NTEST;i++)
	{
		crypto_kem_keypair(pk, sk);
	}
	
	//Chronométrage !
	for(int i=0;i<NSAMPLES;i++)
	{
		crypto_kem_keypair(pk, sk);
		crypto_kem_enc(ct, key1, pk);
		timer = (unsigned long long int)0x1<<63;
		for(int j=0;j<NTEST;j++)
		{
			STAMP(START)
			crypto_kem_keypair(pk, sk);
			STAMP(STOP)

			if(timer>STOP-START) timer = STOP-START;
		}
		
		meanTimer += timer;
		printf("-");
		fflush(stdout);
	}
		
	printf("\n\nmeanTimer      : %llu et %fms                           \n\n",meanTimer/NSAMPLES, (meanTimer/NSAMPLES)/CPU_SPEED);
  
	
	printf("Chronométrage crypto_kem_enc(ct, key1, pk) !\n\n");
	meanTimer = 0;
	// chauffer les caches
	for(int i=0;i<NTEST;i++)
	{
		crypto_kem_enc(ct, key1, pk);
	}
	
	//Chronométrage !
	for(int i=0;i<NSAMPLES;i++)
	{
		crypto_kem_keypair(pk, sk);
		timer = (unsigned long long int)0x1<<63;
		for(int j=0;j<NTEST;j++)
		{
			STAMP(START)
			crypto_kem_enc(ct, key1, pk);
			STAMP(STOP)

			if(timer>STOP-START) timer = STOP-START;
		}
		
		meanTimer += timer;
		printf("-");
		fflush(stdout);
	}
		
	printf("\n\nmeanTimer      : %llu et %fms                           \n\n",meanTimer/NSAMPLES, (meanTimer/NSAMPLES)/CPU_SPEED);
  
	printf("Chronométrage crypto_kem_dec(key2, ct, sk) !\n\n");
	meanTimer = 0;
	// chauffer les caches
	for(int i=0;i<NTEST;i++)
	{
		crypto_kem_dec(key2, ct, sk);
	}
	
	int flag=0;
	//Chronométrage !
	for(int i=0;i<NSAMPLES;i++)
	{
		crypto_kem_keypair(pk, sk);
		crypto_kem_enc(ct, key1, pk);
		if (crypto_kem_dec(key2, ct, sk)) flag++;
		timer = (unsigned long long int)0x1<<63;
		for(int j=0;j<NTEST;j++)
		{
			STAMP(START)
			crypto_kem_dec(key2, ct, sk);
			STAMP(STOP)

			if(timer>STOP-START) timer = STOP-START;
		}
		
		meanTimer += timer;
		printf("-");
		fflush(stdout);
	}
	
	if(flag) printf("%d aborts !!!!\n", flag);
		
	printf("\n\nmeanTimer      : %llu et %fms                           \n\n",meanTimer/NSAMPLES, (meanTimer/NSAMPLES)/CPU_SPEED);
  fin:
	printf("Au revoir et merci !\n");


}
