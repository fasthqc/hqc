#define _GNU_SOURCE

#include <unistd.h>
#include <stdio.h>
#include <sys/syscall.h>

#include "api.h"
#include "parameters.h"

#define NTEST 1000
#define NSAMPLES 50


/**** Measurements procedures according to INTEL white paper

 "How to benchmark code execution times on INTEL IA-32 and IA-64" 
 
 *****/

inline static uint64_t cpucyclesStart (void) {

	unsigned hi, lo;
	__asm__ __volatile__ (	"CPUID\n\t"
				"RDTSC\n\t"
				"mov %%edx, %0\n\t"
				"mov %%eax, %1\n\t"
				: "=r" (hi), "=r" (lo)
				:
				: "%rax", "%rbx", "%rcx", "%rdx");

	return ((uint64_t)lo)^(((uint64_t)hi)<<32);


}

inline static uint64_t cpucyclesStop (void) {

	unsigned hi, lo;
	__asm__ __volatile__(	"RDTSCP\n\t"
				"mov %%edx, %0\n\t"
				"mov %%eax, %1\n\t"
				"CPUID\n\t"
				: "=r" (hi), "=r" (lo)
				:
				: "%rax", "%rbx", "%rcx", "%rdx");
				
	return ((uint64_t)lo)^(((uint64_t)hi)<<32);


}


int main() {

  unsigned char pk[PUBLIC_KEY_BYTES];
  unsigned char sk[SECRET_KEY_BYTES];

  unsigned char ct[CIPHERTEXT_BYTES];
  unsigned char ss1[SHARED_SECRET_BYTES];
  unsigned char ss2[SHARED_SECRET_BYTES];
  
  unsigned char seed[48];
  syscall(SYS_getrandom, seed, 48, 0);
  randombytes_init(seed, NULL, 256);
  
   	unsigned long long timer , meanTimer1 =0, meanTimer2 =0, meanTimer3 =0, t1,t2;


#ifndef VALGRIND
  
	// cache memory heating
	for(int i=0;i<NTEST;i++)
	{
		crypto_kem_keypair(pk, sk);
	}
  
	// timing
	for(int i=0;i<NSAMPLES;i++)
	{
		crypto_kem_keypair(pk, sk);
		timer = (unsigned long long int)0x1<<63;
		for(int j=0;j<NTEST;j++)
		{
			t1 = cpucyclesStart();
			crypto_kem_keypair(pk, sk);
			t2 = cpucyclesStop();

			if(timer>t2-t1) timer = t2-t1;
		}
		
		meanTimer1 += timer;
	}


	// cache memory heating
	for(int i=0;i<NTEST;i++)
	{
		crypto_kem_enc(ct, ss1, pk);
	}
	
	// timing
	for(int i=0;i<NSAMPLES;i++)
	{
		crypto_kem_keypair(pk, sk);
		timer = (unsigned long long int)0x1<<63;
		for(int j=0;j<NTEST;j++)
		{
			t1 = cpucyclesStart();
			crypto_kem_enc(ct, ss1, pk);
			t2 = cpucyclesStop();

			if(timer>t2-t1) timer = t2-t1;
		}
		
		meanTimer2 += timer;
	}//*/

	// cache memory heating
	for(int i=0;i<NTEST;i++)
	{
		crypto_kem_dec(ss2, ct, sk);
	}
	
	int flag=0;
	// timing
	for(int i=0;i<NSAMPLES;i++)
	{
		crypto_kem_keypair(pk, sk);
		crypto_kem_enc(ct, ss1, pk);
		if (crypto_kem_dec(ss2, ct, sk)) flag++;
		timer = (unsigned long long int)0x1<<63;

	  if(memcmp(ss1,ss2,SHARED_SECRET_BYTES)) {
		printf("ERROR\n");
		exit(0);
		
	  }
	  //else   printf("i = %d\n",i);
		for(int j=0;j<NTEST;j++)
		{
			t1 = cpucyclesStart();
			crypto_kem_dec(ss2, ct, sk);
			t2 = cpucyclesStop();
			if(timer>t2-t1) timer = t2-t1;

		}
		
		meanTimer3 += timer;
	}
	
  if(flag) printf("%d aborts !!!!\n", flag);
  printf("\nKeygen: %lld CPU cycles", meanTimer1/NSAMPLES);
  printf("\nEncaps: %lld CPU cycles", meanTimer2/NSAMPLES);
  printf("\nDecaps: %lld CPU cycles", meanTimer3/NSAMPLES);
  printf("\n");

#else

  crypto_kem_keypair(pk, sk);
  
  crypto_kem_enc(ct, ss1, pk);

  crypto_kem_dec(ss2, ct, sk);
  
  if(memcmp(ss1,ss2,SHARED_SECRET_BYTES)) {
	printf("ERROR\n");
	exit(0);
  }

#endif

  return 0;
}




