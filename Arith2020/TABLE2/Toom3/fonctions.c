/*****************************************************************



*****************************************************************/



#include "fonctions.h"




/***************************************************************

	Fonctions d'affichage

*/
  

void afficheVect(unsigned long int *A, char *var, int size)
{
	int i;
	unsigned long int tmp;
	printf("%s := ",var);
	
	for(i=0;i<size;i++){
		tmp=0;
		for(int j=0;j<WORD;j++) tmp^= ((A[i]>>j)&1UL)<<(WORD-1-j);
		printf("%16.16lX ",tmp);
	}
	printf("\n");
}



#include "ToomCookMult.c"

