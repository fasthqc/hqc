OPTIONS= -O3 -mavx2 -mpclmul -std=c99 -g -lm 
WARNINGS= -Wall -pedantic -Wextra
MACROS= #-DWNAF#-DFLAG#-DWNAF#-DPROBE# -DREALTIME -DPROC1  -DPROC0 #-DPROBEFILE -DREALTIME
OPTIONS+= $(WARNINGS)
OPTIONS+= $(MACROS)
FastConvolution : FastConvolution.o fonctions.o
	gcc -o FastConvolution FastConvolution.o fonctions.o $(OPTIONS)

fonctions.o : fonctions.c fonctions.h
	gcc -c fonctions.c $(OPTIONS)

FastConvolution.o : fonctions.o
	gcc -c FastConvolution.c $(OPTIONS)

all : FastConvolution


clean :
	rm *.o -fv *~

cleanall : clean
	rm FastConvolution -f

