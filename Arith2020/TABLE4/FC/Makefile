OPTIONS= -O3 -g -mavx2 -mpclmul -lm  #

ifdef SIZE
        MACROS= -DPARAM_N=$(SIZE) -DPARAM_OMEGA=$(OMEGA) -DP$(SIZE)=1
else      
        MACROS= -DPARAM_N=24677 -DPARAM_OMEGA=67 -DP$(SIZE)=1
endif


OPTIONS+=$(MACROS)

FastConvolution : FastConvolution.o fonctions.o
	gcc -o FastConvolution FastConvolution.o fonctions.o $(OPTIONS)

fonctions.o : fonctions.c fonctions.h gf2x_avx2.h
	gcc -c fonctions.c $(OPTIONS)

FastConvolution.o : fonctions.o
	gcc -c FastConvolution.c $(OPTIONS)

test:
	make -B SIZE=24677 OMEGA=67
	./FastConvolution
	make -B SIZE=43669 OMEGA=101
	./FastConvolution
	make -B SIZE=46747 OMEGA=101
	./FastConvolution
	make -B SIZE=63587 OMEGA=133
	./FastConvolution
	make -B SIZE=67699 OMEGA=133
	./FastConvolution
	make -B SIZE=70853 OMEGA=133
	./FastConvolution
	make -B SIZE=24677 OMEGA=77
	./FastConvolution
	make -B SIZE=43669 OMEGA=117
	./FastConvolution
	make -B SIZE=46747 OMEGA=117
	./FastConvolution
	make -B SIZE=63587 OMEGA=153
	./FastConvolution
	make -B SIZE=67699 OMEGA=153
	./FastConvolution
	make -B SIZE=70853 OMEGA=153
	./FastConvolution

	
all : FastConvolution


clean :
	rm *.o -fv *~

cleanall : clean
	rm FastConvolution -f

