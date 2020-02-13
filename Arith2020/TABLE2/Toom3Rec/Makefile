OPTIONS= -O3 -mavx2 -mpclmul -std=c99 -lm -lgf2x

ifdef SIZE
        MACROS= -DPARAM_N=$(SIZE)
else      
        MACROS= -DPARAM_N=34176  
endif


OPTIONS+=$(MACROS)

ToomCook : ToomCook.o fonctions.o
	gcc -o ToomCook ToomCook.o fonctions.o $(OPTIONS)

fonctions.o : fonctions.c fonctions.h
	gcc -c fonctions.c $(OPTIONS)

ToomCook.o : fonctions.o
	gcc -c ToomCook.c $(OPTIONS)

test:
	make -B SIZE=6528
	./ToomCook
	make -B SIZE=15744
	./ToomCook
	make -B SIZE=34176
	./ToomCook
	make -B SIZE=71040
	./ToomCook
	make -B SIZE=144768
	./ToomCook
	
all : ToomCook


clean :
	rm *.o -fv *~

cleanall : clean
	rm ToomCook -f

