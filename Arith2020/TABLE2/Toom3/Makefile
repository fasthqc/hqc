OPTIONS= -O3 -mavx2 -mpclmul -std=c99 -lm -lgf2x

ifdef SIZE
        MACROS= -DPARAM_N=$(SIZE)
else      
        MACROS= -DPARAM_N=48768  
endif


OPTIONS+=$(MACROS)

ToomCook : ToomCook.o fonctions.o
	gcc -o ToomCook ToomCook.o fonctions.o $(OPTIONS)

fonctions.o : fonctions.c fonctions.h
	gcc -c fonctions.c $(OPTIONS)

ToomCook.o : fonctions.o
	gcc -c ToomCook.c $(OPTIONS)

test:
	make -B SIZE=2688
	./ToomCook
	make -B SIZE=5760
	./ToomCook
	make -B SIZE=11904
	./ToomCook
	make -B SIZE=24192
	./ToomCook
	make -B SIZE=48768
	./ToomCook
	make -B SIZE=97920
	./ToomCook
	
all : ToomCook


clean :
	rm *.o -fv *~

cleanall : clean
	rm ToomCook -f

