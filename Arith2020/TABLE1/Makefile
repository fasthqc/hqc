OPTIONS= -O3 -mavx2 -mpclmul -std=c99 -lm -lgf2x


KaratMult : KaratMult.o fonctions.o
	gcc -o KaratMult KaratMult.o fonctions.o $(OPTIONS)

fonctions.o : fonctions.c fonctions.h
	gcc -c fonctions.c $(OPTIONS)

KaratMult.o : fonctions.o
	gcc -c KaratMult.c $(OPTIONS)

test: KaratMult
	./KaratMult 512
	./KaratMult 1024
	./KaratMult 2048
	./KaratMult 4096
	./KaratMult 8192
	./KaratMult 16384
	./KaratMult 32768
	./KaratMult 65536
	./KaratMult 131072

all : KaratMult


clean :
	rm *.o -fv *~

cleanall : clean
	rm KaratMult -f

