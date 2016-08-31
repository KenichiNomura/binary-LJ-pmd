CC=mpicc

pmd: pmd.c pmd.h
	$(CC) -lm -std=c99 -O3 pmd.c -o pmd	

clean: 
	rm -rfv gr pmd pmd.d *.xyz pmd.dSYM *.d
