CC=mpicc

pmd: pmd.c pmd.h
	$(CC) pmd.c -lm -std=c99 -O3 -o pmd	

clean: 
	rm -rfv gr pmd pmd.d *.xyz pmd.dSYM *.d pmd.png
