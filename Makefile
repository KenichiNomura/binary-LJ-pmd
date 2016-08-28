CC=mpicc

pmd: pmd.c pmd.h
	$(CC) pmd.c -o pmd	

clean: 
	rm -rfv gr pmd pmd.d *.xyz pmd.dSYM	
