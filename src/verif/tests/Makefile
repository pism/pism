# makefile merely for the local task of building the simple programs that 
# evaluate the exact solution codes themselves; no verification of a numerical
# scheme, much less PISM runs, here!

.PHONY : all
all : simpleE simpleH simpleL simpleM simpleN

# TYPICAL BY-HAND COMPILE-LINK-RUN EXAMPLE:
#   $ gcc -c simpleH.c -o simpleH.o
#   $ gcc -c exactTestH.c -o exactTestH.o
#   $ gcc -lm exactTestH.o simpleH.o -o simpleH
#   $ ./simpleH

simpleH : simpleH.o exactTestH.o
	$(CC) $^ -o $@ -lm

simpleL : simpleL.o exactTestL.o
	$(CC) $^ -o $@ -lm -lgsl -lgslcblas

simpleM : simpleM.o exactTestM.o
	$(CC) $^ -o $@ -lm -lgsl -lgslcblas

simpleN : simpleN.o exactTestN.o
	$(CC) $^ -o $@ -lm

simpleO : simpleO.o exactTestO.o
	$(CC) $^ -o $@

simpleP : simpleP.o exactTestP.o
	$(CC) $^ -o $@ -lm -lgsl -lgslcblas

.PHONY: clean
clean:
	@rm -f *.o *~ simpleH simpleL simpleM simpleN simpleO simpleP
