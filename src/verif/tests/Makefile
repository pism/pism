# makefile merely for the local task of building the simple programs that 
# evaluate the exact solution codes themselves; no verification of a numerical
# scheme, much less PISM runs, here!

.PHONY : all
all : simpleABCD simpleE simpleFG simpleH simpleI simpleJ simpleK simpleL simpleM simpleN

# TYPICAL BY-HAND COMPILE-LINK-RUN EXAMPLE:
#   $ gcc -c simpleABCD.c -o simpleABCD.o
#   $ gcc -c exactTestsABCDE.c -o exactTestsABCDE.o
#   $ gcc -lm exactTestsABCDE.o simpleABCD.o -o simpleABCD
#   $ ./simpleABCD

simpleABCD : simpleABCD.o exactTestsABCDE.o
	$(CC) $^ -o $@ -lm

simpleE : simpleE.o exactTestsABCDE.o
	$(CC) $^ -o $@ -lm

simpleFG : simpleFG.o exactTestsFG.o
	$(CC) $^ -o $@ -lm

simpleH : simpleH.o exactTestH.o
	$(CC) $^ -o $@ -lm

simpleI : simpleI.o exactTestsIJ.o
	$(CC) $^ -o $@ -lm 

simpleJ : simpleJ.o exactTestsIJ.o
	$(CC) $^ -o $@ -lm

simpleK : simpleK.o exactTestK.o
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
	@rm -f *.o *~ simpleABCD simpleE simpleFG simpleH simpleI simpleJ simpleK \
	              simpleL simpleM simpleN simpleO simpleP

