.PHONY : all
all : simpleABCD simpleE simpleFG simpleH simpleI simpleJ simpleK simpleL simpleM simpleN

# TYPICAL BY-HAND COMPILE-LINK-RUN EXAMPLE:
#   $ gcc -c simpleABCD.c -o simpleABCD.o
#   $ gcc -c exactTestsABCDE.c -o exactTestsABCDE.o
#   $ gcc -lm exactTestsABCDE.o simpleABCD.o -o simpleABCD
#   $ ./simpleABCD

simpleABCD : simpleABCD.o exactTestsABCDE.o
	$(CC) -lm $^ -o $@

simpleE : simpleE.o exactTestsABCDE.o
	$(CC) -lm $^ -o $@

simpleFG : simpleFG.o exactTestsFG.o
	$(CC) -lm $^ -o $@

simpleH : simpleH.o exactTestH.o
	$(CC) -lm $^ -o $@

simpleI : simpleI.o exactTestsIJ.o
	$(CC) -lm $^ -o $@

simpleJ : simpleJ.o exactTestsIJ.o
	$(CC) -lm $^ -o $@

simpleK : simpleK.o exactTestK.o
	$(CC) -lm $^ -o $@

simpleL : simpleL.o exactTestL.o
	$(CC) -lm -lgsl -lgslcblas $^ -o $@

simpleM : simpleM.o exactTestM.o
	$(CC) -lm -lgsl -lgslcblas $^ -o $@

simpleN : simpleN.o exactTestN.o
	$(CC) -lm $^ -o $@

.PHONY: clean
clean:
	@rm -f *.o simpleABCD simpleE simpleFG simpleH simpleI simpleJ simpleK \
	           simpleL simpleM simpleN

