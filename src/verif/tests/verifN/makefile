ALL: bodvardsson

include ${PETSC_DIR}/conf/rules
include ${PETSC_DIR}/conf/variables

bodvardsson : bodvardsson.o ../exactTestN.o chkopts
	-${CLINKER} -o $@ $< ${PETSC_SNES_LIB} ../exactTestN.o
	rm -f bodvardsson.o

