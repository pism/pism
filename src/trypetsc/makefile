# this makefile is similar to that in src/snes/examples/tutorials/makefile

all: printsome lapbasic lapdof2 poisson localVecMax mgpoisson

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

printsome: printsome.o  chkopts
	-${CLINKER} -o printsome printsome.o ${PETSC_SYS_LIB}
	${RM} printsome.o

localVecMax: localVecMax.o  chkopts
	-${CLINKER} -o localVecMax localVecMax.o ${PETSC_SYS_LIB}
	${RM} localVecMax.o

lapbasic: lapbasic.o  chkopts
	-${CLINKER} -o lapbasic lapbasic.o ${PETSC_KSP_LIB}
	${RM} lapbasic.o

lapdof2: lapdof2.o  chkopts
	-${CLINKER} -o lapdof2 lapdof2.o ${PETSC_KSP_LIB}
	${RM} lapdof2.o

mgpoisson: mgpoisson.o  chkopts
	-${CLINKER} -o mgpoisson mgpoisson.o ${PETSC_KSP_LIB}
	${RM} mgpoisson.o

poisson: poisson.o  chkopts
	-${CLINKER} -o poisson poisson.o ${PETSC_SNES_LIB}
	${RM} poisson.o

