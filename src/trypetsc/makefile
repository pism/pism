# this makefile is an abbreviated form of src/snes/examples/tutorials/makefile

#MANSEC          = KSP SNES

PETSC3_BASE := $(wildcard ${PETSC_DIR}/conf/base)
# PETSC3_BASE is empty string if not 3.0.0
ifeq ($(strip $(PETSC3_BASE)),) 
  # include file which is only present in 2.3.3 and earlier
  HAVE_PETSC3 ?= 0
  include ${PETSC_DIR}/bmake/common/base
else
  # include file which is only present in 3.0.0 and later
  HAVE_PETSC3 ?= 1
  include ${PETSC_DIR}/conf/base
endif

# Compiler flags:
CFLAGS += -DHAVE_PETSC3=${HAVE_PETSC3} -g


lapbasic: lapbasic.o  chkopts
	-${CLINKER} -o lapbasic lapbasic.o ${PETSC_KSP_LIB}
	${RM} lapbasic.o

lapdof2: lapdof2.o  chkopts
	-${CLINKER} -o lapdof2 lapdof2.o ${PETSC_KSP_LIB}
	${RM} lapdof2.o

poisson: poisson.o  chkopts
	-${CLINKER} -o poisson poisson.o ${PETSC_SNES_LIB}
	${RM} poisson.o

localVecMax: localVecMax.o  chkopts
	-${CLINKER} -o localVecMax localVecMax.o ${PETSC_KSP_LIB}
	${RM} localVecMax.o

ex5_petsc2: ex5_petsc2.o chkopts
	-${CLINKER} -o ex5_petsc2 ex5_petsc2.o ${PETSC_SNES_LIB}
	${RM} ex5_petsc2.o

ex5_petsc3: ex5_petsc3.o chkopts
	-${CLINKER} -o ex5_petsc3 ex5_petsc3.o ${PETSC_SNES_LIB}
	${RM} ex5_petsc3.o

ex2_printsome: ex2_printsome.o   chkopts
	-${CLINKER} -o ex2_printsome ex2_printsome.o  ${PETSC_SYS_LIB}
	${RM} -f ex2_printsome.o

