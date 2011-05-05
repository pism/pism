# this makefile is an abbreviated form of src/snes/examples/tutorials/makefile

ALL: lapbasic lapdof2 poisson localVecMax ex2_printsome

# Get PETSc environment and rules:
PETSC3_VARIABLES := $(wildcard ${PETSC_DIR}/conf/variables)
# PETSC3_VARIABLES is empty string if not 3.1
ifeq ($(strip $(PETSC3_VARIABLES)),) 
  # try file which is only present in 3.0.0
  PETSC3_BASE := $(wildcard ${PETSC_DIR}/conf/base)
  ifeq ($(strip $(PETSC3_BASE)),) 
    # case of 2.3.3, etc.
    HAVE_PETSC3 ?= 0
    include ${PETSC_DIR}/bmake/common/base
  else
    # include file which is only present in 3.0.0
    HAVE_PETSC3 ?= 1
    include ${PETSC_DIR}/conf/base
  endif
else
  # include files which are only present in 3.1 and later
  HAVE_PETSC3 ?= 1
  include ${PETSC_DIR}/conf/variables
  include ${PETSC_DIR}/conf/rules
endif

# Compiler flags:
CFLAGS += -DHAVE_PETSC3=${HAVE_PETSC3} -g


# these examples are related to examples in the src/ksp/ksp/examples/tutorials/ and
# src/snes/examples/tutorials/ directories in a PETSC source release

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

ex2_printsome: ex2_printsome.o   chkopts
	-${CLINKER} -o ex2_printsome ex2_printsome.o  ${PETSC_SYS_LIB}
	${RM} -f ex2_printsome.o

