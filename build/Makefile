ifndef PISM_PREFIX
  $(error Please run make from the top-level PISM directory)
endif

# Get the current PISM revision:
url = $(shell svn info | grep URL)
rev = $(shell svnversion ../)
ifeq ($(findstring trunk,$(url)),trunk)
  PISM_REVISION = trunk: revision $(rev)
endif
ifeq ($(findstring branches,$(url)),branches)
  PISM_REVISION = branch $(lastword $(subst /, ,$(url))): revision $(rev)
endif
ifndef PISM_REVISION
  $(error Can't determine PISM's revision number. \
  Please define PISM_REVISION in the top-level Makefile)
endif

DISPLAY = :0.0

# Get PETSc environment and rules:
PETSC3_BASE := $(wildcard ${PETSC_DIR}/conf/base)
# PETSC3_BASE is empty string if not 3.0.0
ifeq ($(strip $(PETSC3_BASE)),) 
  # include file which is only present in 2.3.3 and earlier
  PISM_HAVE_PETSC3 ?= 0
  include ${PETSC_DIR}/bmake/common/base
else
  # include file which is only present in 3.0.0 and later
  PISM_HAVE_PETSC3 ?= 1
  include ${PETSC_DIR}/conf/base
endif

# Compiler flags:
CFLAGS += -DPISM_HAVE_FFTW=${PISM_HAVE_FFTW} -DPISM_LOG_EVENTS=${PISM_LOG_EVENTS} \
          -DPISM_REVISION='"$(PISM_REVISION)"' \
          -DPISM_HAVE_PETSC3=${PISM_HAVE_PETSC3} \
          -Woverloaded-virtual -pipe

SHARED = -shared # default (for Linux)

# Linker flags:
TESTS_LIB_FLAGS := -L../lib -Wl,-rpath,$(PISM_PREFIX)/lib \
                   -ltests -lm -lgsl -lgslcblas
ICE_LIB_FLAGS := -lpism ${TESTS_LIB_FLAGS} ${PETSC_LIB} -lnetcdf 
ifeq (${PISM_HAVE_FFTW}, 1)
  ICE_LIB_FLAGS += -lfftw3
endif

# Override PETSc's choice of the linker (if asked to):
ifeq ($(USE_MPICXX), 1)
  CXXLINKER = $(subst mpicc,mpicxx,$(CLINKER))
endif

SHELL := /bin/sh
VPATH := $(subst src,../src,src/:src/base/:src/num/:src/verif/:src/eismint/:src/ismip/)

# Variables:
executables := pismr pismd pismv pisms
extras := simpleABCD simpleE simpleFG simpleH simpleI simpleJ \
	simpleL simpleM gridL flowTable tryLCbd
exec_default := $(addprefix ../bin/, $(executables))
exec_extras := $(addprefix ../bin/, $(extras))

ice_sources := extrasGSL.cc grid.cc materials.cc nc_util.cc beddefLC.cc		\
	forcing.cc iMadaptive.cc iMbasal.cc iMbeddef.cc iMbootstrap.cc		\
	iMdefaults.cc iMforcing.cc iMgeometry.cc iMgrainsize.cc iMinverse.cc    \
	iMinverseMat.cc iMIO.cc iMmatlab.cc iMnames.cc iMoptions.cc iMpdd.cc    \
	iMreport.cc iMssa.cc iMsia.cc iMtemp.cc iMtests.cc iMutil.cc iMvelocity.cc \
	iMviewers.cc iceModelVec.cc iceModelVec2.cc iceModelVec3.cc iceModel.cc	\
	pism_const.cc LocalInterpCtx.cc

ice_csources := cubature.c pism_signal.c

tests_sources := exactTestsABCDE.c exactTestsFG.c exactTestH.c exactTestsIJ.c \
	exactTestK.c exactTestL.c exactTestM.c

other_sources := pismr.cc pismd.cc pismv.cc pisms.cc iceEISModel.cc	\
	iceMISMIPModel.cc iceROSSModel.cc icePSTexModel.cc	\
	iceCompModel.cc iceExactSSAModel.cc iceCalvBCModel.cc iCMthermo.cc flowTable.cc		\
	tryLCbd.cc
other_csources := simpleABCD.c simpleE.c simpleFG.c simpleH.c simpleI.c	\
	simpleJ.c simpleK.c simpleL.c simpleM.c

# this is the main derived class example so it is structured like an example
#  derived class, including this make piece
include ../src/eismint/grn_make

# Import additional make include files:
ifdef CONFIG
  include $(addprefix ../,$(CONFIG))
endif

all_sources := $(ice_sources) $(ice_csources) $(test_sources) $(other_sources)	\
		$(other_csources)

TESTS_OBJS := $(tests_sources:.c=.o)

ICE_OBJS := $(ice_sources:.cc=.o) $(ice_csources:.c=.o)

depfiles := $(patsubst %.c,%.d,$(patsubst %.cc,%.d,$(all_sources)))


all: ${LIBPISMNAME} ${LIBTESTSNAME} $(exec_default)
	@svnversion ../src/ > ../src/revision

extras : ${exec_extras}

# build libraries; shared by default, static if PISM_STATIC=1.

ifneq ($(PISM_STATIC),1)
LIBPISMNAME := ../lib/libpism.so
LIBTESTSNAME := ../lib/libtests.so
../lib/libpism.so : ${ICE_OBJS}
	${CXXLINKER} $(SHARED) ${ICE_OBJS} -o $@
../lib/libtests.so : ${TESTS_OBJS}
	${CLINKER} $(SHARED) ${TESTS_OBJS} -o $@
else
LIBPISMNAME := ../lib/libpism.a
LIBTESTSNAME := ../lib/libtests.a
../lib/libpism.a : ${ICE_OBJS}
	ar r $@ ${ICE_OBJS} 
../lib/libtests.a : ${TESTS_OBJS}
	ar r $@ ${TESTS_OBJS} 
endif

#build default executables

../bin/pismr : pismr.o ${LIBPISMNAME}
	${CXXLINKER} $< ${ICE_LIB_FLAGS} -o $@

../bin/pismd : pismd.o iceROSSModel.o ${LIBPISMNAME}
	${CXXLINKER} iceROSSModel.o pismd.o ${ICE_LIB_FLAGS} -o $@

../bin/pisms : iceEISModel.o iceMISMIPModel.o icePSTexModel.o pisms.o ${LIBPISMNAME}
	${CXXLINKER} iceEISModel.o iceMISMIPModel.o icePSTexModel.o pisms.o \
	${ICE_LIB_FLAGS} -o $@

../bin/pismv : iCMthermo.o iceCompModel.o iceExactSSAModel.o iceCalvBCModel.o \
 pismv.o ${LIBPISMNAME} ${LIBTESTSNAME}
	${CXXLINKER} iCMthermo.o iceCompModel.o \
	iceExactSSAModel.o iceCalvBCModel.o pismv.o ${ICE_LIB_FLAGS} -o $@

#build extra executables

../bin/flowTable : flowTable.o materials.o ${LIBPISMNAME}
	${CXXLINKER} flowTable.o materials.o ${ICE_LIB_FLAGS} -o $@

../bin/tryLCbd : tryLCbd.o beddefLC.o materials.o ${LIBPISMNAME}
	${CXXLINKER} tryLCbd.o beddefLC.o materials.o ${ICE_LIB_FLAGS} -o $@

../bin/simpleABCD : simpleABCD.o ${LIBTESTSNAME}
	${CLINKER} $< ${TESTS_LIB_FLAGS} -o $@

../bin/simpleE : simpleE.o ${LIBTESTSNAME}
	${CLINKER} $< ${TESTS_LIB_FLAGS} -o $@

../bin/simpleFG : simpleFG.o ${LIBTESTSNAME}
	${CLINKER} $< ${TESTS_LIB_FLAGS} -o $@

../bin/simpleH : simpleH.o ${LIBTESTSNAME}
	${CLINKER} $< ${TESTS_LIB_FLAGS} -o $@

../bin/simpleI : simpleI.o ${LIBTESTSNAME}
	${CLINKER} $< ${TESTS_LIB_FLAGS} -o $@

../bin/simpleJ : simpleJ.o ${LIBTESTSNAME}
	${CLINKER} $< ${TESTS_LIB_FLAGS} -o $@

../bin/simpleK : simpleK.o ${LIBTESTSNAME}
	${CLINKER} $< ${TESTS_LIB_FLAGS} -o $@

../bin/simpleL : simpleL.o ${LIBTESTSNAME}
	${CLINKER} $< ${TESTS_LIB_FLAGS} -o $@

../bin/gridL : gridL.o ${LIBTESTSNAME}
	${CLINKER} $< ${TESTS_LIB_FLAGS} -o $@

../bin/simpleM : simpleM.o ${LIBTESTSNAME}
	${CLINKER} $< ${TESTS_LIB_FLAGS} -o $@

# The GNU recommended procedure for automatically generating prerequisites.
# (The obscure sed command comes straight from the GNU Make manual.)

# This rule updates the `*.d' to reflect changes in `*.cc' files
%.d : %.cc
	@echo "Prerequisites from" $< "-->" $@
	@set -e; rm -f $@; \
	 $(CC) -w -c -MM $< > $@.$$$$; \
	 sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	 rm -f $@.$$$$

# This rule updates the `*.d' to reflect changes in `*.c' files
%.d : %.c
	@echo "Prerequisites from" $< "-->" $@
	@set -e; rm -f $@; \
	 $(CC) -w -c -MM $< > $@.$$$$; \
	 sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	 rm -f $@.$$$$

showEnv:
	@echo 'VPATH = ' ${VPATH}
	@echo 'CC = ' ${CC}
	@echo 'CLINKER = ' ${CLINKER}
	@echo 'CXXLINKER = ' ${CXXLINKER}
	@echo 'PETSC_DIR = ' ${PETSC_DIR}
	@echo 'PETSC_ARCH = ' ${PETSC_ARCH}
	@echo 'PETSC_LIB = ' ${PETSC_LIB}
	@echo 'CFLAGS = ' ${CFLAGS}
	@echo 'ICE_LIB_FLAGS = ' ${ICE_LIB_FLAGS}
	@echo 'TESTS_LIB_FLAGS = ' ${TESTS_LIB_FLAGS}
	@echo 'LIBPISMNAME = ' ${LIBPISMNAME}
	@echo 'LIBTESTSNAME = ' ${LIBTESTSNAME}
	@echo 'ICE_OBJS = ' ${ICE_OBJS}
	@echo 'PISM_REVISION = ' "$(PISM_REVISION)"
	@echo 'PISM_HAVE_PETSC3 = ' ${PISM_HAVE_PETSC3}

depend: ;

.PHONY: clean distclean depclean

ifeq (${PISM_HAVE_PETSC3}, 1)
# PETSc 3.0.0 apparently adds a "clean::" target
clean:: depclean
else
clean: depclean
endif

depclean:
	@rm -f $(depfiles)

distclean: clean
	@rm -f ${LIBPISMNAME} ${LIBTESTSNAME} $(exec_default) $(exec_extras)

# This is to avoid recreating .d files if you all you need is 'make clean' or
# 'make distclean':
ifneq ($(MAKECMDGOALS),clean)
  ifneq ($(MAKECMDGOALS),distclean)
    include $(depfiles)
  endif
endif

