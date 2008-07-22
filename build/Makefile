ifndef PISM_PREFIX
  $(error Please run make from the top-level PISM directory)
endif

# Get the current PISM revision:
url = $(shell svn info | grep URL)
rev = $(shell svnversion ../)
ifeq ($(findstring trunk,$(url)),trunk)
  PISM_REVISION = (trunk: revision $(rev))
endif
ifeq ($(findstring branches,$(url)),branches)
  PISM_REVISION = (branch $(lastword $(subst /, ,$(url))): revision $(rev))
endif
ifndef PISM_REVISION
  $(error Can't determine PISM's revision number. \
  Please define PISM_REVISION in the top-level Makefile)
endif

DISPLAY = :0.0

# Get PETSc environment and rules:
include $(PETSC_DIR)/bmake/common/base

# Compiler flags:
CFLAGS += -DWITH_FFTW=${WITH_FFTW} -DLOG_PISM_EVENTS=${LOG_PISM_EVENTS} \
          -DMISMIP_PLAY=${MISMIP_PLAY} -DPISM_REVISION='"$(PISM_REVISION)"' \
          -pipe

SHARED = -shared # default (for Linux)

# Linker flags:
TESTS_LIB_FLAGS := -L. -Wl,-rpath,$(PISM_PREFIX) -Wl,-rpath,$(PISM_PREFIX)/lib \
                   -ltests -lm -lgsl -lgslcblas

ICE_LIB_FLAGS := -lpism ${TESTS_LIB_FLAGS} ${PETSC_LIB} -lnetcdf 
ifeq (${WITH_FFTW}, 1)
  ICE_LIB_FLAGS += -lfftw3
endif

# Override PETSc's choice of the linker (if asked to):
ifeq ($(USE_MPICXX), 1)
  CXXLINKER = $(subst mpicc,mpicxx,$(CLINKER))
endif

SHELL := /bin/sh
VPATH := $(subst src,../src,src/:src/base/:src/num/:src/verif/:src/eismint/:src/ismip/)

# Variables:
executables := pismr pismd pismv pisms pgrn
extra_execs := simpleABCD simpleE simpleFG simpleH simpleI simpleJ \
   simpleL gridL flowTable tryLCbd pant

ice_sources := extrasGSL.cc grid.cc materials.cc nc_util.cc beddefLC.cc \
	forcing.cc iMadaptive.cc iMbasal.cc iMbeddef.cc iMbootstrap.cc \
	iMdefaults.cc iMforcing.cc iMgrainsize.cc iMIO.cc iMinverse.cc \
	iMmatlab.cc iMnames.cc iMoptions.cc iMpdd.cc iMreport.cc iMssa.cc \
	iMsia.cc iMtemp.cc iMtests.cc iMutil.cc iMvelocity.cc iMviewers.cc \
	iceModelVec.cc iceModelVec3.cc iceModel.cc pism_const.cc

ice_csources := cubature.c pism_signal.c

tests_sources := exactTestsABCDE.c exactTestsFG.c exactTestH.c exactTestsIJ.c \
   exactTestK.c exactTestL.c

other_sources := pismr.cc pismd.cc pismv.cc pisms.cc pgrn.cc \
	iceEISModel.cc iceMISMIPModel.cc iceROSSModel.cc iceGRNModel.cc \
	icePSTexModel.cc iceCompModel.cc \
	iceExactSSAModel.cc iCMthermo.cc flowTable.cc tryLCbd.cc
other_csources := simpleABCD.c simpleE.c simpleFG.c simpleH.c simpleI.c \
	simpleJ.c simpleK.c simpleL.c
all_sources :=  $(ice_sources) $(ice_csources) $(test_sources) \
		$(other_sources) $(other_csources)

TESTS_OBJS := $(tests_sources:.c=.o)

ICE_OBJS := $(ice_sources:.cc=.o) $(ice_csources:.c=.o)

depfiles := $(patsubst %.c,%.d,$(patsubst %.cc,%.d,$(all_sources)))

# Import additional make include files:
ifdef CONFIG
  include $(addprefix ../,$(CONFIG))
endif

all: libpism.so libtests.so $(executables)
	@svnversion ../src/ > ../src/revision
ifeq ($(findstring install,$(GOALS)),)
	@echo "All the PISM libraries and executables are in $(PWD)."
	@echo "Please don't forget to run 'make install'."
endif

local_install: libpism.so libtests.so $(executables)
	@cp libpism.so ${PISM_PREFIX}/lib/
	@cp libtests.so ${PISM_PREFIX}/lib/
	@echo 'PISM libraries installed in ' ${PISM_PREFIX}'/lib/'
	@cp $(executables) ${PISM_PREFIX}/bin/
	@rm -f libpism.so
	@rm -f libtests.so
	@rm -f $(executables)
	@echo 'PISM executables installed in ' ${PISM_PREFIX}'/bin/'
	@echo 'Note: Object files created during the build process were not deleted.'
	@echo "      Please run 'make clean' if space is an issue."

libpism.so : ${ICE_OBJS}
	${CXXLINKER} $(SHARED) ${ICE_OBJS} -o $@
#for static-linking:  ar cru -s libpism.a ${ICE_OBJS}   etc

libtests.so : ${TESTS_OBJS}
	${CLINKER} $(SHARED) ${TESTS_OBJS} -o $@

pismr : pismr.o libpism.so
	${CXXLINKER} $< ${ICE_LIB_FLAGS} -o $@

pismd : pismd.o iceROSSModel.o libpism.so
	${CXXLINKER} iceROSSModel.o pismd.o ${ICE_LIB_FLAGS} -o $@

pisms : iceEISModel.o iceMISMIPModel.o icePSTexModel.o pisms.o libpism.so
	${CXXLINKER} iceEISModel.o iceMISMIPModel.o icePSTexModel.o pisms.o \
	${ICE_LIB_FLAGS} -o $@

pismv : iCMthermo.o iceCompModel.o iceExactSSAModel.o pismv.o libpism.so libtests.so
	${CXXLINKER} iCMthermo.o iceCompModel.o \
	iceExactSSAModel.o pismv.o ${ICE_LIB_FLAGS} -o $@

pgrn : iceGRNModel.o pgrn.o libpism.so
	${CXXLINKER} iceGRNModel.o pgrn.o ${ICE_LIB_FLAGS} -o $@

flowTable : flowTable.o materials.o
	${CXXLINKER} flowTable.o materials.o ${ICE_LIB_FLAGS} -o $@

tryLCbd : tryLCbd.o beddefLC.o materials.o
	${CXXLINKER} tryLCbd.o beddefLC.o materials.o ${ICE_LIB_FLAGS} -o $@

simpleABCD : simpleABCD.o libtests.so
	${CLINKER} $< ${TESTS_LIB_FLAGS} -o $@

simpleE : simpleE.o libtests.so
	${CLINKER} $< ${TESTS_LIB_FLAGS} -o $@

simpleFG : simpleFG.o libtests.so
	${CLINKER} $< ${TESTS_LIB_FLAGS} -o $@

simpleH : simpleH.o libtests.so
	${CLINKER} $< ${TESTS_LIB_FLAGS} -o $@

simpleI : simpleI.o libtests.so
	${CLINKER} $< ${TESTS_LIB_FLAGS} -o $@

simpleJ : simpleJ.o libtests.so
	${CLINKER} $< ${TESTS_LIB_FLAGS} -o $@

simpleL : simpleL.o libtests.so
	${CLINKER} $< ${TESTS_LIB_FLAGS} -o $@

simpleK : simpleK.o libtests.so
	${CLINKER} $< ${TESTS_LIB_FLAGS} -o $@

gridL : gridL.o libtests.so
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
	@echo 'ICE_OBJS = ' ${ICE_OBJS}
	@echo 'PISM_REVISION = ' "$(PISM_REVISION)"

depend: ;

depclean:
	@rm -f $(depfiles)

clean: depclean

.PHONY: clean distclean depclean

distclean: clean
	@rm -f \
	libpism.so libtests.so \
	$(PISM_PREFIX)/lib/libpism.so $(PISM_PREFIX)/lib/libtests.so \
	$(executables) $(extra_execs) \
	$(addprefix $(PISM_PREFIX)/bin/, $(executables) $(extra_execs))

# This is to avoid recreating .d files if you all you need is 'make clean' or
# 'make distclean':
ifneq ($(MAKECMDGOALS),clean)
  ifneq ($(MAKECMDGOALS),distclean)
    include $(depfiles)
  endif
endif

