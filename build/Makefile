ifndef PISM_PREFIX
  $(error Please run make from the top-level PISM directory)
endif

# Get the current PISM revision:
url = $(shell svn info | grep URL)
rev = $(shell svnversion ../)
ifeq ($(findstring trunk,$(url)),trunk)
  PISM_REVISION = revision $(rev) trunk
  # PISM_REVISION = trunk: revision $(rev)
endif
ifeq ($(findstring tags,$(url)),tags)
  PISM_REVISION = revision $(rev) tags $(lastword $(subst /, ,$(url)))
  # PISM_REVISION = tags $(lastword $(subst /, ,$(url))): revision $(rev)
endif
ifndef PISM_REVISION
  $(error Can not determine revision number.  Please define PISM_REVISION in the top-level Makefile)
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
CFLAGS += -DPISM_HAVE_FFTW=${PISM_HAVE_FFTW} -DPISM_REVISION='"$(PISM_REVISION)"' \
          -I../src/udunits -I../src/base -DPISM_DEFAULT_CONFIG_FILE='"$(PISM_PREFIX)/lib/pism_config.nc"'
ifeq (${PISM_USE_GNU_FLAGS}, 1)
  CFLAGS += -Woverloaded-virtual -pipe -O2
endif
ifeq (${PISM_USE_DEBUG}, 1)
  CFLAGS += -g -Wall -Wextra
endif

# Flags needed by UDUNITS:
UDUNITS_CFLAGS   = -c -O2 -fPIC
UDUNITS_CPPFLAGS = -DUT_INSTALL_PATH='"$(PISM_PREFIX)/lib/pismudunits.dat"' \
	-DUT_SOURCE_PATH='"$(PISM_PREFIX)/lib/pismudunits.dat"' -DNDEBUG -Df2cFortran

# Flags related to building libraries; defaults are for GNU/Linux
SHARED = -shared

# Linker flags:
BASE_LIB_FLAGS := ${PISM_EXTRA_LIBS} \
                  -L../lib -Wl,-rpath,$(PISM_PREFIX)/lib \
                  -lm -lgsl -lgslcblas
ifeq (${PISM_HAVE_FFTW}, 1)
  BASE_LIB_FLAGS += -lfftw3
endif
TESTS_LIB_FLAGS := -ltests ${BASE_LIB_FLAGS}
ICE_LIB_FLAGS := -lpism -lnetcdf -lpismudunits ${PETSC_LIB} ${BASE_LIB_FLAGS}

# set names for libraries
ifneq ($(PISM_STATIC),1)
LIBPISMNAME := ../lib/libpism.so
LIBTESTSNAME := ../lib/libtests.so
LIBUDUNITS = ../lib/libpismudunits.so
else
LIBPISMNAME := ../lib/libpism.a
LIBTESTSNAME := ../lib/libtests.a
LIBUDUNITS = ../lib/libpismudunits.a
endif

ifeq ($(PISM_USE_MPICXX), 1)
  # Override PETSc's choice of CXX and of linker
  CXX = $(subst mpicc,mpicxx,$(CC))
  CXXLINKER = $(subst mpicc,mpicxx,$(CLINKER))
endif

SHELL := /bin/sh
VPATH := $(subst src,../src,src/:src/base/:src/earth/:src/verif/:src/verif/tests/:src/software_tests:src/eismint/:src/ismip/:src/udunits/:src/coupler/):../lib

# Variables:
executables := penth pismr pismd pismv pisms pcctest pismtests
extras := flowTable tryLCbd
exec_default := $(addprefix ../bin/, $(executables))
exec_extras := $(addprefix ../bin/, $(extras))
config_file := ../lib/pism_config.nc

# ice_sources and ice_csources go in libpism.[so|a]
ice_sources := grid.cc LocalInterpCtx.cc materials.cc nc_util.cc \
	pism_default_config.cc pism_revision.cc columnSystem.cc \
	pism_const.cc iceModelVec.cc iceModelVec2.cc iceModelVec3.cc NCVariable.cc \
	matlablike.cc greens.cc deformation.cc \
	localMassBalance.cc monthlyDataMaps.cc pccoupler.cc pGreenlandAtmosCoupler.cc \
	iceModel.cc iMadaptive.cc iMbasal.cc iMbeddef.cc iMbootstrap.cc \
	iMdefaults.cc iMgeometry.cc iMgrainsize.cc iMinit.cc iMinverse.cc \
	iMinverseMat.cc iMIO.cc iMmatlab.cc iMoptions.cc \
	iMreport.cc iMssa.cc iMssaSNES.cc iMsia.cc iMtemp.cc \
	iMutil.cc iMvelocity.cc iMviewers.cc Timeseries.cc PISMVars.cc iMtimeseries.cc \
	enthalpyConverter.cc
ice_csources := cubature.c pism_signal.c

# tests_sources go in libtests.[so|a]
tests_sources := exactTestsABCDE.c exactTestsFG.c exactTestH.c exactTestsIJ.c \
	exactTestK.c exactTestL.c exactTestM.c

other_sources := penth.cc pismr.cc pismd.cc pismv.cc pisms.cc \
	iceEISModel.cc iceMISMIPModel.cc iceROSSModel.cc icePSTexModel.cc \
	iceCompModel.cc iCMthermo.cc iceExactSSAModel.cc iceCalvBCModel.cc \
	iceEnthalpyModel.cc enthColumnSystem.cc \
	pcctest.cc flowTable.cc tryLCbd.cc iMtests.cc pismtests.cc

# this is the primary derived class example so it is structured as an example
#   derived class; as an example it includes "grn_make"; creates executable
#   "pgrn"
include ../src/eismint/grn_make

# Import additional make include files:
ifdef CONFIG
  include $(addprefix ../,$(CONFIG))
endif

all_sources := $(ice_sources) $(ice_csources) $(test_sources) $(other_sources)	\
		$(other_csources)

TESTS_OBJS := $(addprefix ../src/verif/tests/, $(tests_sources:.c=.o))

ICE_OBJS := $(ice_sources:.cc=.o) $(ice_csources:.c=.o)

depfiles := $(patsubst %.c,%.d,$(patsubst %.cc,%.d,$(all_sources)))

all: ${LIBPISMNAME} ${LIBTESTSNAME} $(LIBUDUNITS) $(exec_default) $(config_file)
	@svnversion ../src/ > ../src/revision

extras : ${exec_extras}

# how to build libraries; shared by default, static if PISM_STATIC=1.
ifneq ($(PISM_STATIC),1)
  BUILD_LIB := ${CXXLINKER} $(SHARED) -o
else
  BUILD_LIB := ar r
endif

# $@ should go first here:
$(LIBPISMNAME)  : ${ICE_OBJS}
	$(BUILD_LIB) $@ $^
$(LIBTESTSNAME) : ${TESTS_OBJS}
	$(BUILD_LIB) $@ $^

# Build UDUNITS:
ud_objects := utlib.o utparse.o utscan.o udalloc.o
# This ensures that UDUNITS is built using its special flags:
$(ud_objects): %.o: %.c
	$(CC) $(UDUNITS_CPPFLAGS) $(UDUNITS_CFLAGS) $< -o $@

$(LIBUDUNITS): $(ud_objects)
	$(BUILD_LIB) $@ $^
	cp ../src/udunits/pismudunits.dat ../lib/

#build default executables

$(config_file): pism_config.cdl
	ncgen -o $@ $^

$(exec_default) : $(LIBPISMNAME) $(LIBUDUNITS)

../bin/penth : penth.o iceEnthalpyModel.o enthColumnSystem.o enthalpyConverter.o ${LIBPISMNAME}
	${CXXLINKER} penth.o iceEnthalpyModel.o enthColumnSystem.o enthalpyConverter.o ${ICE_LIB_FLAGS} -o $@

../bin/pismr : pismr.o ${LIBPISMNAME}
	${CXXLINKER} pismr.o ${ICE_LIB_FLAGS} -o $@

../bin/pismd : pismd.o iceROSSModel.o ${LIBPISMNAME}
	${CXXLINKER} iceROSSModel.o pismd.o ${ICE_LIB_FLAGS} -o $@

../bin/pisms : iceEISModel.o iceMISMIPModel.o icePSTexModel.o pisms.o ${LIBPISMNAME}
	${CXXLINKER} iceEISModel.o iceMISMIPModel.o icePSTexModel.o pisms.o \
	${ICE_LIB_FLAGS} -o $@

../bin/pismv : iCMthermo.o iceCompModel.o iceExactSSAModel.o iceCalvBCModel.o \
 pismv.o ${LIBTESTSNAME} ${LIBPISMNAME}
	${CXXLINKER} iCMthermo.o iceCompModel.o iceExactSSAModel.o \
	iceCalvBCModel.o pismv.o ${ICE_LIB_FLAGS} ${TESTS_LIB_FLAGS} -o $@

../bin/pcctest : pcctest.o ${LIBPISMNAME}
	${CXXLINKER} pcctest.o ${ICE_LIB_FLAGS} -o $@

../bin/pismtests: pismtests.o iMtests.o ${LIBPISMNAME}
	${CXXLINKER} pismtests.o iMtests.o ${ICE_LIB_FLAGS} -o $@


#build extra executables

../bin/flowTable : flowTable.o materials.o ${LIBPISMNAME}
	${CXXLINKER} flowTable.o materials.o ${ICE_LIB_FLAGS} -o $@

../bin/tryLCbd : tryLCbd.o deformation.o greens.o materials.o pism_const.o ${LIBPISMNAME}
	${CXXLINKER} tryLCbd.o deformation.o greens.o materials.o pism_const.o ${ICE_LIB_FLAGS} -o $@

# The GNU recommended procedure for automatically generating prerequisites.
# (The obscure sed command comes straight from the GNU Make manual.)

# This rule updates the `*.d' to reflect changes in `*.cc' files
%.d : %.cc
	@echo "Prerequisites from" $< "-->" $@
	@set -e; rm -f $@; \
	 $(CXX) $(CFLAGS) -w -c -MM $< > $@.$$$$; \
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
	@echo 'ICE_OBJS = ' ${ICE_OBJS}
	@echo 'ICE_LIB_FLAGS = ' ${ICE_LIB_FLAGS}
	@echo 'TESTS_LIB_FLAGS = ' ${TESTS_LIB_FLAGS}
	@echo 'PETSC_LIB = ' ${PETSC_LIB}
	@echo 'CFLAGS = ' ${CFLAGS}
	@echo 'VPATH = ' ${VPATH}
	@echo 'CC = ' ${CC}
	@echo 'CXX = ' ${CXX}
	@echo 'CFLAGS = ' ${CFLAGS}
	@echo 'CLINKER = ' ${CLINKER}
	@echo 'CXXLINKER = ' ${CXXLINKER}
	@echo 'PETSC_DIR = ' ${PETSC_DIR}
	@echo 'PETSC_ARCH = ' ${PETSC_ARCH}
	@echo 'LIBPISMNAME = ' ${LIBPISMNAME}
	@echo 'LIBTESTSNAME = ' ${LIBTESTSNAME}
	@echo 'PISM_REVISION = ' ${PISM_REVISION}
	@echo 'PISM_PREFIX = ' ${PISM_PREFIX}
	@echo 'PISM_HAVE_PETSC3 = ' ${PISM_HAVE_PETSC3}
	@echo 'PISM_HAVE_FFTW = ' ${PISM_HAVE_FFTW}
	@echo 'PISM_STATIC = ' ${PISM_STATIC}
	@echo 'PISM_USE_DEBUG = ' ${PISM_USE_DEBUG}
	@echo 'PISM_USE_MPICXX = ' ${PISM_USE_MPICXX}
	@echo 'PISM_USE_GNU_FLAGS = ' ${PISM_USE_GNU_FLAGS}
	@echo 'PISM_EXTRA_LIBS = ' ${PISM_EXTRA_LIBS}
	@echo 'CONFIG = ' ${CONFIG}
	@echo 'UDUNITS_CFLAGS = ' $(UDUNITS_CFLAGS)
	@echo 'UDUNITS_CPPFLAGS = ' $(UDUNITS_CPPFLAGS)

depend: ;

.PHONY: clean distclean depclean

ifeq (${PISM_HAVE_PETSC3}, 1)
# PETSc 3.0.0 apparently adds a "clean::" target
clean:: depclean
else
clean: depclean
endif

depclean:
	@rm -f $(depfiles) src/revision

distclean: clean
	@rm -f ${LIBPISMNAME} ${LIBTESTSNAME} $(LIBUDUNITS) $(exec_default) $(exec_extras) \
            ../lib/pismudunits.dat

# This is to avoid recreating .d files if you all you need is 'make clean' or
# 'make distclean':
ifneq ($(MAKECMDGOALS),clean)
  ifneq ($(MAKECMDGOALS),distclean)
    include $(depfiles)
  endif
endif

