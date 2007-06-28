SHELL = /bin/sh
VPATH = src:src/exact
ALL : all

# get PETSc environment, rules:
include ${PETSC_DIR}/bmake/common/base

#FLAGS:

MARGIN_TRICK?=0
MARGIN_TRICK_TWO?=0
WITH_FFTW?=1
WITH_GSL?=1
CFLAGS+= -DWITH_FFTW=${WITH_FFTW} -DWITH_GSL=${WITH_GSL}\
   -DMARGIN_TRICK=${MARGIN_TRICK} -DMARGIN_TRICK_TWO=${MARGIN_TRICK_TWO} -pipe

ICE_LIB_FLAGS= -L`pwd`/lib -Wl,-rpath,`pwd`/lib -lpism -ltests ${PETSC_LIB}\
   -lnetcdf_c++ -lnetcdf
ifeq (${WITH_FFTW}, 1)
	ICE_LIB_FLAGS+= -lfftw3
endif
ifeq (${WITH_GSL}, 1)
	ICE_LIB_FLAGS+= -lgsl -lgslcblas
endif

#VARIABLES:

executables= pismr pismv pisms pant
extra_execs= pgrn simpleISO simpleFG simpleI flowTable

ice_sources= extrasGSL.cc grid.cc iMbasal.cc iMbeddef.cc iMdefaults.cc\
	iMgrainsize.cc iMIO.cc iMIOnetcdf.cc iMmacayeal.cc iMoptions.cc\
	iMregrid.cc iMregrid_netCDF.cc iMtemp.cc iMutil.cc iMvelocity.cc\
	 iMviewers.cc iceModel.cc materials.cc nc_util.cc
ice_csources= cubature.c pism_signal.c
ICE_OBJS= $(ice_sources:.cc=.o) $(ice_csources:.c=.o)

tests_sources= exactTestsABCDE.c exactTestsFG.c exactTestH.c exactTestI.c
TESTS_OBJS= $(tests_sources:.c=.o)

other_sources= pismr.cc pismv.cc pisms.cc pant.cc\
	flowTable.cc iceEISModel.cc iceHEINOModel.cc\
	iceROSSModel.cc iceCompModel.cc shelf.cc iceGRNModel.cc pgrn.cc
other_csources= simpleISO.c simpleFG.c simpleI.c

depfiles= $(ice_sources:.cc=.d) $(ice_csources:.c=.d) $(tests_sources:.c=.d)\
	$(other_sources:.cc=.d) $(other_csources:.c=.d)

all : depend libpism libtests $(executables)

libpism : ${ICE_OBJS}
	${CLINKER} -shared -o lib/libpism.so ${ICE_OBJS}
#for static-linking:  ar cru -s lib/libpism.a ${ICE_OBJS}   etc

libtests : ${TESTS_OBJS}
	${CLINKER} -shared -o lib/libtests.so ${TESTS_OBJS}

pismr : pismr.o lib/libpism.so
	${CLINKER} $< ${ICE_LIB_FLAGS} -o bin/pismr

pisms : iceEISModel.o iceHEINOModel.o iceROSSModel.o pisms.o lib/libpism.so
	${CLINKER} iceEISModel.o iceHEINOModel.o iceROSSModel.o pisms.o ${ICE_LIB_FLAGS} -o bin/pisms

pismv : iceCompModel.o iceExactStreamModel.o pismv.o lib/libpism.so lib/libtests.so
	${CLINKER} iceCompModel.o iceExactStreamModel.o pismv.o ${ICE_LIB_FLAGS} -o bin/pismv

pant : pant.o lib/libpism.so
	${CLINKER} $< ${ICE_LIB_FLAGS} -o bin/pant

pgrn : iceGRNModel.o pgrn.o lib/libpism.so
	${CLINKER} iceGRNModel.o pgrn.o ${ICE_LIB_FLAGS} -o bin/pgrn

#shelf : shelf.o lib/libpism.so
#	${CLINKER} $< ${ICE_LIB_FLAGS} -o bin/shelf

flowTable : flowTable.o lib/libpism.so
	${CLINKER} $< ${ICE_LIB_FLAGS} -o bin/flowTable

simpleISO : simpleISO.o lib/libtests.so
	${CLINKER} $< -lm -L`pwd`/lib -Wl,-rpath,`pwd`/lib -ltests -o bin/simpleISO

simpleFG : simpleFG.o lib/libtests.so
	${CLINKER} $< -lm -L`pwd`/lib -Wl,-rpath,`pwd`/lib -ltests -o bin/simpleFG

simpleI : simpleI.o lib/libtests.so
	${CLINKER} $< -lm -L`pwd`/lib -Wl,-rpath,`pwd`/lib -ltests -o bin/simpleI

# Cancel the implicit rules
% : %.cc
% : %.c

# Emacs style tags
.PHONY: tags TAGS
tags TAGS :
	etags *.cc *.hh *.c

# The GNU recommended proceedure for automatically generating dependencies.
# This rule updates the `*.d' to reflect changes in `*.cc' files
%.d : %.cc
	@echo "Dependencies from" $< "-->" $@
	@set -e; rm -f $@; \
	 $(CC) -MM $< > $@.$$$$; \
	 sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	 rm -f $@.$$$$

# This rule updates the `*.d' to reflect changes in `*.c' files
%.d : %.c
	@echo "Dependencies from" $< "-->" $@
	@set -e; rm -f $@; \
	 $(CC) -MM $< > $@.$$$$; \
	 sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	 rm -f $@.$$$$

depend : $(depfiles)

depclean :
	@rm -f *.d

clean : depclean

distclean : clean
	rm -f TAGS lib/libpism.so lib/libtests.so \
	 $(patsubst %, bin/%, ${executables})\
	 $(patsubst %, bin/%, ${extra_execs})

.PHONY: clean

include $(depfiles)

