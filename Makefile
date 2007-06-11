SHELL = /bin/sh
VPATH = src:src/exact
ALL : all

include ${PETSC_DIR}/bmake/common/base

#FLAGS:

WITH_NETCDF?=1
WITH_FFTW?=1
WITH_GSL?=1
CFLAGS+= -DWITH_NETCDF=${WITH_NETCDF} -DWITH_FFTW=${WITH_FFTW} -DWITH_GSL=${WITH_GSL} -pipe

ICE_LIB_FLAGS= -L`pwd`/obj -Wl,-rpath,`pwd`/obj -lpism -ltests ${PETSC_LIB}\
   -lnetcdf_c++ -lnetcdf
ifeq (${WITH_FFTW}, 1)
	ICE_LIB_FLAGS+= -lfftw3
endif
ifeq (${WITH_GSL}, 1)
	ICE_LIB_FLAGS+= -lgsl -lgslcblas
endif

#VARIABLES:

executables= pismr pismv pisms pant

ice_sources= extrasGSL.cc grid.cc iMbasal.cc iMbeddef.cc iMdefaults.cc\
	iMgrainsize.cc iMIO.cc iMIOnetcdf.cc iMmacayeal.cc iMoptions.cc\
	iMregrid.cc iMregrid_netCDF.cc iMtemp.cc iMutil.cc iMvelocity.cc\
	 iMviewers.cc iceModel.cc materials.cc nc_util.cc
ice_csources= cubature.c pism_signal.c
ICE_OBJS= $(ice_sources:.cc=.o) $(ice_csources:.c=.o)

tests_sources= exactTestsABCDE.c exactTestsFG.c exactTestH.c exactTestI.c
TESTS_OBJS= $(tests_sources:.c=.o)

other_sources= flowTable.cc simplify.cc iceEISModel.cc iceHEINOModel.cc\
	iceROSSModel.cc run.cc verify.cc iceCompModel.cc shelf.cc pant.cc
other_csources= simpleISO.c simpleFG.c simpleI.c

depfiles= $(ice_sources:.cc=.d) $(ice_csources:.c=.d) $(tests_sources:.c=.d)\
	$(other_sources:.cc=.d) $(other_csources:.c=.d)

all : depend libpism libtests $(executables)

# no longer works correctly with petsc 2.3.3:
#CCLINKER=`echo ${CLINKER} | sed 's/mpicc/mpicxx/'`

libpism : ${ICE_OBJS}
	${CLINKER} -shared -o obj/libpism.so ${ICE_OBJS}
#	ar cru -s obj/libpism.a ${ICE_OBJS}

libtests : ${TESTS_OBJS}
	${CLINKER} -shared -o obj/libtests.so ${TESTS_OBJS}

pismr : run.o obj/libpism.so
	${CLINKER} $< ${ICE_LIB_FLAGS} -o bin/pismr

pisms : iceEISModel.o iceHEINOModel.o iceROSSModel.o simplify.o obj/libpism.so
	${CLINKER} iceEISModel.o iceHEINOModel.o iceROSSModel.o simplify.o ${ICE_LIB_FLAGS} -o bin/pisms

pismv : iceCompModel.o iceExactStreamModel.o verify.o obj/libpism.so obj/libtests.so
	${CLINKER} iceCompModel.o iceExactStreamModel.o verify.o ${ICE_LIB_FLAGS} -o bin/pismv

pant : pant.o obj/libpism.so
	${CLINKER} $< ${ICE_LIB_FLAGS} -o bin/pant

#shelf : shelf.o obj/libpism.so
#	${CLINKER} $< ${ICE_LIB_FLAGS} -o bin/shelf

flowTable : flowTable.o obj/libpism.so
	${CLINKER} $< ${ICE_LIB_FLAGS} -o bin/flowTable

simpleISO : simpleISO.o obj/libtests.so
	${CLINKER} $< -lm -L`pwd`/obj -Wl,-rpath,`pwd`/obj -ltests -o bin/simpleISO

simpleFG : simpleFG.o obj/libtests.so
	${CLINKER} $< -lm -L`pwd`/obj -Wl,-rpath,`pwd`/obj -ltests -o bin/simpleFG

simpleI : simpleI.o obj/libtests.so
	${CLINKER} $< -lm -L`pwd`/obj -Wl,-rpath,`pwd`/obj -ltests -o bin/simpleI

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
	rm -f TAGS obj/libpism.so obj/libtests.so \
	 $(patsubst %, bin/%, ${executables})

.PHONY: clean

include $(depfiles)

