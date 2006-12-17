SHELL = /bin/sh
VPATH = src:src/exact
ALL : all

#FLAGS:

WITH_NETCDF?=1
WITH_FFTW?=1
WITH_GSL?=1
CFLAGS+= -DWITH_NETCDF=${WITH_NETCDF} -DWITH_FFTW=${WITH_FFTW}\
	-DWITH_GSL=${WITH_GSL} -pipe

ICE_LIB_FLAGS= -L`pwd`/obj -Wl,-rpath,`pwd`/obj -lpism -ltests ${PETSC_LIB}
ifeq (${WITH_NETCDF}, 1)
	ICE_LIB_FLAGS+= -lnetcdf_c++ -lnetcdf
endif
ifeq (${WITH_FFTW}, 1)
	ICE_LIB_FLAGS+= -lfftw3
endif
ifeq (${WITH_GSL}, 1)
	ICE_LIB_FLAGS+= -lgsl -lgslcblas
endif

#VARIABLES:

executables= flowTable pismr pismv pisms simpleISO simpleFG simpleI shelf get_drag

ice_sources= extrasGSL.cc grid.cc iMbasal.cc iMbeddef.cc iMdefaults.cc\
	iMgrainsize.cc iMIO.cc iMIOnetcdf.cc iMmacayeal.cc iMoptions.cc\
	iMregrid.cc iMtemp.cc iMutil.cc iMvelocity.cc iMviewers.cc\
	iceModel.cc materials.cc 
ice_csources= cubature.c
ICE_OBJS= $(ice_sources:.cc=.o) cubature.o

tests_sources= exactTestsABCDE.c exactTestsFG.c exactTestH.c exactTestI.c
TESTS_OBJS= $(tests_sources:.c=.o)

other_sources= flowTable.cc simplify.cc iceEISModel.cc iceHEINOModel.cc\
	iceROSSModel.cc run.cc verify.cc iceCompModel.cc get_drag.cc shelf.cc
other_csources= simpleISO.c simpleFG.c simpleI.c

depfiles= $(ice_sources:.cc=.d) $(ice_csources:.c=.d) $(tests_sources:.c=.d)\
	$(other_sources:.cc=.d) $(other_csources:.c=.d)

include ${PETSC_DIR}/bmake/common/base

#TARGETS:

all : depend libpism libtests $(executables)

libpism : ${ICE_OBJS}
	${CLINKER} -shared -o obj/libpism.so ${ICE_OBJS}
#	ar cru -s obj/libpism.a ${ICE_OBJS}

libtests : ${TESTS_OBJS}
	${CLINKER} -shared -o obj/libtests.so ${TESTS_OBJS}

flowTable : obj/libpism.so flowTable.o
	${CLINKER} $^ ${ICE_LIB_FLAGS} -o obj/flowTable

get_drag : obj/libpism.so get_drag.o
	${CLINKER} $^ ${ICE_LIB_FLAGS} -o obj/get_drag

pismr : obj/libpism.so run.o
	${CLINKER} $^ ${ICE_LIB_FLAGS} -o obj/pismr

pisms : obj/libpism.so iceEISModel.o iceHEINOModel.o iceROSSModel.o simplify.o
	${CLINKER} $^ ${ICE_LIB_FLAGS} -o obj/pisms

pismv : obj/libpism.so obj/libtests.so iceCompModel.o iceExactStream.o verify.o
	${CLINKER} $^ ${ICE_LIB_FLAGS} -o obj/pismv

shelf : obj/libpism.so shelf.o
	${CLINKER} $^ ${ICE_LIB_FLAGS} -o obj/shelf

simpleISO : obj/libtests.so simpleISO.o
	${CLINKER} $^ -lm -L`pwd`/obj -Wl,-rpath,`pwd`/obj -ltests \
	 -o obj/simpleISO

simpleFG : obj/libtests.so simpleFG.o
	${CLINKER} $^ -lm -L`pwd`/obj -Wl,-rpath,`pwd`/obj -ltests \
	 -o obj/simpleFG

simpleI : obj/libtests.so simpleI.o
	${CLINKER} $^ -lm -L`pwd`/obj -Wl,-rpath,`pwd`/obj -ltests \
	 -o obj/simpleI

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
	 $(PETSC_COMPILE_SINGLE) -MM $< > $@.$$$$; \
	 sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	 rm -f $@.$$$$

# This rule updates the `*.d' to reflect changes in `*.c' files
%.d : %.c
	@echo "Dependencies from" $< "-->" $@
	@set -e; rm -f $@; \
	 $(PETSC_COMPILE_SINGLE) -MM $< > $@.$$$$; \
	 sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	 rm -f $@.$$$$

depend : $(depfiles)

depclean :
	@rm -f *.d

clean : depclean

distclean : clean
	rm -f TAGS obj/libpism.so obj/libtests.so \
	 $(patsubst %, obj/%, ${executables})

.PHONY: clean

include $(depfiles)
