SHELL = /bin/sh
VPATH = src:src/exact
ALL : all

# get PETSc environment, rules:
include ${PETSC_DIR}/bmake/common/base


#FLAGS:
MARGIN_TRICK?=0
MARGIN_TRICK_TWO?=0
WITH_FFTW?=1
CFLAGS+= -DWITH_FFTW=${WITH_FFTW} \
   -DMARGIN_TRICK=${MARGIN_TRICK} -DMARGIN_TRICK_TWO=${MARGIN_TRICK_TWO} -pipe

TESTS_LIB_FLAGS= -L`pwd` -L`pwd`/lib -Wl,-rpath,`pwd` -Wl,-rpath,`pwd`/lib -lm\
   -lgsl -lgslcblas -ltests
ICE_LIB_FLAGS= ${TESTS_LIB_FLAGS} -lnetcdf_c++ -lnetcdf -lpism ${PETSC_LIB}
ifeq (${WITH_FFTW}, 1)
	ICE_LIB_FLAGS+= -lfftw3
endif


#VARIABLES:
executables= pismr pismv pisms pgrn pant
extra_execs= simpleABCD simpleE simpleFG simpleH simpleI simpleL gridL flowTable

ice_sources= extrasGSL.cc grid.cc iMbasal.cc iMbeddef.cc iMdefaults.cc\
	iMgrainsize.cc iMIO.cc iMIOnetcdf.cc iMmacayeal.cc iMoptions.cc iMpdd.cc\
	iMregrid.cc iMtemp.cc iMutil.cc iMvelocity.cc\
	iMviewers.cc iceModel.cc materials.cc nc_util.cc
ice_csources= cubature.c pism_signal.c
ICE_OBJS= $(ice_sources:.cc=.o) $(ice_csources:.c=.o)

tests_sources= exactTestsABCDE.c exactTestsFG.c exactTestH.c exactTestI.c exactTestL.c
TESTS_OBJS= $(tests_sources:.c=.o)

other_sources= pismr.cc pismv.cc pisms.cc pant.cc\
	flowTable.cc iceEISModel.cc iceHEINOModel.cc\
	iceROSSModel.cc iceCompModel.cc iCMthermo.cc shelf.cc iceGRNModel.cc pgrn.cc
other_csources= simpleABCD.c simpleE.c simpleFG.c simpleH.c simpleI.c simpleL.c

depfiles= $(ice_sources:.cc=.d) $(ice_csources:.c=.d) $(tests_sources:.c=.d)\
	$(other_sources:.cc=.d) $(other_csources:.c=.d)

all : depend libpism.so libtests.so $(executables) .makeremind

install : local_install clean

local_install : depend libpism.so libtests.so $(executables)
	@cp libpism.so lib/
	@cp libtests.so lib/
	@echo 'PISM libraries installed in lib/'
	@cp $(executables) bin/
	@rm -f libpism.so
	@rm -f libtests.so
	@rm -f $(executables)
	@echo 'PISM executables installed in bin/'
	@rm .makeremind

libpism.so : ${ICE_OBJS}
	${CLINKER} -shared ${ICE_OBJS} -o $@
#for static-linking:  ar cru -s libpism.a ${ICE_OBJS}   etc

libtests.so : ${TESTS_OBJS}
	${CLINKER} -shared ${TESTS_OBJS} -o $@

pismr : pismr.o libpism.so
	${CLINKER} $< ${ICE_LIB_FLAGS} -o $@

pisms : iceEISModel.o iceHEINOModel.o iceROSSModel.o pisms.o libpism.so
	${CLINKER} iceEISModel.o iceHEINOModel.o iceROSSModel.o pisms.o ${ICE_LIB_FLAGS} -o $@

pismv : iCMthermo.o iceCompModel.o iceExactStreamModel.o pismv.o libpism.so libtests.so
	${CLINKER} iCMthermo.o iceCompModel.o iceExactStreamModel.o pismv.o ${ICE_LIB_FLAGS} -o $@

pant : pant.o libpism.so
	${CLINKER} $< ${ICE_LIB_FLAGS} -o $@

pgrn : iceGRNModel.o pgrn.o libpism.so
	${CLINKER} iceGRNModel.o pgrn.o ${ICE_LIB_FLAGS} -o $@

#shelf : shelf.o libpism.so
#	${CLINKER} $< ${ICE_LIB_FLAGS} -o $@

flowTable : flowTable.o libpism.so
	${CLINKER} $< ${ICE_LIB_FLAGS} -o $@

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

simpleL : simpleL.o libtests.so
	${CLINKER} $< ${TESTS_LIB_FLAGS} -o $@

gridL : gridL.o libtests.so
	${CLINKER} $< ${TESTS_LIB_FLAGS} -o $@

# Cancel the implicit rules  WHY?
#% : %.cc
#% : %.c

.makeremind :
	@touch .makeremind
	@echo '*** Remember to "make install".  For now, new executables are "./pismv", etc. ****'

showEnv :
	@echo ${CLINKER}
	@echo ${PETSC_DIR}
	@echo ${PETSC_LIB}
	@echo ${ICE_LIB_FLAGS}
	@echo ${TESTS_LIB_FLAGS}

# Emacs style tags
.PHONY: tags TAGS
tags TAGS :
	etags *.cc *.hh *.c

# The GNU recommended procedure for automatically generating dependencies.
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
	rm -f TAGS \
	libpism.so libtests.so lib/libpism.so lib/libtests.so \
	${executables} $(patsubst %, bin/%, ${executables})\
	${extra_execs} $(patsubst %, bin/%, ${extra_execs})

.PHONY: clean

include $(depfiles)

