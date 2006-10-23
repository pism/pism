SHELL = /bin/sh
VPATH = src
ALL : all

#FLAGS:

WITH_NETCDF?=1
CFLAGS+= -DWITH_NETCDF=${WITH_NETCDF} -pipe

WITH_FFTW?=1
CFLAGS+= -DWITH_FFTW=${WITH_FFTW} -pipe

WITH_GSL?=1
CFLAGS+= -DWITH_GSL=${WITH_GSL} -pipe

#VARIABLES:

petsc_executables= flowTable pismr pismv pisms shelf get_drag
executables= $(petsc_executables) simpleFG simpleBCD

ice_sources= extrasGSL.cc grid.cc iMbasal.cc iMbeddef.cc iMdefaults.cc\
	iMgrainsize.cc iMIO.cc iMIOnetcdf.cc iMmacayeal.cc iMoptions.cc\
	iMtemp.cc iMutil.cc iMvelocity.cc iMviewers.cc iceCompModel.cc\
	iceModel.cc materials.cc 
csources= exactTestsBCD.c exactTestsFG.c cubature.c simpleBCD.c simpleFG.c
exec_sources= flowTable.cc simplify.cc run.cc verify.cc get_drag.cc shelf.cc

depfiles= $(ice_sources:.cc=.d) $(csources:.c=.d) $(exec_sources:.cc=.d)

ICE_OBJS= $(patsubst %.cc, %.o, ${ice_sources}) cubature.o

include ${PETSC_DIR}/bmake/common/base

ICE_LIB_FLAGS= -L`pwd`/obj -Wl,-rpath,`pwd`/obj -lpism -ltests ${PETSC_LIB}

ifneq (${WITH_NETCDF}, 0)
	ICE_LIB_FLAGS+= -lnetcdf_c++ -lnetcdf
endif
ifneq (${WITH_FFTW}, 0)
	ICE_LIB_FLAGS+= -lfftw3
endif
ifneq (${WITH_GSL}, 0)
	ICE_LIB_FLAGS+= -lgsl -lgslcblas
endif

#general TARGETS:

all : depend libpism libtests $(executables)

libpism : ${ICE_OBJS}
	${CLINKER} -shared -o obj/libpism.so ${ICE_OBJS} ${PETSC_LIB}

libtests : exactTestsBCD.o exactTestsFG.o
	${CLINKER} -shared -o obj/libtests.so $^

#libpism : ${ICE_OBJS}
#	ar r libpism.a ${ICE_OBJS}


#executable TARGETS:

flowTable : obj/libpism.so flowTable.o
	${CLINKER} $^ ${ICE_LIB_FLAGS} -o obj/flowTable

get_drag : obj/libpism.so get_drag.o
	${CLINKER} $^ ${ICE_LIB_FLAGS} -o obj/get_drag

pismr : obj/libpism.so run.o
	${CLINKER} $^ ${ICE_LIB_FLAGS} -o obj/pismr

pisms : obj/libpism.so simplify.o
	${CLINKER} $^ ${ICE_LIB_FLAGS} -o obj/pisms

pismv : obj/libpism.so obj/libtests.so iceCompModel.o verify.o
	${CLINKER} $^ ${ICE_LIB_FLAGS} -o obj/pismv

shelf : obj/libpism.so shelf.o
	${CLINKER} $^ ${ICE_LIB_FLAGS} -o obj/shelf

simpleBCD : obj/libpism.so simpleBCD.o
	${CLINKER} $^ -lm -L`pwd`/obj -ltests -o obj/simpleBCD

simpleFG : obj/libpism.so simpleFG.o
	${CLINKER} $^ -lm -L`pwd`/obj -ltests -o obj/simpleFG

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
	rm -f obj/libpism.so obj/libpism.a obj/libtests.so\
	      $(patsubst %, obj/%, ${executables}) *.d TAGS

.PHONY: clean

include $(depfiles)
