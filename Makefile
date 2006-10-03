SHELL = /bin/sh
VPATH = src
ALL : all

WITH_NETCDF?=1
CFLAGS+= -DWITH_NETCDF=${WITH_NETCDF} -pipe

WITH_FFTW?=1
CFLAGS+= -DWITH_FFTW=${WITH_FFTW} -pipe

WITH_GSL?=1
CFLAGS+= -DWITH_GSL=${WITH_GSL} -pipe


petsc_executables= flowTable run_ice simplify verify shelf get_drag
executables= $(petsc_executables) simpleFG
sources= exactTestsBCD.cc exactTestsFG.cc flowTable.cc grid.cc  \
	iMbasal.cc iMbeddef.cc iMdefaults.cc iMgrainsize.cc iMIO.cc iMIOnetcdf.cc iMmacayeal.cc  \
	iMoptions.cc iMtemp.cc iMutil.cc iMvelocity.cc iMviewers.cc  \
	iceCompModel.cc iceModel.cc materials.cc extrasGSL.cc \
	simplify.cc run_ice.cc simpleFG.cc verify.cc get_drag.cc shelf.cc
csources= cubature.c

depfiles= $(sources:.cc=.d) $(csources:.c=.d)
ICE_OBJS= grid.o materials.o iceModel.o iMbasal.o iMbeddef.o iMgrainsize.o iMIO.o iMIOnetcdf.o \
	iMmacayeal.o iMtemp.o iMutil.o iMvelocity.o iMviewers.o iMdefaults.o \
	iMoptions.o extrasGSL.o cubature.o

include ${PETSC_DIR}/bmake/common/base

LIBICE= libice.a
ICE_LIB= -L`pwd` -Wl,-rpath,`pwd` -lice ${PETSC_LIB}
ifneq (${WITH_NETCDF}, 0)
	ICE_LIB+= -lnetcdf_c++ -lnetcdf
endif
ifneq (${WITH_FFTW}, 0)
	ICE_LIB+= -lfftw3
endif
ifneq (${WITH_GSL}, 0)
	ICE_LIB+= -lgsl -lgslcblas
endif

all : depend $(executables)
$(petsc_executables) : ${LIBICE}

libice.so : ${ICE_OBJS}
	${CLINKER} -shared -o libice.so ${ICE_OBJS} ${PETSC_LIB}

libice.a : ${ICE_OBJS}
	ar cr libice.a ${ICE_OBJS}

# If we are using GNU Make, this is sufficient
simpleFG : simpleFG.o exactTestsFG.o
	${LINK.cc} $^ $(LOADLIBES) $(LDLIBS) -o $@

verify : exactTestsBCD.o exactTestsFG.o iceCompModel.o verify.o
	${CLINKER} $^ ${ICE_LIB} -o $@

# The general rule for Petsc executables
% : %.o
	${CLINKER} $^ ${ICE_LIB} -o $@

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

distclean : clean
	rm -f libice.so libice.a $(executables) *.d TAGS

.PHONY: clean

include $(depfiles)
