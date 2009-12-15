# Configuration flags:
PISM_PREFIX ?= $(PWD)
PISM_HAVE_FFTW ?= 1
PISM_STATIC ?= 0
# if ==1 then adds -g
PISM_USE_DEBUG ?= 0
# PETSc has troubles choosing a linker which can link C++. PISM is C++. Setting
# this to zero would let PETSc choose a linker.
PISM_USE_MPICXX ?= 1
# if ==1 then adds -Woverloaded-virtual -pipe to CFLAGS
PISM_USE_GNU_FLAGS ?= 0

# Put additional make include files here: 
#CONFIG = config/macosx_macports

# Miscellaneous variables:
BUILD_DIR = $(PWD)/build
GOALS = $(MAKECMDGOALS)

ALL: all

update: svn_update
	@$(MAKE) rebuild

rebuild:
ifeq ($(shell touch src/revision; svnversion src/ | diff src/revision -),)
	@echo "src/ directory is up to date."
else
	@echo "Rebuilding PISM..."
	@$(MAKE) all
endif

#FIXME: this has undesirable effect that "make clean && make"  does not rebuild executables
#pismr pismv pismd pgrn pcctest flowTable tryLCbd gridL simple%:
#	$(MAKE) -C build ../bin/$@

depclean:
	@rm -f $(BUILD_DIR)/*.d

svn_update:
	@echo "Running 'svn update' ($(shell svn info |grep 'Repository Root'))..."
	@svn update

clean:
	@$(MAKE) -C $(BUILD_DIR) clean
	@$(MAKE) -C doc clean

userman refman browser summary fullbib installation:
	@cd doc && $(MAKE) $@

.DEFAULT:
	@cd $(BUILD_DIR) && $(MAKE) $@

.EXPORT_ALL_VARIABLES: ;
