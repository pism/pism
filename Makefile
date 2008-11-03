# Configuration flags:
PISM_PREFIX ?= $(PWD)
PISM_HAVE_FFTW ?= 1
PISM_LOG_EVENTS ?= 0
# PETSc has troubles choosing a linker which can link C++. PISM is C++. Setting
# this to zero would let PETSc choose a linker.
USE_MPICXX = 1

# Put additional make include files here: 
#CONFIG = config/ryan_make
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
	@$(MAKE) all install
endif

install: local_install

depclean:
	@rm -f $(BUILD_DIR)/*.d

svn_update:
	@echo "Running 'svn update' ($(shell svn info |grep 'Repository Root'))..."
	@svn update

clean:
	@make -C $(BUILD_DIR) clean
	@make -C doc clean

userman refman browser summary fullbib installation:
	@cd doc && $(MAKE) $@

.DEFAULT:
	@cd $(BUILD_DIR) && $(MAKE) $@

.EXPORT_ALL_VARIABLES: ;

