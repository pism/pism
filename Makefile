# Configuration flags:
PISM_PREFIX ?= $(PWD)
WITH_FFTW ?= 1
LOG_PISM_EVENTS ?= 0
MISMIP_PLAY ?= 0
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
ifeq ($(shell touch src/revision; svnversion src/ | diff src/revision -),)
	@echo "src/ directory is up to date."
else
	@echo "Rebuilding PISM..."
	$(MAKE) all install
endif

install: local_install

depclean:
	@rm -f $(BUILD_DIR)/*.d

svn_update:
	@echo "Running 'svn update' ($(shell svn info |grep 'Repository Root'))..."
	@svn update

userman refman browser summary fullbib:
	@cd doc && $(MAKE) $@

.DEFAULT:
	@cd $(BUILD_DIR) && $(MAKE) $@

.EXPORT_ALL_VARIABLES: ;