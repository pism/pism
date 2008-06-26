export PISM_PREFIX ?= $(PWD)
BUILD_DIR = $(PWD)/build

# Flags:
export WITH_FFTW ?= 1
export LOG_PISM_EVENTS ?= 0
export MISMIP_PLAY ?= 0

# PETSc has trouble choosing a linker which can link C++.  PISM is C++.
# Set this to zero to let PETSc choose a linker.
export USE_MPICXX = 1

# Put additional make include files here: 
#export CONFIG = config/ryan_make
#export CONFIG = config/macosx_macports

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

.DEFAULT:
	@cd $(BUILD_DIR) && $(MAKE) $@

