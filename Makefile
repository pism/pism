PISM_INSTALL_PREFIX ?= $(PWD)
BUILD_DIR=build
PISM_STATIC ?= 0

ALL: install

install:
ifeq ($(PISM_STATIC),0)
	cd build && PISM_INSTALL_PREFIX=$(PISM_INSTALL_PREFIX) cmake ..
else
	cd build && PISM_INSTALL_PREFIX=$(PISM_INSTALL_PREFIX) cmake -DBUILD_SHARED_LIBS=OFF ..
endif
	$(MAKE) -C $(BUILD_DIR) install
.PHONY: install

update:
	@svn up
	$(MAKE) install
.PHONY: update

clean:
	@$(MAKE) -C $(BUILD_DIR) clean

.DEFAULT:
	$(MAKE) -C $(BUILD_DIR) $@
