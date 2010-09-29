PISM_INSTALL_PREFIX ?= $(PWD)
BUILD_DIR=build

ALL: install

install:
	cd build && PISM_INSTALL_PREFIX=$(PISM_INSTALL_PREFIX) cmake ..
	$(MAKE) -C $(BUILD_DIR) install
.PHONY: install

userman browser installation:
	@cd doc && $(MAKE) $@

update:
	@svn up
	$(MAKE) install
.PHONY: update

clean:
	@$(MAKE) -C $(BUILD_DIR) clean
	@$(MAKE) -C doc clean

.DEFAULT:
	$(MAKE) -C $(BUILD_DIR) $@
