DX ?= 25
N ?= 4
DURATION ?= 300

all: figure.png

input.nc debris_input.nc: create_input.py
	python3 $< --dx $(DX) --debris-file debris_input.nc input.nc

debris.nc debris_spatial.nc debris_scalar.nc: input.nc debris_input.nc run.sh
	N=$(N) DX=$(DX) bash -x ./run.sh input.nc $(DURATION) debris

figure.png: plot_results.py debris_spatial.nc debris_scalar.nc
	python3 $< debris_spatial.nc debris_scalar.nc $@

.PHONY: clean
clean:
	rm -f *.nc *.log *.png *~
