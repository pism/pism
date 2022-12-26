all: out.nc

out.nc: run.sh slab.nc slab_dT.nc
	bash -x run.sh

slab.nc slab_dT.nc: generate_inputs.py
	python3 $^

.PHONY: clean
clean:
	@rm -f *.nc
