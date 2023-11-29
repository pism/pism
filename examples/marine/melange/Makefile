all: delta_MBP.nc input.nc
	./run.sh

delta_MBP.nc: generate_forcing.py
	PYTHONPATH=../../preprocessing:${PYTHONPATH} python3 $^ $@

input.nc:
	../circular/circular_dirichlet.py -Mx 151 -My 151 -o ./input.nc
