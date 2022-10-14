Mx=201
My=101

input=-i flood.nc -bootstrap -Mx ${Mx} -My ${My} -Mz 21 -Lz 1500
ssa=-stress_balance ssa+sia -plastic_phi 15
ocean=-sea_level constant,delta_sl -ocean_delta_sl_file flood.nc
misc=-energy none
extra=-spatial_vars thk,mask,velbar_mag -spatial_times 10

options=${input} ${ssa} ${ocean} ${misc} ${extra}

all: no_cfbc cfbc cfbc_kill_icebergs

# this works for the wrong reason
no_cfbc: flood.nc
	pismr ${options} -o $@.nc -spatial_file ex_$@.nc

# this fails
cfbc: flood.nc
	pismr ${options} -cfbc -o $@.nc -spatial_file ex_$@.nc || exit 0

# this works the way it should
cfbc_kill_icebergs: flood.nc
	mpiexec -n 4 pismr ${options} -cfbc -kill_icebergs -o $@.nc -spatial_file ex_$@.nc

flood.nc: generate_input.py
	PYTHONPATH=../../preprocessing:${PYTHONPATH} python3 $^ -M 1001 -o $@

clean:
	rm -f *.nc *.nc~
