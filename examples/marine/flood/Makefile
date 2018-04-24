input=-i flood.nc -bootstrap -Mx 201 -My 101 -Mz 21 -Lz 1500
ssa=-stress_balance ssa+sia -plastic_phi 15
ocean=-sea_level constant,delta_sl -ocean_delta_sl_file flood.nc
misc=-energy none
extra=-extra_vars thk,mask,velbar_mag -extra_times 10

options=${input} ${ssa} ${ocean} ${misc} ${extra}

all: no_cfbc cfbc cfbc_kill_icebergs

# this works for the wrong reason
no_cfbc: flood.nc
	pismr ${options} -o $@.nc -extra_file ex_$@.nc

# this fails (but should fail more gracefully)
cfbc: flood.nc
	pismr ${options} -cfbc -o $@.nc -extra_file ex_$@.nc

# this works the way it should
cfbc_kill_icebergs: flood.nc
	mpiexec -n 4 pismr ${options} -cfbc -kill_icebergs -o $@.nc -extra_file ex_$@.nc

flood.nc: generate_input.py
	python generate_input.py -M 1001 -o $@

clean:
	rm -f *.nc *.nc~
