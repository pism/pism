input=-boot_file flood.nc -Mx 101 -My 101 -Mz 21 -Lz 1500
ssa=-ssa_sliding -plastic_phi 15
ocean=-ocean constant,delta_SL -ocean_delta_SL_file flood.nc
misc=-no_energy
extra=-extra_vars thk,mask,cbar -extra_times 10

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
	pismr ${options} -cfbc -kill_icebergs -o $@.nc -extra_file ex_$@.nc

flood.nc: generate_input.py
	python generate_input.py -o $@

clean:
	rm -f *.nc *.nc~
