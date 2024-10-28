M=101
grid=-Mx ${M} -My ${M} -Mz 3 -Lz 200
ssa=-stress_balance ssa+sia -ssa_dirichlet_bc -cfbc
mass=-part_grid
misc=-energy none -o output.nc
extra=-extra_times 5 -extra_file ex.nc -extra_vars velbar_mag,thk,mask,taud,velbar
length=-y 1500
ts=-ts_file ts.nc -ts_times 1

options=${grid} ${ssa} ${mass} ${misc} ${length} ${ts}

default: tongues.nc
	mpiexec -n 4 pism -i tongues.nc -bootstrap ${options} ${extra}

extra_calving=-extra_times 5 -extra_file ex.nc -extra_vars taud,velbar_mag,thk,mask,mass_fluxes,velbar
eigen_calving: tongues.nc
	mpiexec -n 4 pism -i tongues.nc -bootstrap -calving eigen_calving -eigen_calving_K 2e18 \
                           ${options} ${extra_calving},eigen_calving_rate

vonmises_calving: tongues.nc
	mpiexec -n 4 pism -i tongues.nc -bootstrap -calving vonmises_calving -vonmises_calving_sigma_max 1e6\
                           ${options} ${extra_calving},vonmises_calving_rate

eigen_plus_thickness_calving: tongues.nc
	mpiexec -n 4 pism -i tongues.nc -bootstrap -calving eigen_calving,thickness_calving \
	                   -eigen_calving_K 2e18 -thickness_calving_threshold 50 \
                           ${options} ${extra_calving},eigen_calving_rate

tongues.nc: generate_input.py
	PYTHONPATH=../../preprocessing:${PYTHONPATH} python3 ./generate_input.py -o $@ -M ${M}

clean:
	rm -f *.nc *~
