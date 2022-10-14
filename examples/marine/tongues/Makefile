M=101
grid=-Mx ${M} -My ${M} -Mz 3 -Lz 200
ssa=-stress_balance ssa+sia -ssa_dirichlet_bc -cfbc
mass=-part_grid
misc=-o_order zyx -energy none -o output.nc
extra=-spatial_times 5 -spatial_file ex.nc -spatial_vars velbar_mag,thk,mask,taud,velbar
length=-y 1500
ts=-scalar_file ts.nc -scalar_times 1

options=${grid} ${ssa} ${mass} ${misc} ${length} ${ts}

default: tongues.nc
	mpiexec -n 4 pismr -i tongues.nc -bootstrap ${options} ${extra}

extra_calving=-spatial_times 5 -spatial_file ex.nc -spatial_vars taud,velbar_mag,thk,mask,discharge_flux_cumulative,velbar,discharge_flux,vonmises_calving_rate,eigen_calving_rate
eigen_calving: tongues.nc
	mpiexec -n 4 pismr -i tongues.nc -bootstrap -calving eigen_calving -eigen_calving_K 2e18 \
                           ${options} ${extra_calving}

vonmises_calving: tongues.nc
	mpiexec -n 4 pismr -i tongues.nc -bootstrap -calving vonmises_calving -vonmises_calving_sigma_max 1e6\
                           ${options} ${extra_calving}

eigen_plus_thickness_calving: tongues.nc
	mpiexec -n 4 pismr -i tongues.nc -bootstrap -calving eigen_calving,thickness_calving \
	                   -eigen_calving_K 2e18 -thickness_calving_threshold 50 \
                           ${options} ${extra_calving}

tongues.nc: generate_input.py
	python ./generate_input.py -o $@ -M ${M}

clean:
	rm -f *.nc *~
