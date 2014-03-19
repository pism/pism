M=101
grid=-Mx ${M} -My ${M} -Mz 3 -Lz 200
ssa=-stress_balance ssa+sia -ssa_dirichlet_bc -cfbc
mass=-part_grid -part_redist
misc=-o_order zyx -energy none -o output.nc
extra=-extra_times 5 -extra_file ex.nc -extra_vars cbar,thk,mask,taud,velbar
length=-y 1500

options=${grid} ${ssa} ${mass} ${misc} ${length}

default: tongues.nc
	mpiexec -n 4 pismr -boot_file tongues.nc ${options} ${extra}

extra_calving=-extra_times 5 -extra_file ex.nc -extra_vars taud,cbar,thk,mask,discharge_flux_cumulative,velbar
eigen_calving: tongues.nc
	mpiexec -n 4 pismr -boot_file tongues.nc -calving eigen_calving -eigen_calving_K 2e18 \
                           ${options} ${extra_calving}

eigen_plus_thickness_calving: tongues.nc
	mpiexec -n 4 pismr -boot_file tongues.nc -calving eigen_calving,thickness_calving \
	                   -eigen_calving_K 2e18 -thickness_calving_threshold 50 \
                           ${options} ${extra_calving}

tongues.nc: generate_input.py
	python ./generate_input.py -o $@ -M ${M}

clean:
	rm -f *.nc *~
