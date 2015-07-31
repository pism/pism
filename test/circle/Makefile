


circle_test: 
	ncgen -o circle_config.nc circle_config.cdl
	mpirun -n 2 ~/pism/install/bin/pisms -eisII A -Mx 101 -My 101 -Mz 41 -Lz 5000 -y 1e-4 -o circle_template.nc  -extra_vars thk,velsurf_mag,h_x,h_y,rank -extra_times 8e-5  -o_size big -extra_file stress_file.nc
	ncap2 -vs  'climatic_mass_balance=climatic_mass_balance*0+100;thk=thk*0+1000;ice_surface_temp=ice_surface_temp*0+240;topg=topg*0' circle_template.nc input_circle.nc
	./generate_input
	mpirun -np 2 ~/pism/install/bin/pismr \
	-y 1e-6 -ys 0 \
	-i input_circle.nc  -bootstrap  -Mz 41 -Lz 1100  -config_override circle_config.nc -o circle_init.nc \
	-extra_file circle_ex -extra_times 100  -extra_vars thk,velsurf_mag,h_x,h_y,rank,diffusivity,uvelsurf,vvelsurf  -extra_split  \
	-gradient prescribed -gradient_file stress_file.nc -seed_times 10000 |tee  circle_init.log

	mpirun -np 2 ~/pism/install/bin/pismr \
	-y 1e4 \
	-i circle_init.nc  -config_override circle_config.nc -o circle.nc \
	-extra_file circle_ex -extra_times 100  -extra_vars thk,velsurf_mag,h_x,h_y,rank,diffusivity,uvelsurf,vvelsurf  -extra_split  \
	-gradient prescribed -gradient_file stress_file.nc -seed_times 10000 |tee circle.log

	select_tracers_by_region circle_ex-00000100.000.nc  -x 295000 300000 -y -30000 -28000 |head -n 1     | trace_tracer circle_ex-0000*.nc && plot_tracers tracers_traced.nc
	./compare_tracer.py

all: circle_test
