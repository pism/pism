common_options=-bootstrap -Mx 360 -My 561 -y 1day -atmosphere uniform -surface simple  -regrid_vars topg,thk

std_greenland=~/github/pism/pism/examples/std-greenland/pism_Greenland_5km_v1.1.nc
result_std.nc: grid.nc
	mpiexec -n 1 pismr -o $@ -i $^ ${common_options} -regrid_file ${std_greenland}

highres=~/github/ckhroulev/pism-flexible-interpolation/data/gris_900m.nc
result_highres.nc: grid.nc
	mpiexec -n 8 pismr -o $@ -i $^ ${common_options} -regrid_file ${highres}

grid.nc: grid_definition.cdl
	ncgen -o $@ $^
