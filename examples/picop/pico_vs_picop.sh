#!/bin/bash

set -x
set -u
set -e

NN=8
if [ $# -gt 0 ] ; then
  NN="$1" # first arg is number of processes
  echo "Using $NN cores"
fi

run_length=1s
spatial_vars="beta,bmelt,mask,topg,usurf,thk,velsurf_mag,velbase_mag,climatic_mass_balance,taub_mag"

# prepare the input file so we can test PICOP without running the SSA solver every time:
# mpiexec -n $NN pism \
# 	-bootstrap \
# 	-i ../antarctica/pism_Antarctica_5km.nc \
# 	-ys 0 \
# 	-y $run_length \
# 	-basal_resistance.pseudo_plastic.enabled yes \
# 	-basal_yield_stress.mohr_coulomb.till_phi_default 45 \
# 	-Lz 4500 \
# 	-Mz 5 \
# 	-y 1s \
# 	-atmosphere uniform \
# 	-surface simple \
# 	-atmosphere.uniform.precipitation 0 \
# 	-stress_balance ssa+sia \
# 	-cfbc \
# 	-ssafd_ksp_monitor \
# 	-o input.nc

# ncrename -v u_ssa,ubar -v v_ssa,vbar -O input.nc input.nc

# Run PISM with PICO:
out=pico_${run_length}.nc
mpiexec -n $NN pism \
	-i input.nc \
	-ys 0 \
	-y $run_length \
	-atmosphere uniform \
	-surface simple \
	-ocean pico \
	-ocean.pico.file ../../test/regression/pico_split/bedmap2_schmidtko14_50km.nc \
	-atmosphere.uniform.precipitation 0 \
	-o_size medium \
	-output.sizes.medium uvel,vvel,sftgif,velsurf_mag,mask,usurf,bmelt,velbar \
	-output.extra.file spatial_$out \
	-output.extra.times 1s \
	-output.extra.vars pico_temperature,pico_salinity,pico_basal_melt_rate,$spatial_vars \
	-output.file $out \
	-stress_balance prescribed_sliding+sia \
	-stress_balance.prescribed_sliding.file input.nc

# Run PISM with PICOP:
out=picop_${run_length}.nc
mpiexec -n $NN pism \
	-i input.nc \
	-ys 0 \
	-y $run_length \
	-atmosphere uniform \
	-surface simple \
	-ocean picop \
	-ocean.pico.file ../../test/regression/pico_split/bedmap2_schmidtko14_50km.nc \
	-atmosphere.uniform.precipitation 0 \
	-o_size medium \
	-output.sizes.medium uvel,vvel,sftgif,velsurf_mag,mask,usurf,bmelt,velbar \
	-output.extra.file spatial_$out \
	-output.extra.times 1s \
	-output.extra.vars picop_temperature,picop_salinity,picop_basal_melt_rate,picop_grounding_line_elevation,picop_grounding_line_slope,picop_shelf_base_elevation,$spatial_vars \
	-output.file $out \
	-stress_balance prescribed_sliding+sia \
	-stress_balance.prescribed_sliding.file input.nc
