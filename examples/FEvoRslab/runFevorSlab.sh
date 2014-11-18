#!/bin/bash

if (( $# != 2 )); then
  echo "Usage: $0 PATH_TO_FEvoRslab.py PATH_TO_pismr"
  exit 1
fi
 
set -x

gridMx=11
gridMz=11

nP=$(( (${gridMx}-1) * (${gridMz}-1) ))

echo "The grid is ${gridMx} by ${gridMz} (X by Z)"
echo "There are ${nP} FEvoR particles"

${1%/}/FEvoRslab.py 

flowline.py -o pism-in.nc --expand -d y fevor-slab-in.nc

${2%/}/pismr -surface given -boot_file pism-in.nc -periodicity xy \
 -Mx ${gridMx} -My ${gridMx} -Lx 50.05 -Ly 50.5 -Mz ${gridMz} -Lz 1000 -y 150 \
 -bed_smoother_range 0 -z_spacing equal \
 -stress_balance sia_fevor \
 -energy none \
 -sia_fevor_use_constant_slope -sia_fevor_bed_slope_degrees -3 \
 -extra_file ex.nc -extra_times 10 -extra_vars distributions,enhancement_factor,thk,velsurf,flux_divergence,diffusivity_staggered,taud_mag,enhancement_factor,h_x,enthalpysurf,ice_surface_temp,climatic_mass_balance,enthalpy \
 -fevor_n_particles ${nP} \
 -o pism-out.nc -o_order zyx

flowline.py -o fevor-slab-out.nc --collapse -d y pism-out.nc

exit 0
