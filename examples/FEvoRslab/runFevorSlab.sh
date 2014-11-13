#!/bin/bash

set -x

gridMx=10
gridMz=7

nP=$(( (${gridMx}-1) * (${gridMz}-1) ))

echo "The grid is ${gridMx} by ${gridMz} (X by Z)"
echo "There are ${nP} FEvoR particles"

./FEvoRslab.py 

flowline.py -o pism-in.nc --expand -d y fevor-slab-in.nc

${1%/}/pismr -surface given -boot_file pism-in.nc -periodicity xy \
 -Mx ${gridMx} -My 6 -Lx 50.05 -Mz ${gridMz} -Lz 550 -y 1000 \
 -bed_smoother_range 0 \
 -stress_balance sia_fevor \
 -sia_fevor_use_constant_slope -sia_fevor_surface_slope_degrees -3 \
 -fevor_n_particles ${nP} \
 -extra_file ex.nc -extra_times 10 \
 -extra_vars distributions,recrystallizations,thk,flux_mag,velsurf_mag,h_x,h_y,taud_mag,enhancement_factor,tauxz,tauyz,pressure \
 -o pism-out.nc

flowline.py -o fevor-slab-out.nc --collapse -d y pism-out.nc
