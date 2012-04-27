#!/bin/bash

#ncrename -v ice_surface_temp,artm -O Ross_combined.nc
#ncatted -a standard_name,climatic_mass_balance,a,c,land_ice_surface_specific_mass_balance -O Ross_combined.nc
#ncatted -a standard_name,thk,a,c,land_ice_thickness -O Ross_combined.nc
#ncatted -a standard_name,topg,a,c,bedrock_altitude -O Ross_combined.nc

mpiexec -n 1 pismr -boot_file Ross_combined.nc -Mx 211 -My 211 -Mz 21 -Lz 3000 -ssa_floating_only -pik -ssa_dirichlet_bc -y 0 -o out.nc -o_order zyx -surface given
