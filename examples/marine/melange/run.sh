#!/bin/bash

xx=151
yy=$xx
length=400

input="-i input.nc -bootstrap"

grid="-Mx $xx -My $yy -Mz 11 -Mbz 1 -Lz 1500 -Lbz 0 -y $length"

physics="-stress_balance ssa+sia -ssa_dirichlet_bc -cfbc -part_grid"

spatial_output="-spatial_vars thk,mask,velbar_mag,velbar,ice_area_specific_volume"

spatial_output="$spatial_output -spatial_times 5 -spatial_file spatial.nc"
scalar_output="-scalar_file scalar.nc -scalar_times 1"

output="-o o.nc $spatial_output $scalar_output"

ocean="-ssa_method fd -ocean constant,delta_MBP -ocean.delta_MBP.file delta_MBP.nc"
mpiexec -n 4 pism $input $grid $physics $ocean $output

# Cut out a slice with x == 0 (through the center of the domain).
ncks -O -d x,$(( $xx / 2 )) spatial.nc center-mbp.nc

ocean="-ssa_method fd -ocean constant"
mpiexec -n 4 pism $input $grid $physics $ocean $output

# Cut out a slice with x == 0 (through the center of the domain).
ncks -O -d x,$(( $xx / 2 )) spatial.nc center-no-mbp.nc

ncdiff -O center-mbp.nc center-no-mbp.nc diff.nc
