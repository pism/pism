#!/bin/bash

xx=151
yy=$xx
length=400

input="-i input.nc -bootstrap"

grid="-Mx $xx -My $yy -Mz 11 -Mbz 1 -Lz 1500 -Lbz 0 -y $length"

physics="-stress_balance ssa+sia -ssa_dirichlet_bc -cfbc -part_grid"

extra="-extra_vars thk,mask,velbar_mag,velbar,ice_area_specific_volume"

extra="$extra -extra_times 5 -extra_file ex.nc"
ts="-ts_file ts.nc -ts_times 1"

output="-o o.nc $extra $ts"

ocean="-ssa_method fd -ocean constant,delta_MBP -ocean.delta_MBP.file delta_MBP.nc"
mpiexec -n 4 pism $input $grid $physics $ocean $output

# Cut out a slice with x == 0 (through the center of the domain).
ncks -O -d x,$(( $xx / 2 )) ex.nc center-mbp.nc

ocean="-ssa_method fd -ocean constant"
mpiexec -n 4 pism $input $grid $physics $ocean $output

# Cut out a slice with x == 0 (through the center of the domain).
ncks -O -d x,$(( $xx / 2 )) ex.nc center-no-mbp.nc

ncdiff -O center-mbp.nc center-no-mbp.nc diff.nc
