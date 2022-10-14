#!/bin/bash

xx=151
yy=$xx
length=400

input="-i input.nc -bootstrap"

grid="-Mx $xx -My $yy -Mz 11 -Mbz 1 -Lz 1500 -Lbz 0 -y $length"

physics="-stress_balance ssa+sia -ssa_dirichlet_bc -cfbc -part_grid"

extra="-spatial_vars thk,mask,velbar_mag,Href,velbar"

extra="$extra -spatial_times 5 -spatial_file ex.nc"
ts="-ts_file ts.nc -ts_times 1"

output="-o o.nc -o_order zyx $extra $ts"

ocean="-ssa_method fd -ocean constant,delta_MBP -ocean_delta_MBP_file delta_MBP.nc"
mpiexec -n 4 pismr $input $grid $physics $ocean $output

# Cut out a slice with x == 0 (through the center of the domain).
ncks -O -d x,$(( $xx / 2 )) ex.nc center-mbp.nc

ocean="-ssa_method fd -ocean constant"
mpiexec -n 4 pismr $input $grid $physics $ocean $output

# Cut out a slice with x == 0 (through the center of the domain).
ncks -O -d x,$(( $xx / 2 )) ex.nc center-no-mbp.nc

ncdiff -O center-mbp.nc center-no-mbp.nc diff.nc
