#!/bin/bash

N=4
xx=101
yy=$xx
length=400 

infile="circular_noshelf.nc"

./circular_dirichlet.py -Mx $xx -My $yy -o $infile

grid="-Mx $xx -My $yy -Mz 31 -Mbz 5 -Lz 1500 -Lbz 1000"

pismopts="-i $infile -bootstrap $grid -stress_balance ssa+sia -ssa_dirichlet_bc -o_order zyx"

doit="mpiexec -n $N pism $pismopts"

spatial_output="-spatial_times 10 -spatial_vars thk,mask,velbar_mag,ice_area_specific_volume,velbar,usurf,mass_fluxes -spatial_file ns_spatial.nc"
scalar_output="-scalar_file ns_scalar.nc -scalar_times 1"

$doit $pismopts -y $length -ssa_method fd -cfbc -part_grid -o ns_o.nc $spatial_output $scalar_output
