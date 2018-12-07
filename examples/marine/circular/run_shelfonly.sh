#!/bin/bash

N=4
xx=151
yy=151
length=10

infile="circular_shelfonly.nc"

./circular_dirichlet.py -Mx $xx -My $yy -o $infile -shelf

grid="-Mx $xx -My $yy -Mz 31 -Mbz 1 -Lz 1500 -Lbz 1000"

extra="-extra_times 1 -extra_vars thk,mask,velbar_mag,ice_area_specific_volume,velbar -extra_file "

pismopts="-i $infile -bootstrap $grid -stress_balance ssa+sia -ssa_dirichlet_bc"

doit="mpiexec -n $N pismr $pismopts"

# run with strength extension, the old PISM method
#$doit $pismopts -y $length -o so_old.nc

# check that this result is similar
#$doit -y $length -energy cold -o so_old_cold.nc


# run with strength extension and part_grid but no CFBC
# this could be a regression for -part_grid
$doit -y $length -part_grid -o so_part.nc $extra so_ex_part.nc

# run with CFBC but no part_grid
# this could be a regression for -ssa_method fd_pik only
$doit -y $length -ssa_method fd -cfbc -o so_cfbc.nc $extra so_ex_cfbc.nc

# run with CFBC and part_grid
# this could be a regression for -ssa_method fd_pik only
$doit $pismopts -y $length -ssa_method fd -cfbc -part_grid -o so_partcfbc.nc $extra so_ex_partcfbc.nc
