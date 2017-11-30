#!/bin/bash

N=4
xx=151
yy=151
length=10

infile=circular_withshelf.nc

./circular_ice_sheet.py -Mx $xx -My $yy -o $infile

grid="-Mx $xx -My $yy -Mz 31 -Mbz 1 -Lz 4500 -Lbz 1000"

pismopts="-y $length -i $infile -bootstrap $grid -stress_balance ssa+sia -ssa_method fd"

doit="mpiexec -n $N pismr $pismopts"

extra="-extra_times 1 -extra_vars thk,mask,velbar_mag,ice_area_specific_volume,velbar -extra_file"

# run with strength extension and part_grid but no CFBC
# this could be a regression for the option combination "-part_grid"
$doit -part_grid -o ws_part.nc $extra ws_ex_part.nc

# run with CFBC but no part_grid
# this could be a regression for the option "-cfbc"
$doit -cfbc -o ws_cfbc.nc $extra ws_ex_cfbc.nc

# run with CFBC and part_grid
# this could be a regression for the option combination "-part_grid -cfbc"
# FIXME: with N=4 processors and xx=yy=301 I observe slight asymmetry in velbase_mag here?
$doit -part_grid -cfbc -o ws_partcfbc.nc $extra ws_ex_partcfbc.nc
