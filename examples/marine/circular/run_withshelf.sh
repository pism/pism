#!/bin/bash

N=4
xx=151
yy=151
length=10

infile=circular_withshelf.nc

if [[ ! -r $infile ]]
then
    echo "generating the input file..."
    ./circular_ice_sheet.py -o $infile
fi

grid="-Mx $xx -My $yy -Mz 31 -Mbz 5 -Lz 4500 -Lbz 1000"

pismopts="-y $length -boot_file $infile $grid -stress_balance ssa+sia -ssa_method fd"

doit="mpiexec -n $N pismr $pismopts"

extra="-extra_times 1 -extra_vars thk,mask,cbar,Href,velbar -extra_file"

# run with strength extension and part_grid but no CFBC
# this could be a regression for the option combination "-part_grid -part_redist"
$doit -part_grid -part_redist -o ws_part.nc $extra ws_ex_part.nc

# run with CFBC but no part_grid
# this could be a regression for the option "-cfbc"
$doit -cfbc -o ws_cfbc.nc $extra ws_ex_cfbc.nc

# run with CFBC and part_grid
# this could be a regression for the option combination "-part_grid -part_redist -cfbc"
# FIXME: with N=4 processors and xx=yy=301 I observe slight asymmetry in cbase here?
$doit -part_grid -part_redist -cfbc -o ws_partcfbc.nc $extra ws_ex_partcfbc.nc

