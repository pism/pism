#!/bin/bash

N=4
xx=151
yy=151
length=10

infile="circular_shelfonly.nc"

if [[ ! -r $infile ]]
then
    echo "generating the input file..."
    ./circular_dirichlet.py -o $infile -shelf
fi

grid="-Mx $xx -My $yy -Mz 31 -Mbz 5 -Lz 1500 -Lbz 1000"

extra="-extra_times 1 -extra_vars thk,mask,cbar,Href,velbar -extra_file "

pismopts="-boot_file $infile $grid -verbose 3 -ssa_sliding -ssa_dirichlet_bc"

doit="mpiexec -n $N pismr $pismopts"

# run with strength extension, the old PISM method
#$doit $pismopts -y $length -o so_old.nc

# check that this result is similar
#$doit -y $length -cold -o so_old_cold.nc


# run with strength extension and part_grid but no CFBC
# this could be a regression for -part_grid and -part_redist
$doit -y $length -part_grid -part_redist -o so_part.nc $extra so_ex_part.nc

# run with CFBC but no part_grid
# this could be a regression for -ssa_method fd_pik only
$doit -y $length -ssa_method fd -cfbc -o so_cfbc.nc $extra so_ex_cfbc.nc

# run with CFBC and part_grid
# this could be a regression for -ssa_method fd_pik only
$doit $pismopts -y $length -ssa_method fd -cfbc -part_grid -part_redist -o so_partcfbc.nc $extra so_ex_partcfbc.nc

