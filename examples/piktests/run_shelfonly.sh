#!/bin/bash

N=2
xx=151   # FIXME:  for quick regression, try xx=51,yy=51
yy=151
length=1 # FIXME:  for quick regression, make long enough to cause multiple time steps

infile="circular_shelfonly.nc"

grid="-Mx $xx -My $yy -Mz 31 -Mbz 5 -Lz 1500 -Lbz 1000"

pismopts="-boot_file $infile $grid -verbose 3 -ssa_sliding -ssa_dirichlet_bc "

#doit="mpiexec -n $N ../../bin/pismr $pismopts"
doit="mpiexec -n $N pismr $pismopts"

# run with strength extension, the old PISM method
#$doit $pismopts -y $length -o old_so.nc

# check that this result is similar
#$doit -y $length -cold -o old_cold_so.nc

# run with strength extension and part_grid but no CFBC
# this could be a regression for -part_grid and -part_redist
$doit -y $length -part_grid -part_redist -o part_so.nc

# run with CFBC but no part_grid
# this could be a regression for -ssa_method fd_pik only
$doit -y $length -ssa_method fd -cfbc -o cfbc_so.nc

# run with CFBC and part_grid
# this could be a regression for -ssa_method fd_pik only
$doit $pismopts -y $length -ssa_method fd -cfbc -part_grid -part_redist -o partcfbc_so.nc

