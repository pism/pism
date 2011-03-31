#!/bin/bash

N=2
xx=151   # FIXME:  for quick regression, try xx=51,yy=51
yy=151
length=1 # FIXME:  for quick regression, make long enough to cause multiple time steps
infile=circular_with_shelf_12km.nc

grid="-Mx $xx -My $yy -Mz 31 -Mbz 5 -Lz 4500 -Lbz 1000"

pismopts="-boot_file $infile $grid -surface constant -ssa_sliding"

doit="mpiexec -n $N pismr $pismopts"

# run circular_with_shelf.py first to generate input file circular_with_shelf_12km.nc

# run with strength extension, the old PISM method
$doit -y $length -o old.nc

# check that this result is similar
$doit -y $length -cold -o old_cold.nc

# run with strength extension and part_grid but no CFBC
# this could be a regression for -part_grid and -part_redist
$doit -y $length -part_grid -part_redist -o part.nc

# run with CFBC but no part_grid
# this could be a regression for -ssa_method fd_pik only
$doit -y $length -ssa_method fd_pik -o cfbc.nc

# run with CFBC and part_grid
# this could be a regression for -ssa_method fd_pik only
# FIXME: with N=4 processors and xx=yy=301 I observe slight asymmetry in cbase here?
$doit -y $length -part_grid -part_redist -ssa_method fd_pik -o partcfbc.nc  

# note FIXME: I think "-ssa_method fd_pik" in above should be "-ssa_method fd -cfbc"
