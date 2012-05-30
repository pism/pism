#!/bin/bash

N=2
xx=151   # FIXME:  for quick regression, try xx=51,yy=51
yy=151
length=1 # FIXME:  for quick regression, make sure the run is long enough to cause multiple time steps
infile=circular_withshelf_12km.nc

grid="-Mx $xx -My $yy -Mz 31 -Mbz 5 -Lz 4500 -Lbz 1000"

pismopts="-boot_file $infile $grid -ssa_sliding -ssa_method fd"

#doit="mpiexec -n $N ../../bin/pismr $pismopts"
doit="mpiexec -n $N pismr $pismopts"

# run circular_with_shelf.py first to generate input file circular_with_shelf_12km.nc

# run with strength extension and part_grid but no CFBC
# this could be a regression for the option combination "-part_grid -part_redist"
$doit -y $length -part_grid -part_redist -o part_ws.nc

# run with CFBC but no part_grid
# this could be a regression for the option "-cfbc"
$doit -y $length -cfbc -o cfbc_ws.nc

# run with CFBC and part_grid
# this could be a regression for the option combination "-part_grid -part_redist -cfbc"
# FIXME: with N=4 processors and xx=yy=301 I observe slight asymmetry in cbase here?
$doit -y $length -part_grid -part_redist -cfbc -o partcfbc_ws.nc  

