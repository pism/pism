#!/bin/bash

N=2

# run circular_with_shelf.py first

# FIXME: works for option "-cold", but not in default case
#        (see bug #17950 at gna.org)
mpiexec -n $N pismr -boot_file circular_with_shelf_12km.nc -surface constant \
  -Mx 151 -My 151 -Mz 31 -Mbz 5 -Lz 4500 -Lbz 2000 \
  -ssa_sliding -y 1 -cold -o result1cold.nc

mpiexec -n $N pismr -boot_file circular_with_shelf_12km.nc -surface constant \
  -Mx 151 -My 151 -Mz 31 -Mbz 5 -Lz 4500 -Lbz 2000 \
  -ssa_sliding -y 1 -o result1enth.nc

# this could be a regression for -part_grid and -part_redist
# just to better distinguish between icefree ocean and floating ice
mpiexec -n $N pismr -boot_file circular_with_shelf_12km.nc -surface constant \
  -Mx 151 -My 151 -Mz 31 -Mbz 5 -Lz 4500 -Lbz 2000 \
  -ssa_sliding -y 1 -part_grid -part_redist -o result1_part.nc

