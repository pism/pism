#!/bin/bash

# compare IceModel and IceEnthalpyModel result for short runs from green20km_y1.nc

NN=8  # number of processors
set -e  # exit on error

#mpiexec -n $NN pgrn -boot_from eis_green_smoothed.nc \
#  -Mx 83 -My 141 -Lz 4000 -Mz 51 -Lbz 2000 -skip 1 -y 1 -o green20km_y1.nc

mpiexec -n 8 pismr -e 3 -ocean_kill -surface pdd -atmosphere searise_greenland -skip 10 -y 500 -i green20km_y1.nc -o pismr_y501.nc

mpiexec -n 8 penth -e 3 -ocean_kill -init_from_temp -surface pdd -atmosphere searise_greenland -skip 10 -y 500 -i green20km_y1.nc -o penth_y501.nc

