#!/bin/bash

# compare polythermal and cold ice methods, using short runs based
# on a recomputed "big" version of green20km_y1.nc which has variable 'temp'

set -e  # exit on error
set -x  # show commands as they go

NN=2  # number of processors

LENGTH=500
OPTIONS="-e 3 -ocean_kill -surface pdd -atmosphere searise_greenland -skip 10"

IN=green20km_y1_big.nc

mpiexec -n $NN pgrn -boot_file eis_green_smoothed.nc \
  -Mx 83 -My 141 -Lz 4000 -Mz 51 -Lbz 2000 -skip 1 -y 1 -o $IN -o_size big

# polythermal version
mpiexec -n $NN pismr $OPTIONS -y $LENGTH -i $IN -o poly_y$LENGTH.nc

# cold ice version
mpiexec -n $NN pismr -cold $OPTIONS -y $LENGTH -i $IN -o cold_y$LENGTH.nc

