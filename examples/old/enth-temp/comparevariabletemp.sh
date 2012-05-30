#!/bin/bash

# go to examples/searise-greenland/ and preprocess first, then link to the
# relevant files:
#   $ ln -s ../searise-greenland/searise_config.nc
#   $ ln -s ../searise-greenland/pism_Greenland_5km_v1.1.nc

NN=2

COMMON="-config_override searise_config.nc -title 'SeaRISE Greenland Spinup' -skip 10 -boot_file pism_Greenland_5km_v1.1.nc -Mx 76 -My 141 -Lz 4000 -Lbz 2000 -Mz 101 -Mbz 11 -z_spacing equal -atmosphere searise_greenland -surface pdd -ocean_kill -y 100"

# constant conductivity and constant specific heat:
mpiexec -n $NN pismr           $COMMON -o g20km_pre100.nc

# variable conductivity and constant specific heat:
mpiexec -n $NN varkpismr       $COMMON -o vark_g20km_pre100.nc

# constant conductivity and variable specific heat:
mpiexec -n $NN pismr -varc     $COMMON -o varc_g20km_pre100.nc

# variable conductivity and variable specific heat:
mpiexec -n $NN varkpismr -varc $COMMON -o both_g20km_pre100.nc

