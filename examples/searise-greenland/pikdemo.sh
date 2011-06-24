#!/bin/bash

# Copyright (C) 2011 The PISM Authors

# a tiny bit of the SeaRISE-Greenland spinup in a version that demonstrates
#   PISM-PIK handling of ice shelves

# on 10km grid; about 3.5 Gb memory used

set -e  # exit on error

NN=2  # default number of processors
if [ $# -gt 0 ] ; then  # if user says "spinup.sh 8" then NN = 8
  NN="$1"
fi

#(spinup.sh)  bootstrapping plus short smoothing run (for 100a)
mpiexec -n $NN pismr -config_override searise_config.nc -skip 20 \
  -boot_file pism_Greenland_5km_v1.1.nc \
  -Mx 151 -My 281 -Lz 4000 -Lbz 2000 -Mz 201 -Mbz 21 -z_spacing equal \
  -atmosphere searise_greenland -surface pdd \
  -y 100 -o pik_g10km_pre100.nc \
  -pik -eigen_calving 2.0e18 -calving_at_thickness 50.0

# now turn on ssa
mpiexec -n $NN pismr -config_override searise_config.nc -skip 20 \
  -i pik_g10km_pre100.nc \
  -ssa_sliding -thk_eff -topg_to_phi 5.0,20.0,-300.0,700.0,10.0 \
  -atmosphere searise_greenland -surface pdd \
  -y 50 -o pik_g10km_ssa50.nc \
  -pik -eigen_calving 2.0e18 -calving_at_thickness 50.0
