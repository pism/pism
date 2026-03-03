#!/bin/bash

set -u
set -x
set -e

pism_dir=${pism_dir:-$HOME/github/pism/pism}
antarctica=${pism_dir}/test/regression/pico_split/bedmap2_schmidtko14_50km.nc

# The output should contain:
#
# - coordinates attribute when lon,lat are present
# - bounds of lon and lat when lon_bnds and lat_bnds re present
# - grid_mapping when PISM has projection information
pism \
  -Lz 4500 \
  -Mx 60 \
  -My 60 \
  -Mz 3 \
  -atmosphere uniform \
  -bootstrap \
  -spatial_file ex.nc \
  -spatial_times 100 \
  -spatial_vars thk,lon,lat,lon_bnds,lat_bnds \
  -i ${antarctica} \
  -o_size none \
  -surface simple \
  -y 200 \
   ""
