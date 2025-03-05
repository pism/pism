#!/bin/bash

set -x
set -e
set -u

run_length=300000
output_file=out.nc
grid="-Mx 3 -My 3 -Mz 100 -Lz 1000 -Mbz 1 -Lbz 0 -z_spacing equal"

pism \
  $grid \
  -atmosphere uniform,delta_T \
  -atmosphere.delta_T.file slab_dT.nc \
  -atmosphere.uniform.precipitation 0 \
  -atmosphere.uniform.temperature 243.15 \
  -bootstrap \
  -bootstrapping.defaults.geothermal_flux 0.042 \
  -extra_file ex.nc \
  -extra_times 1000 \
  -extra_vars bmelt,tillwat,temppabase,tempbase,tempsurf \
  -hydrology.tillwat_max 1000 \
  -i slab.nc \
  -max_dt 10 \
  -no_mass \
  -o $output_file \
  -surface simple \
  -y $run_length \
  ;
