#!/bin/bash

set -u
set -x
set -e

pism -eisII A \
     -Lz 5000 \
     -Mx 30 \
     -My 30 \
     -Mz 10 \
     -o input-1.nc \
     -y 1000 \
     ""

pism \
  -i input-1.nc \
  -max_dt 50 \
  -o input-2.nc \
  -output.scalar.buffer_size 100 \
  -output.scalar.filename ts.nc \
  -output.scalar.times 1 \
  -output.scalar.variables ice_volume,tendency_of_ice_mass \
  -y 1000 \
  ""

pism \
  -i input-2.nc \
  -max_dt 50 \
  -o result.nc \
  -output.scalar.append \
  -output.scalar.buffer_size 100 \
  -output.scalar.filename ts.nc \
  -output.scalar.times 1 \
  -output.scalar.variables ice_volume,tendency_of_ice_mass \
  -y 1000 \
  ""
