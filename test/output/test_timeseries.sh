#!/bin/bash

set -u
set -x
set -e

pism -eisII A \
   -Lz 5000 \
   -Mx 30 \
   -My 30 \
   -Mz 10 \
   -o_size none \
   -output.scalar.buffer_size 100 \
   -output.scalar.filename ts.nc \
   -output.scalar.times 1 \
   -output.scalar.variables ice_volume,tendency_of_ice_mass \
   -y 1000 \
   ""
