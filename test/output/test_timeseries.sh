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
   -output.timeseries.buffer_size 100 \
   -output.timeseries.filename ts.nc \
   -output.timeseries.times 1 \
   -output.timeseries.variables ice_volume \
   -y 1000 \
   ""
