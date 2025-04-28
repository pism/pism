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
  -output.timeseries.buffer_size 100 \
  -output.timeseries.filename ts.nc \
  -output.timeseries.times 1 \
  -output.timeseries.variables ice_volume \
  -y 1000 \
  ""

pism \
  -i input-2.nc \
  -max_dt 50 \
  -o result.nc \
  -output.timeseries.append \
  -output.timeseries.buffer_size 100 \
  -output.timeseries.filename ts.nc \
  -output.timeseries.times 1 \
  -output.timeseries.variables ice_volume \
  -y 1000 \
  ""
