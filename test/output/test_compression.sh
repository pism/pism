#!/bin/bash

set -u
set -x
set -e

eis2_grid="-Mx 3 -My 3 -Mz 3 -Lz 5000"

# choosing compression level
pism \
   -eisII A \
   -o result_compressed.nc \
   -output.compression_level 4 \
   -output.format netcdf4_serial \
   -y 1000 ${eis2_grid} \
   ""
