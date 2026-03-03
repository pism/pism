#!/bin/bash

set -u
set -x
set -e

eis2_grid="-Mx 20 -My 20 -Mz 3 -Lz 5000"

spatial="-o_size none -spatial_file ${filename:-ex.nc} -spatial_times 2000 -spatial_split ${split:-no}"

# spatial file
pism -eisII A -y 1e4 ${eis2_grid} ${spatial} -spatial_vars thk
