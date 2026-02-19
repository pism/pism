#!/bin/bash

set -u
set -x
set -e

eis2_grid="-Mx 20 -My 20 -Mz 3 -Lz 5000"

extra="-o_size none -spatial_file ${filename:-ex.nc} -spatial_times 2000 -extra_split ${split:-no}"

# extra file
pism -eisII A -y 1e4 ${eis2_grid} ${extra} -spatial_vars thk
