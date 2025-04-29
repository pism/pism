#!/bin/bash

set -u
set -x
set -e

eis2_grid="-Mx 3 -My 3 -Mz 3 -Lz 5000"

# set output.use_MKS just to exercise one more code branch
pism -eisII A -y 1000 ${eis2_grid} -output.use_MKS true -o result.nc
