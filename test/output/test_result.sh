#!/bin/bash

set -u
set -x
set -e

eis2_grid="-Mx 3 -My 3 -Mz 3 -Lz 5000"

pism -eisII A -y 1000 ${eis2_grid} -o result.nc
