#!/bin/bash

set -u
set -x
set -e

pism -eisII A -Mx 20 -My 20 -Mz 11 -y 1000 -o input.nc

# The state of the Lingle-Clark bed deformation model contains a variable
# (viscous_bed_displacement) that uses a different grid, so this output file will contain
# *two* grids.
pism -i input.nc -bootstrap -bed_def lc -y 100 -o result.nc
