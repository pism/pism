#!/bin/bash

set -u
set -x
set -e

# appending to an spatial file
spatial="-spatial_file ex.nc -spatial_times 100 -spatial_vars thk"

pism -eisII A -y 1 -Mx 3 -My 3 -Mz 3 -Lz 5000 -o input-1.nc

pism -i input-1.nc -y 1000 -o input-2.nc ${spatial}

pism -i input-2.nc -y 1000 ${spatial} -spatial_append
