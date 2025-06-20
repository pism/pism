#!/bin/bash

set -u
set -x
set -e

# appending to an extra file
extra="-extra_file ex.nc -extra_times 100 -extra_vars thk"

pism -eisII A -y 1 -Mx 3 -My 3 -Mz 3 -Lz 5000 -o input-1.nc

pism -i input-1.nc -y 1000 -o input-2.nc ${extra}

pism -i input-2.nc -y 1000 ${extra} -extra_append
