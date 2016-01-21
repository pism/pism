#!/bin/bash

set -x
set -e

# preprocess stage 1: create config file
cat <<EOF | ncgen -o conf.nc 
netcdf conf {
variables:
byte pism_overrides;
pism_overrides:hydrology_tillwat_max = 1000.0;
}
EOF

# preprocess stage 2: create bootstrap file
./generate_inputs.py

# run stage
pismexec="pismr"
ts=0
te=300000
oname=out_box.nc
grid="-Mx 3 -My 3 -Mz 100 -Lz 1000 -Mbz 1 -Lbz 0 -z_spacing equal -max_dt 10"
physics="-config_override conf.nc -no_mass"
surface="-surface given,delta_T -surface_delta_T_file box_dT.nc"
extra="-extra_file ex_box.nc -extra_times $ts:1000:$te -extra_vars bmelt,tillwat,temppabase,tempbase,tempsurf"

$pismexec -bootstrap -i box.nc $grid $physics $surface $extra -ys $ts -ye $te -o $oname
