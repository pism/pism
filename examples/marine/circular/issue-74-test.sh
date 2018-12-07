#!/bin/bash

# This script uses the circular very similar to the one in run_noshelf.sh to confirm the
# issue #74 (https://github.com/pism/pism/issues/74), i.e. that when -calving float_kill
# is used the sub-shelf ice flux reported by PISM (the variable basal_mass_flux_floating)
# can be non-zero even when the total area of the floating ice (variable
# ice_area_glacierized_floating) is zero.

# This is due to the fact that during time-stepping the sub-shelf ice flux is computed
# before calving is applied. This error has an O(dt) character. This will be fixed once
# the 2D mass transport code gets an overhaul, but for now basal_mass_flux_floating should
# be attributed to calving. (This is consistent with what the model does right now.)

set -x
set -e

N=4
xx=101
yy=$xx
length=50

infile="circular_noshelf.nc"

./circular_dirichlet.py -Mx $xx -My $yy -o $infile

grid="-Mx $xx -My $yy -Mz 31 -Mbz 1 -Lz 1500 -Lbz 0"

pismopts="-i $infile -bootstrap $grid -stress_balance ssa -ssa_dirichlet_bc -o_order zyx -energy none -ssa_method fd -cfbc -part_grid"

doit="mpiexec -n $N pismr $pismopts"

extra="-extra_times 1 -extra_vars thk,mask,velbar_mag,ice_area_specific_volume,velbar,usurf,mass_fluxes -extra_file issue-74_ex.nc"
ts="-ts_file issue-74_ts.nc -ts_times 1"

$doit -y $length -o issue-74_o.nc $extra $ts -calving float_kill
