#!/bin/bash

# This script uses the circular very similar to the one in
# run_noshelf.sh to confirm the issue #74
# (https://github.com/pism/pism/issues/74), i.e. that when -float_kill
# is used the sub-shelf ice flux reported by PISM (the
# sub_shelf_ice_flux variable) can be non-zero even when the total
# area of the floating ice (variable iareaf) is zero.

# This is due to the fact that during time-stepping the sub-shelf ice
# flux is computed before calving is applied. This error has an O(dt)
# character. This will be fixed once the 2D mass transport code gets
# an overhaul, but for now sub_shelf_ice_flux should be attributed to
# calving. (This is consistent with what the model does right now.)

N=4
xx=101
yy=$xx
length=400

infile="circular_noshelf.nc"

if [[ ! -r $infile ]]
then
    echo "generating the input file..."
    ./circular_dirichlet.py -o $infile
fi

grid="-Mx $xx -My $yy -Mz 31 -Mbz 1 -Lz 1500 -Lbz 0"

pismopts="-boot_file $infile $grid -ssa_sliding -ssa_dirichlet_bc -o_order zyx -energy none -no_sia"

doit="mpiexec -n $N pismr $pismopts"

extra="-extra_times 10 -extra_vars thk,mask,cbar,Href,velbar,usurf -extra_file issue-74_ex.nc"
ts="-ts_file issue-74_ts.nc -ts_times 1"

$doit $pismopts -y $length -ssa_method fd -cfbc -part_grid -part_redist -o issue-74_o.nc $extra $ts -float_kill
