#!/bin/bash

N=4
M=201
xx=$M
yy=$M
length=400

infile="circular_noshelf.nc"

./circular_dirichlet.py -Mx $M -My $M -o $infile

grid="-Mx $xx -My $yy -Mz 31 -Mbz 1 -Lz 1500 -Lbz 1000"

stressbalance="-ssa_method fd -stress_balance ssa -ssa_dirichlet_bc"

output_basename="test_eigencalving"

calving="-calving eigen_calving -eigen_calving_K 1e15"

diagnostics="thk,mask,velbar_mag,ice_area_specific_volume,velbar,tendency_of_ice_mass_due_to_discharge"

viewers="-view $diagnostics"

extra="-extra_times 10 -extra_vars $diagnostics -extra_file ${output_basename}_ex.nc"

ts="-ts_times 10 -ts_file ${output_basename}_ts.nc"

misc_options="-cfbc -part_grid -o_order zyx -energy none"

pismopts="-i $infile -bootstrap $grid $stressbalance $calving $viewers $extra $ts $misc_options"

doit="mpiexec -n $N pismr $pismopts"

set -x
# run with CFBC and part_grid
$doit $pismopts -y $length -o ${output_basename}_o.nc
