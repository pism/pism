#!/bin/bash

N=4
M=101
xx=$M
yy=$M
length=400

infile="circular_noshelf.nc"

if [[ ! -r $infile ]]
then
    echo "generating the input file..."
    ./circular_dirichlet.py -o $infile
fi

grid="-Mx $xx -My $yy -Mz 31 -Mbz 1 -Lz 1500 -Lbz 1000"

stressbalance="-ssa_method fd -stress_balance ssa -ssa_dirichlet_bc"

output_basename="test_eigencalving"

calving="-calving eigen_calving -eigen_calving_K 1e15"

diagnostics="thk,mask,cbar,Href,velbar,discharge_flux_cumulative"

viewers="-view_map $diagnostics"

extra="-extra_times 10 -extra_vars $diagnostics -extra_file ${output_basename}_ex.nc"

ts="-ts_times 10 -ts_file ${output_basename}_ts.nc"

misc_options="-cfbc -part_grid -part_redist -o_order zyx"

pismopts="-boot_file $infile $grid $stressbalance $calving $viewers $extra $ts $misc_options"

doit="mpiexec -n $N pismr $pismopts"

# run with CFBC and part_grid
$doit $pismopts -y $length -o ${output_basename}_o.nc
