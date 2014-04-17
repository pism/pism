#!/bin/bash

N=4
M=101
xx=$M
yy=$M
length=600

infile="circular_noshelf.nc"

output_basename="test_calving_at_thickness"

./circular_dirichlet.py -o $infile

grid="-Mx $xx -My $yy -Mz 3 -Mbz 1 -Lz 1500 -Lbz 0"

stressbalance="-ssa_method fd -stress_balance ssa -ssa_dirichlet_bc $pc"

calving="-calving thickness_calving -thickness_calving_threshold 300"

diagnostics="thk,mask,velbar_mag,Href,velbar,discharge_flux_cumulative"

viewers="-view_map $diagnostics"

extra="-extra_times 10 -extra_vars $diagnostics -extra_file ${output_basename}_ex.nc"

misc_options="-cfbc -part_grid -part_redist -o_order zyx -energy none"

pismopts="-boot_file $infile $grid $stressbalance $calving $viewers $extra $misc_options"

doit="mpiexec -n $N pismr $pismopts"

set -x

$doit $pismopts -y $length -o ${output_basename}_o.nc
