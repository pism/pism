#!/bin/bash

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

grid="-Mx $xx -My $yy -Mz 31 -Mbz 5 -Lz 1500 -Lbz 1000"

pismopts="-boot_file $infile $grid -stress_balance ssa+sia -ssa_dirichlet_bc -o_order zyx"

doit="mpiexec -n $N pismr $pismopts"

extra="-extra_times 10 -extra_vars thk,mask,velbar_mag,Href,velbar,usurf -extra_file ns_ex.nc"
ts="-ts_file ns_ts.nc -ts_times 1"

$doit $pismopts -y $length -ssa_method fd -cfbc -part_grid -part_redist -o ns_o.nc $extra $ts
