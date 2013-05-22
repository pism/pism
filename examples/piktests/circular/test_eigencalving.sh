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

grid="-Mx $xx -My $yy -Mz 31 -Mbz 5 -Lz 1500 -Lbz 1000"

pismopts="-boot_file $infile $grid -ssa_sliding -ssa_dirichlet_bc -eigen_calving -eigen_calving_K 1e14 -view_map mask,thk -o_order zyx -no_sia -ksp_type gmres -ksp_norm_type unpreconditioned -ksp_pc_side right -pc_type asm -sub_pc_type lu"

doit="mpiexec -n $N pismr $pismopts"

extra="-extra_times 10 -extra_vars thk,mask,cbar,Href,velbar -extra_file ns_ex.nc"

# run with CFBC and part_grid
$doit $pismopts -y $length -ssa_method fd -cfbc -part_grid -part_redist -o ns_partcfbc.nc $extra
