#!/bin/bash

N=4
xx=151
yy=151
length=1

infile="test_iceberg_removal.nc"

./test_iceberg_removal.py -Mx $xx -My $yy -o $infile -shelf -square

grid="-Mx $xx -My $yy -Mz 31 -Mbz 5 -Lz 1500 -Lbz 1000"

pc_asm=" -ssafd_ksp_type gmres -ssafd_ksp_norm_type unpreconditioned -ssafd_ksp_pc_side right -ssafd_pc_type asm -ssafd_sub_pc_type lu "

pismopts="-i $infile -bootstrap $grid -stress_balance ssa+sia -ssa_dirichlet_bc -ssa_view_nuh -view mask,thk -calving eigen_calving,thickness_calving -eigen_calving_K 1e18 -thickness_calving_threshold 300 -o_order zyx -ssafd_ksp_max_it 75 $pc_asm"


doit="mpiexec -n $N pismr $pismopts"

extra="-extra_times 0.05 -extra_vars thk,mask,velbar_mag,ice_area_specific_volume,velbar,mass_fluxes -extra_file iceberg_ex.nc"

# run with CFBC and part_grid
$doit -y $length -ssa_method fd -cfbc -part_grid -o iceberg_o.nc $extra
