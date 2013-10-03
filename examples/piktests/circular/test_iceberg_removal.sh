#!/bin/bash

N=4
xx=151
yy=151
length=1

infile="test_iceberg_removal.nc"

if [[ ! -r $infile ]]
then
    echo "generating the input file..."
    ./test_iceberg_removal.py -o $infile -shelf -square
fi

grid="-Mx $xx -My $yy -Mz 31 -Mbz 5 -Lz 1500 -Lbz 1000"

pc_asm=" -ksp_type gmres -ksp_norm_type unpreconditioned -ksp_pc_side right -pc_type asm -sub_pc_type lu "

pismopts="-boot_file $infile $grid -ssa_sliding -ssa_dirichlet_bc -ssa_view_nuh -kill_icebergs -view_map mask,thk -eigen_calving -eigen_calving_K 1e18 -thickness_calving -calving_at_thickness 300 -o_order zyx -ssafd_ksp_max_it 100 -ssafd_ksp_monitor"


doit="mpiexec -n $N pismr $pismopts"

extra="-extra_times 0.05 -extra_vars thk,mask,cbar,Href,velbar,IcebergMask -extra_file iceberg_ex.nc"

# run with CFBC and part_grid
$doit $pismopts -y $length -ssa_method fd -cfbc -part_grid -part_redist -o iceberg_o.nc $extra
