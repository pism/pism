#!/bin/bash

N=2
xx=101   # FIXME:  for quick regression, try xx=51,yy=51
yy=101
length=400 

infile="circular_noshelf.nc"

grid="-Mx $xx -My $yy -Mz 31 -Mbz 5 -Lz 1500 -Lbz 1000"


pismopts="-boot_file $infile $grid -verbose 3 -ssa_sliding -ssa_dirichlet_bc"

#doit="mpiexec -n $N ../../bin/pismr $pismopts"
doit="mpiexec -n $N pismr $pismopts"

# run with strength extension, the old PISM method, FIXME: does not work yet with defined boundary conditions
#$doit $pismopts -y $length -o old_ns.nc

# run with CFBC and part_grid
$doit $pismopts -y $length -ssa_method fd -cfbc -part_grid -part_redist -o partcfbc_ns.nc
