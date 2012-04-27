#!/bin/bash

mpiexec -n 1 pismr -boot_file Ross_combined.nc -Mx 211 -My 211 -Mz 21 -Lz 3000 -ssa_floating_only -pik -ssa_dirichlet_bc -y 0 -o out.nc -o_order zyx