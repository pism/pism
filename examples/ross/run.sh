#!/bin/bash

# -no_sia takes care of high velocities at inland boundaries (SIA "sees" discontinuities of surface and bedrock elevation fields).
# -ssa_sliding has a twofold effect: 1) essentially prescribing the 0 velocity inland and 2) (potentially) allowing some sliding over the Roosevelt island
# -pik sets the calving front boundary condition
mpiexec -n 4 pismr -boot_file Ross_combined.nc -Mx 211 -My 211 -Mz 21 -Lz 3000 \
                   -no_sia -ssa_sliding -cfbc -kill_icebergs -ssa_dirichlet_bc -y 0 -o out.nc -no_energy \
                   -o_order zyx -surface given
