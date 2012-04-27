#!/bin/bash

# -no_sia takes care of high velocities at inland boundaries (SIA "sees"
#     discontinuities of surface and bedrock elevation fields).
# -ssa_sliding has a twofold effect: 1) essentially prescribing the zero
#     velocity inland and 2) (potentially) allowing some sliding over
#     Roosevelt island
# -pik sets the calving front boundary condition
#     (and can be replaced by '-cfbc -part_grid -kill_icebergs' but
#     NOT(?) '-cfbc' alone or '-cfbc -part_grid'j or '-cfbc -kill_icebergs')
# -ssa_method fem seems not to work?
mpiexec -n 4 pismr -boot_file Ross_combined.nc -Mx 211 -My 211 -Mz 21 -Lz 3000 \
  -surface given -no_sia -no_energy \
  -ssa_sliding -pik -ssa_dirichlet_bc \
  -y 0 -o out.nc -o_order zyx

