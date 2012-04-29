#!/bin/bash

# -no_sia takes care of high velocities at inland boundaries (SIA "sees"
#     discontinuities of surface and bedrock elevation fields).
# -pik sets the calving front boundary condition
#     (and can be replaced by '-cfbc -part_grid -kill_icebergs' or '-cfbc -kill_icebergs' but
#     NOT(?) '-cfbc' alone or '-cfbc -part_grid')
# -ssa_method fem seems not to work?

N=211   # for 5km
# N=526 for 2km,  N=421 for 2.5km,  N=351 for 3km

mpiexec -n 4 pismr -boot_file Ross_combined.nc -Mx $N -My $N \
  -Mz 21 -Lz 3000 -z_spacing equal -surface given -no_sia -no_energy \
  -ssa_floating_only -pik -ssa_dirichlet_bc -ssa_view_nuh \
  -y 0 -o out_$N.nc -o_order zyx
