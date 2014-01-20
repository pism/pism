#!/bin/bash

# -no_sia takes care of high velocities at inland boundaries (SIA "sees"
#     discontinuities of surface and bedrock elevation fields).
# -pik sets the calving front boundary condition
#     (can be reduced to '-cfbc -kill_icebergs')
# -ssa_method fem seems not to work?

NN=2  # default number of processors
if [ $# -gt 0 ] ; then  # if user says "run.sh 8" then NN = 8
  NN="$1"
fi

M=211   # for 5km
# M=526 for 2km,  M=421 for 2.5km,  M=351 for 3km
if [ $# -gt 1 ] ; then  # if user says "run.sh 8 351" then NN = 8 and M = 351
  M="$2"
fi

SSAE=0.6
if [ $# -gt 2 ] ; then  # if user says "run.sh 8 351 0.7" then ... and -ssa_e 0.7
  SSAE="$3"
fi

PISMPREFIX=""
#PISMPREFIX="../../../bin/"

cmd="mpiexec -n $NN ${PISMPREFIX}pismr -boot_file ../Ross_combined_prog.nc -Mx $M -My $M \
  -Mz 3 -Lz 3000 -z_spacing equal -surface given -no_sia -no_energy \
  -ssa_floating_only -pik -ssa_dirichlet_bc -ssa_view_nuh \
  -y 0 -o Ross_result_diag_Mx${M}.nc -o_order zyx -ssa_e $SSAE -ssafd_ksp_monitor"

echo "running command:"
echo
echo "$cmd"
echo
$cmd
echo
echo "... done"
