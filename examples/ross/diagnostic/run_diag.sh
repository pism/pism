#!/bin/bash

# '-pik' sets the calving front boundary condition
#     (is equivalent to '-cfbc -kill_icebergs')
# '-ssa_method fem' will not work because it lacks '-cfbc'

NN=2  # default number of processors
if [ $# -gt 0 ] ; then  # if user says "run_diag.sh 8" then NN = 8
  NN="$1"
fi

M=211   # for 5km;  note M=526 for 2km,  M=421 for 2.5km,  M=351 for 3km
if [ $# -gt 1 ] ; then  # if user says "run_diag.sh 8 351" then NN = 8 and M = 351
  M="$2"
fi

SSAE=0.6
if [ $# -gt 2 ] ; then  # if user says "run_diag.sh 8 351 0.7" then ... and '-ssa_e 0.7'
  SSAE="$3"
fi

PISMPREFIX=""
#PISMPREFIX="../../../bin/"

cmd="mpiexec -n $NN ${PISMPREFIX}pismr -i ../Ross_combined.nc -bootstrap -Mx $M -My $M \
  -Mz 3 -Lz 3000 -z_spacing equal -surface given -stress_balance ssa -energy none -no_mass \
  -yield_stress constant -tauc 1e6 -pik -ssa_dirichlet_bc \
  -y 1.0 -o diag_Mx${M}.nc -o_order zyx -ssa_e $SSAE -ssafd_ksp_monitor"

echo "running command:"
echo
echo "$cmd"
$cmd
