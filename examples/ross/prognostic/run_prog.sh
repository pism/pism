#!/bin/bash

# -no_sia takes care of high velocities at inland boundaries (SIA "sees"
#     discontinuities of surface and bedrock elevation fields).
# -pik sets the calving front boundary condition
#     (can be reduced to '-cfbc -kill_icebergs')
# -ssa_method fem seems not to work?

NN=4  # default number of processors
if [ $# -gt 0 ] ; then  # if user says "run_prog.sh 8" then NN = 8
  NN="$1"
fi

M=211   # for 5km
# M=526 for 2km,  M=421 for 2.5km,  M=351 for 3km
if [ $# -gt 1 ] ; then  # if user says "run_prog.sh 8 211" then NN = 8 and M = 351
  M="$2"
fi

SSAE=0.6   # ssa enhancement factor
if [ $# -gt 2 ] ; then  # if user says "run_prog.sh 8 211 0.6" then ... and -ssa_e 0.6
  SSAE="$3"
fi

YEARS=500    # integration time
if [ $# -gt 3 ] ; then  # if user says "run_prog.sh 8 211 0.6 500" then ... and -y 500
  YEARS="$4"
fi
interval=2


ECALV=1e17   #  constant for eigencalving parameterization
if [ $# -gt 4 ] ; then  # if user says "run_prog.sh 8 211 0.6 500 7e16" then ... and -eigen_calving_K 7e16
  ECALV="$5"
fi

PISMPREFIX=""
# PISMPREFIX="../../../bin/"

cmd_diag="mpiexec -n $NN ${PISMPREFIX}pismr -boot_file Ross_combined_prog.nc -Mx $M -My $M \
  -Mz 61 -Lz 3000 -z_spacing equal -surface given -no_sia \
  -ssa_floating_only -pik -ssa_dirichlet_bc -ssa_e $SSAE \
  -y 0 -ys 0.0 -o startfile_Mx${M}.nc -o_order zyx "

cmd_prog="mpiexec -n $NN ${PISMPREFIX}pismr -i startfile_Mx${M}.nc \
  -surface given -no_sia -ssa_floating_only -pik -ssa_dirichlet_bc -ssa_e ${SSAE} \
  -y ${YEARS} -o Mx${M}_year-000${YEARS}.nc -o_order zyx -o_size big \
  -eigen_calving -eigen_calving_K ${ECALV} -thickness_calving -calving_at_thickness 50.0 \
  -ts_file ts-prog.nc -ts_times 0:1:${YEARS} \
  -extra_file ex_Mx${M}.nc -extra_times ${interval}:${interval}:${YEARS} -extra_vars thk,mask,csurf,IcebergMask"

# -ssa_rtol 1.0e-3 -ssa_eps 5.0e15
# -cfl_eigencalving
# -nuBedrock 1e15

echo "running command:"
echo
echo "$cmd_diag"
echo
${cmd_diag}
echo
echo "$cmd_prog"
echo
${cmd_prog}
echo
echo "... done"
