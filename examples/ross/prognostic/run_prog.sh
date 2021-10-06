#!/bin/bash

NN=2  # default number of processors
if [ $# -gt 0 ] ; then  # if user says "run_prog.sh 8" then NN = 8
  NN="$1"
fi

M=211   # for 5km;  note M=526 for 2km,  M=421 for 2.5km,  M=351 for 3km
if [ $# -gt 1 ] ; then  # if user says "run_prog.sh 8 211" then NN = 8 and M = 351
  M="$2"
fi

SSAE=0.6   # ssa enhancement factor
if [ $# -gt 2 ] ; then  # if user says "run_prog.sh 8 211 0.7" then ... and -ssa_e 0.7
  SSAE="$3"
fi

YEARS=3000    # integration time
if [ $# -gt 3 ] ; then  # if user says "run_prog.sh 8 211 0.7 200" then ... and -y 200
  YEARS="$4"
fi

PISMPREFIX=""
#PISMPREFIX="../../../bin/"

# use these for more robust SSA solves
# STRONGKSP="-ssafd_ksp_type gmres -ssafd_ksp_norm_type unpreconditioned -ssafd_ksp_pc_side right -ssafd_pc_type asm -ssafd_sub_pc_type lu"
# STRONGKSP=""

# preliminary bootstrap and diagnostic run:
STARTNAME=startfile_Mx${M}.nc
cmd_diag="mpiexec -n $NN ${PISMPREFIX}pismr -regional -i ../Ross_combined.nc -bootstrap -Mx $M -My $M \
  -Mz 61 -Lz 3000 -z_spacing equal -surface given -stress_balance ssa \
  -yield_stress constant -tauc 1e6 -pik -ssa_dirichlet_bc -ssa_e $SSAE \
  $STRONGKSP -y 0.1 -ys 0.0 -o $STARTNAME -o_order zyx "
echo "running command:"
echo
echo "$cmd_diag"
${cmd_diag}

# prognostic run
NAME=prog_Mx${M}_yr${YEARS}.nc
ECALV=1e18   #  constant for eigen_calving parameterization
CTHICK=50.0  #  constant thickness for thickness_calving
exdt=1
exvars="thk,mask,velsurf_mag,strain_rates,tendency_of_ice_mass_due_to_flow,tendency_of_ice_mass_due_to_discharge"
cmd_prog="mpiexec -n $NN ${PISMPREFIX}pismr -regional -i $STARTNAME \
  -surface given -stress_balance ssa -yield_stress constant -tauc 1e6 -pik \
  -calving eigen_calving,thickness_calving -eigen_calving_K $ECALV -front_retreat_cfl \
  -ssa_dirichlet_bc -ssa_e $SSAE -ys 0 -y $YEARS -o $NAME -o_size big \
  -thickness_calving_threshold $CTHICK $STRONGKSP \
  -ts_file ts-${NAME} -ts_times 0:monthly:${YEARS} \
  -extra_file ex-${NAME} -extra_times 0:${exdt}:${YEARS} -extra_vars ${exvars} \
  -options_left"
echo "running command:"
echo
echo "$cmd_prog"
${cmd_prog}

