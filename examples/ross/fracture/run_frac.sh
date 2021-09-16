#!/bin/bash

NN=2  # default number of processors
if [ $# -gt 0 ] ; then  # if user says "run_frac.sh 8" then NN = 8
  NN="$1"
fi

M=211   # for 5km;  note M=526 for 2km,  M=421 for 2.5km,  M=351 for 3km
if [ $# -gt 1 ] ; then  # if user says "run_frac.sh 8 211" then NN = 8 and M = 211
  M="$2"
fi

SSAE=0.6   # ssa enhancement factor
if [ $# -gt 2 ] ; then  # if user says "run_frac.sh 8 211 0.6" then ... and -ssa_e 0.6
  SSAE="$3"
fi

YEARS=3000    # integration time
if [ $# -gt 3 ] ; then  # if user says "run_frac.sh 8 211 0.6 500" then ... and -y 500
  YEARS="$4"
fi
exdt=25 # for the extrafile

FRACTHRESHOLD=1.3e5   #  stress threshold

FRACRATE=0.5   #  fracture rate

HEALTHRESHOLD=2.0e-10   #  healing threshold

HEALRATE=0.15   #  healing rate

SOFTRES=0.01   #  softening residual (avoid viscosity from degeneration), value 1 inhibits softening effect


# options ###############################

PISMPREFIX=""
#PISMPREFIX="../../../bin/"

NAME=frac_Mx${M}_yr-${YEARS}.nc

output="-o $NAME -o_order zyx -o_size big"

ssa="-stress_balance ssa -yield_stress constant -tauc 1e6 -ssa_dirichlet_bc -ssa_e ${SSAE} -part_grid -cfbc "
#-pik:-part_grid -cfbc -kill_icebergs

extra="-extra_file ex-${NAME} -extra_times 0:${exdt}:${YEARS} \
       -extra_vars thk,mask,velsurf_mag,fracture_density,fracture_flow_enhancement,fracture_growth_rate,fracture_healing_rate,fracture_toughness"

timeseries="-ts_file ts-${NAME} -ts_times 0:1:${YEARS}"

criterion=""

boundary="-do_frac_on_grounded"

healing=""
#healing="-constant_healing" #independent of strain rates
#healing="-fracture_weighted_healing"

#softening="-fracture_softening 1.0" #no softening
softening="-fracture_softening ${SOFTRES}" #residual eps=0.001

fractures="-fractures -fracture_parameters  ${FRACRATE},${FRACTHRESHOLD},${HEALRATE},${HEALTHRESHOLD} \
           -scheme_fd2d ${healing} ${boundary} ${criterion} ${softening}"


# run commands #############################################################################

cmd_diag="mpiexec -n $NN ${PISMPREFIX}pismr -i ../Ross_combined.nc -bootstrap -Mx $M -My $M \
  -Mz 61 -Lz 3000 -z_spacing equal -surface given ${ssa} -kill_icebergs \
  -y 0 -ys 0.0 -o startfile_Mx${M}.nc -o_order zyx -fractures -fracture_parameters 0,0,0,0 -write_fd_fields "

# add "-verbose 4" to this command for additional internal info
cmd_frac="mpiexec -n $NN ${PISMPREFIX}pismr -i startfile_Mx${M}.nc -surface given \
  ${ssa} -y ${YEARS} ${output} -front_retreat_file startfile_Mx${M}.nc \
  ${fractures} ${extra} ${timeseries}"

# -ssafd_picard_rtol 1.0e-3 -ssa_eps 5.0e15

echo "running command:"
echo
echo "$cmd_diag"
echo
${cmd_diag}

ncap2 -A \
      -s 'land_ice_area_fraction_retreat=0 * thk' \
      -s 'where(thk > 0 || topg > 0) land_ice_area_fraction_retreat=1' \
      -s 'land_ice_area_fraction_retreat@units="1"' \
      startfile_Mx${M}.nc startfile_Mx${M}.nc
ncatted -a standard_name,land_ice_area_fraction_retreat,d,, startfile_Mx${M}.nc

echo
echo "$cmd_frac"
echo
${cmd_frac}
echo
echo "... done"
