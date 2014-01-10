#!/bin/bash

# -no_sia takes care of high velocities at inland boundaries (SIA "sees"
#     discontinuities of surface and bedrock elevation fields).
# -pik sets the calving front boundary condition
#     (can be reduced to '-cfbc -kill_icebergs')
# -ssa_method fem will not work as '-pik', including '-cfbc', are not implemented in fem

NN=1  # default number of processors
if [ $# -gt 0 ] ; then  # if user says "run_prog.sh 8" then NN = 8
  NN="$1"
fi

M=211   # for 5km
# M=526 for 2km,  M=421 for 2.5km,  M=351 for 3km
if [ $# -gt 1 ] ; then  # if user says "run_prog.sh 8 211" then NN = 8 and M = 211
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
interval=25

ECALV=1e17   #  constant for eigencalving parameterization
if [ $# -gt 4 ] ; then  # if user says "run_prog.sh 8 211 0.6 500 7e16" then ... and -eigen_calving_K 7e16
  ECALV="$5"
fi

THRESHOLD=7.0e4   #  stress threshold
if [ $# -gt 5 ] ; then  # if user says "run_prog.sh 8 211 0.6 500 7e16 7.0e4" then ... and -fractures x,7e16,x,x
  THRESHOLD="$6"
fi

FRACRATE=0.5   #  fracture rate
if [ $# -gt 6 ] ; then  # if user says "run_prog.sh 8 211 0.6 500 7e16 7.0e4 0.5 " then ... and -fractures 0.5,x,x,x
  FRACRATE="$7"
fi

HEALTHRESHOLD=5.0e-10   #  healing threshold
if [ $# -gt 7 ] ; then  # if user says "run_prog.sh 8 211 0.6 500 7e16 0.5 5.0e-10" then ... and -fractures x,x,x,5.0e-10
  HEALTHRESHOLD="$8"
fi

HEALRATE=0.1   #  healing threshold
if [ $# -gt 8 ] ; then  # if user says "run_prog.sh 8 211 0.6 500 7e16 0.5 5.0e-10 0.1" then ... and -fractures x,x,0.1,x
  HEALRATE="$9"
fi

SOFTRES=1.0   #  healing threshold
if [ $# -gt 9 ] ; then  # if user says "run_prog.sh 8 211 0.6 500 7e16 0.5 5.0e-10 0.1 0.001" then ... and -fracture_softening 0.001
  SOFTRES="$10"
fi

# options ###############################

PISMPREFIX=""
#PISMPREFIX="../../../bin/"

output="-o Mx${M}_year-000${YEARS}.nc -o_order zyx -o_size big"

ssa="-no_sia -ssa_floating_only -ssa_dirichlet_bc -ssa_e ${SSAE} -part_grid -cfbc "
#-pik:-part_grid -cfbc -kill_icebergs -part_redist

#calving="-eigen_calving -eigen_calving_K ${ECALV} -thickness_calving -calving_at_thickness 50.0 "
calving="-ocean_kill "

extra="-extra_file ex_Mx${M}.nc -extra_times ${interval}:${interval}:${YEARS} -extra_vars thk,mask,csurf,IcebergMask,fracture_density,fracture_flow_enhancement,fracture_growth_rate,fracture_healing_rate,fracture_toughness"

timeseries="-ts_file ts-prog.nc -ts_times 0:1:${YEARS}"

# fractures ##############################################################################

criterion=""
#criterion="-lefm"
#criterion="-max_shear"

boundary="-do_frac_on_grounded"
#boundary="-phi0"

healing=""
#healing="-constant_healing" #independent of strain rates
#healing="-fracture_weighted_healing"

#softening="-fracture_softening 1.0" #no softening
softening="-fracture_softening ${SOFTRES}" #residual eps=0.001

fractures="-fractures ${FRACRATE},${THRESHOLD},${HEALRATE},${HEALTHRESHOLD} -write_fd_fields -scheme_fd2d ${healing} ${boundary} ${criterion} ${softening}"

#-constant_fd

# run commands #############################################################################

cmd_diag="mpiexec -n $NN ${PISMPREFIX}pismr -boot_file Ross_combined_prog.nc -Mx $M -My $M \
  -Mz 61 -Lz 3000 -z_spacing equal -surface given ${ssa} -kill_icebergs \
  -y 0 -ys 0.0 -o startfile_Mx${M}.nc -o_order zyx -fractures 0,0,0,0 -write_fd_fields "

cmd_prog_frac="mpiexec -n $NN ${PISMPREFIX}pismr -i startfile_Mx${M}.nc -surface given \
  ${ssa} -y ${YEARS} ${output} ${calving} ${fractures} ${extra} ${timeseries} -verbose 4 "



# -ssa_rtol 1.0e-3 -ssa_eps 5.0e15
# -cfl_eigencalving
# -nuBedrock 1e15

echo "running command:"
echo
echo "$cmd_diag"
echo
${cmd_diag}
echo
echo "$cmd_prog_frac"
echo
${cmd_prog_frac}
echo
echo "... done"
