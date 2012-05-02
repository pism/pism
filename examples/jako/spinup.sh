#!/bin/bash

NN=4  # default number of processors
if [ $# -gt 0 ] ; then  # if user says "spinup.sh 8" then NN=8
  NN="$1"
fi

Mx=178
if [ $# -gt 1 ] ; then  # if user says "spinup.sh 8 534" then NN=8 and Mx=534
  Mx="$2"
fi

My=112
if [ $# -gt 2 ] ; then  # if user says "spinup.sh 8 534 335" then NN=8 and Mx=534 and My=335
  My="$2"
fi

CLIMATEFILE=g5km_climate.nc
BCFILE=g5km_bc.nc

CLIMATE="-surface given,forcing -surface_given_file $CLIMATEFILE -force_to_thk jako.nc"

cmd="mpiexec -n $NN pismo -boot_file jako.nc -Mx $Mx -My $My -no_model_strip 10 -Lz 4000 -Lbz 1000 -Mz 201 -Mbz 51 -z_spacing equal -ocean_kill $CLIMATE -regrid_file $BCFILE -regrid_vars bmelt,enthalpy -y 1 -skip -skip_max 10 -o jako3km_y1.nc"
echo "running:   $cmd"
$cmd
echo

#exit  # FIXME: temporary

FINALSPIN=2000
FS_EXDT=20    # about 100 frames

cmd="mpiexec -n $NN pismo -i jako3km_y1.nc -ocean_kill \
  -topg_to_phi 5.0,20.0,-300.0,700.0 -diffuse_bwat -thk_eff \
  -ssa_sliding -plastic_pwfrac 0.95 -pseudo_plastic_q 0.25 \
  -extra_file ex_jako3km_0.nc -extra_times -$FINALSPIN:$FS_EXDT:0 \
  -extra_vars thk,cbase,bwp,tauc,dhdt,hardav,csurf,temppabase,diffusivity,bmelt,tempicethk_basal \
  -ts_file ts_jako3km_0.nc -ts_times -$FINALSPIN:yearly:0 \
  -ssa_dirichlet_bc -regrid_file $BCFILE -regrid_vars vel_ssa_bc \
  $CLIMATE -ys -$FINALSPIN -ye 0 -skip -skip_max 30 -o jako3km_0.nc"
echo "running:   $cmd"
$cmd

