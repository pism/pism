#!/bin/bash

# do a basic Jakoshavn spinup; run as
#   ./spinup NN Mx My >> out.spin &

# the grid DEFAULTS below apply to:
#   $ ncdump -h jako.nc |head -n 4
#   netcdf jako {
#   dimensions:
#       y = 425 ;
#       x = 620 ;

NN=4  # default number of processors
if [ $# -gt 0 ] ; then  # if user says "spinup.sh 8" then NN=8
  NN="$1"
fi

Mx=208
if [ $# -gt 1 ] ; then  # if user says "spinup.sh 8 534" then NN=8 and Mx=620
  Mx="$2"
fi

My=143
if [ $# -gt 2 ] ; then  # if user says "spinup.sh 8 534 335" then NN=8 and Mx=620 and My=425
  My="$2"
fi

CLIMATEFILE=g5km_climate.nc
BCFILE=g5km_bc.nc

CLIMATE="-surface given,forcing -surface_given_file $CLIMATEFILE -force_to_thk jako.nc"
SKIP=5

cmd="mpiexec -n $NN pismo -boot_file jako.nc -no_model_strip 10 \
  -Mx $Mx -My $My -Lz 4000 -Lbz 1000 -Mz 201 -Mbz 51 -z_spacing equal \
  -ocean_kill $CLIMATE -regrid_file $BCFILE -regrid_vars bmelt,bwat,enthalpy \
  -y 1 -skip -skip_max $SKIP -o jako3km_y1.nc"
echo "running:   $cmd"
$cmd
echo

#exit  # FIXME: temporary

FINALSPIN=2000
FS_EXDT=20    # about 100 frames

# OLD:
#  -topg_to_phi 5.0,20.0,-300.0,700.0 -diffuse_bwat -thk_eff \
#  -ssa_sliding -plastic_pwfrac 0.95 -pseudo_plastic_q 0.25 \

# useful diagnostic:  -ssa_view_nuh

cmd="mpiexec -n $NN pismo -i jako3km_y1.nc -no_model_strip 10 \
  -ocean_kill -cfbc -kill_icebergs \
  -topg_to_phi 5.0,30.0,-300.0,700.0 -diffuse_bwat -thk_eff \
  -ssa_sliding -plastic_pwfrac 0.95 -pseudo_plastic_q 0.15 \
  -extra_file ex_jako3km_0.nc -extra_times -$FINALSPIN:$FS_EXDT:0 \
  -extra_vars thk,cbase,bwp,tauc,dhdt,hardav,csurf,temppabase,diffusivity,bmelt,tempicethk_basal \
  -ts_file ts_jako3km_0.nc -ts_times -$FINALSPIN:yearly:0 \
  -ssa_dirichlet_bc -regrid_file $BCFILE -regrid_vars vel_ssa_bc \
  $CLIMATE -ys -$FINALSPIN -ye 0 -skip -skip_max $SKIP -o jako3km_0.nc"
echo "running:   $cmd"
$cmd

# good postprocess:
#ncap -O -s "dtau=taud_mag-tauc" jako3km_0.nc jako3km_0.nc 

