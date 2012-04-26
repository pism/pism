#!/bin/bash

NN=8

CLIMATE="-atmosphere given -atmosphere_given_file jako.nc -surface pdd,forcing -force_to_thk jako.nc"

NOMASSSIA=25000
NMS_EXDT=200  # about 100 frames
FINALSPIN=2000
FS_EXDT=20    # about 100 frames

mpiexec -n $NN pismo -boot_file jako.nc -Mx 178 -My 112 -no_model_strip 10 \
      -Lz 4000 -Lbz 1000 -Mz 201 -Mbz 51 -z_spacing equal \
      $CLIMATE -y 1 -skip 30 -o jako3km_y1.nc

mpiexec -n $NN pismo -i jako3km_y1.nc -no_mass \
      -extra_file ex_jako3km_steadySIA.nc -extra_times 0:$NMS_EXDT:$NOMASSSIA \
      -extra_vars hardav,csurf,temppabase,diffusivity,bmelt,tempicethk_basal \
      $CLIMATE -y $NOMASSSIA -o jako3km_steadySIA.nc

#FIXME: options here harvested from spinup.sh and ctrl_config.cdl in icerepo/UAF-regional/run/
#  (Dani's thesis has minimal comments on these)
mpiexec -n $NN pismo -i jako3km_steadySIA.nc -ocean_kill \
      -topg_to_phi 5.0,20.0,-300.0,700.0 -diffuse_bwat -thk_eff \
      -ssa_sliding -plastic_pwfrac 0.95 -pseudo_plastic_q 0.25 \
      -extra_file ex_jako3km_0.nc -extra_times -$FINALSPIN:$FS_EXDT:0 \
      -extra_vars thk,cbase,bwp,tauc,dhdt,hardav,csurf,temppabase,diffusivity,bmelt,tempicethk_basal \
      -ts_file ts_jako3km_0.nc -ts_times -$FINALSPIN:yearly:0 \
      $CLIMATE -ys -$FINALSPIN -ye 0 -skip 30 -o jako3km_0.nc
