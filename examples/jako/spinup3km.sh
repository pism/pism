#!/bin/bash

NN=4

CLIMATE="-atmosphere given -atmosphere_given_file jako.nc -surface pdd,forcing -force_to_thk jako.nc"

NOMASSSIA=25000
FINALSPIN=2000

mpiexec -n $NN pismo -boot_file jako.nc -Mx 178 -My 112 -no_model_strip 10 \
      -Lz 4000 -Lbz 1000 -Mz 201 -Mbz 51 -z_spacing equal \
      $CLIMATE -sia_e 3.0 -y 1 -skip 30 -o jako3km_y1.nc

mpiexec -n $NN pismo -i jako3km_y1.nc -no_mass -sia_e 3.0 \
      $CLIMATE -y $NOMASSSIA -o jako3km_steadySIA.nc

mpiexec -n $NN pismo -i jako3km_steadySIA.nc -sia_e 1.5 -ocean_kill \
      -topg_to_phi 5.0,20.0,-300.0,700.0 \
      -ssa_sliding -plastic_pwfrac 0.95 -pseudo_plastic_q 0.25 \
      $CLIMATE -ys -$FINALSPIN -ye 0 -skip 30 -o jako3km_0.nc
