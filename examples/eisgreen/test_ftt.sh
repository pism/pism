#!/bin/bash

NN=8

mpiexec -n $NN pismr -ys -500.0 -ye 0 -skip 5 -i green_ssl2_110ka.nc -atmosphere searise_greenland \
  -surface pdd,forcing -pdd_fausto -force_to_thk green20km_y1.nc \
  -o with_force.nc -ts_file ts_with_force.nc -ts_times -500:10:0
  
mpiexec -n $NN pismr -ys -500.0 -ye 0 -skip 5 -i green_ssl2_110ka.nc -atmosphere searise_greenland \
  -surface pdd -pdd_fausto \
  -o no_force.nc -ts_file ts_no_force.nc -ts_times -500:10:0

mpiexec -n $NN pismr -ys -500.0 -ye 0 -skip 5 -i green_ssl2_110ka.nc -atmosphere searise_greenland \
  -surface pdd,forcing -pdd_fausto -force_to_thk green20km_y1.nc -force_to_thk_alpha 0.0002 \
  -o weak_force.nc -ts_file ts_weak_force.nc -ts_times -500:10:0
  
