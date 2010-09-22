#!/bin/bash

NN=2

mpiexec -n $NN pismr -ocean_kill -e 3 -skip 5 -boot_file pism_Greenland_5km_v0.93.nc -Mx 39 -My 71 -Lz 4000 -Lbz 2000 -Mz 21 -Mbz 6 -atmosphere searise_greenland -surface pdd -pdd_fausto -y 100 -o g40km_pre100.nc

mpiexec -n $NN pismr -ocean_kill -e 3 -skip 5 -i g40km_pre100.nc -atmosphere searise_greenland -surface pdd -pdd_fausto -no_mass -y 50000 -extra_file ex_g40km_steady.nc -extra_vars enthalpybase,temppabase,bmelt,bwat,csurf,hardav,mask -extra_times 0:500:50000 -o g40km_steady.n

mpiexec -n $NN pismr -ocean_kill -e 3 -skip 5 -i g40km_steady.nc -atmosphere searise_greenland -surface pdd -pdd_fausto -y 100 -o g40km_SIA.nc

mpiexec -n $NN pismr -ocean_kill -e 3 -skip 5 -i g40km_SIA.nc -topg_to_phi 5.0,20.0,-300.0,700.0,10.0 -ssa_sliding -thk_eff -pseudo_plastic_q 0.25 -plastic_pwfrac 0.98 -atmosphere searise_greenland,forcing -surface pdd -pdd_fausto -dTforcing pism_dT.nc -ocean constant -ys -125000 -y 20 -o g40km_m124980.nc

