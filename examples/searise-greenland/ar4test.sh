#!/bin/bash

# THIS SCIPT IS TEMPORARY IN pism-dev

# use links to get anomalies from icerepo/SeaRISE-Greenland/data/future_forcing/

## compare control run

# mpiexec -n 8 pismr -ocean_kill -skip 5 -ssa_sliding -thk_eff -e 3 -pseudo_plastic_q 0.25 -plastic_pwfrac 0.98 -bed_def lc -i g20km_0.nc -atmosphere searise_greenland -surface pdd -pdd_fausto -ocean constant -ye 500 -extra_file ex_y500.nc -extra_times 0:5:500 -extra_vars csurf,cbase,usurf,topg,thk,bmelt,bwat,bwp,dHdt,mask,uvelsurf,vvelsurf,wvelsurf,uvelbase,vvelbase,wvelbase,tempsurf,tempbase -ts_file ts_y500.nc -ts_times 0:5:500 -o control_y500.nc

mpiexec -n 8 pismr -ocean_kill -e 3 -skip 5 \
  -ssa_sliding -thk_eff -pseudo_plastic_q 0.25 -plastic_pwfrac 0.98 -bed_def lc \
  -i g20km_0.nc -ocean constant \
  -atmosphere searise_greenland,dTforcing -surface pdd -pdd_fausto \
  -anomaly_temp ar4_temp_anomaly.nc \
  -anomaly_precip ar4_precip_anomaly.nc \
  -ys 0 -ye 500 -ts_file ts_ar4_y500.nc -ts_times 0:5:500 \
  -extra_file ex_ar4_y500.nc -extra_times 0:5:500 \
  -extra_vars csurf,cbase,usurf,topg,thk,bmelt,bwat,bwp,dHdt,mask,uvelsurf,vvelsurf,wvelsurf,uvelbase,vvelbase,wvelbase,tempsurf,tempbase \
  -o_size big -o ar4_y500.nc
  
