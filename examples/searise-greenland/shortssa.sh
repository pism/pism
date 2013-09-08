#!/bin/bash

DELTA=0.01
DELNAME=1percent

mpiexec -n 4 pismr -config_override searise_config.nc -till_effective_fraction_overburden $DELTA -skip -skip_max 10 -i g20km_steady.nc -ssa_sliding -topg_to_phi 15.0,40.0,-300.0,700.0 -bed_def lc -atmosphere searise_greenland,delta_T,paleo_precip -surface pdd -atmosphere_paleo_precip_file pism_dT.nc -atmosphere_delta_T_file pism_dT.nc -ocean constant,delta_SL -ocean_delta_SL_file pism_dSL.nc -ocean_kill pism_Greenland_5km_v1.1.nc -ts_file ts_g20km_shortssa_$DELNAME.nc -ts_times -125000:1:-124900 -extra_file ex_g20km_shortssa_$DELNAME.nc -extra_vars diffusivity,temppabase,tempicethk_basal,bmelt,tillwat,csurf,hardav,mask,dHdt,cbase,tauc,thk,topg,usurf,climatic_mass_balance_cumulative -extra_times -125000:1:-124900 -ys -125000 -ye -124900 -o g20km_shortssa_$DELNAME.nc
