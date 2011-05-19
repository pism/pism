#!/bin/bash

# Copyright (C) 2011 The PISM Authors

# a version of SeaRISE-Greenland spinup that demonstrates PISM-PIK handling of
#   ice shelves, using shortened runs; all on 10km grid

#(spinup.sh)              NN = 8
#(spinup.sh)      PISM_MPIDO = mpiexec -n 

#(spinup.sh)     grid = '-Mx 151 -My 281 -Lz 4000 -Lbz 2000 -Mz 101 -Mbz 41' (= 10 km)
#(spinup.sh)       fine grid = '-Mx 151 -My 281 -Lz 4000 -Lbz 2000 -Mz 101 -Mbz 41' (= 10 km)
#(spinup.sh)      executable = 'pismr -e 3'  # no -ocean_kill
#(spinup.sh)    full physics = '-ssa_sliding -thk_eff -pseudo_plastic_q 0.25 -plastic_pwfrac 0.98 -topg_to_phi 5.0,20.0,-300.0,700.0,10.0'
#(spinup.sh)  simple coupler = '-atmosphere searise_greenland -surface pdd -config_override config_269.0_0.001_0.80_-0.500_9.7440.nc'
#(spinup.sh) forcing coupler = '-atmosphere searise_greenland,forcing -surface pdd -config_override config_269.0_0.001_0.80_-0.500_9.7440.nc -paleo_precip -dTforcing pism_dT.nc -ocean constant,forcing -dSLforcing pism_dSL.nc'

#(spinup.sh)  bootstrapping plus short smoothing run (for 100a)
mpiexec -n 8 pismr -e 3 -skip 20 -boot_file pism_Greenland_5km_v1.1.nc -Mx 151 -My 281 -Lz 4000 -Lbz 2000 -Mz 101 -Mbz 41 -atmosphere searise_greenland -surface pdd -config_override config_269.0_0.001_0.80_-0.500_9.7440.nc -y 100 -o g10km_pre100.nc -pik -eigen_calving

#(spinup.sh)  -no_mass (no surface change) SIA run to achieve approximate temperature equilibrium, for 20000a
mpiexec -n 8 pismr -e 3 -skip 20 -i g10km_pre100.nc -atmosphere searise_greenland -surface pdd -config_override config_269.0_0.001_0.80_-0.500_9.7440.nc -no_mass -y 20000 -extra_file ex_g10km_steady.nc -extra_vars diffusivity,temppabase,bmelt,csurf,hardav,mask -extra_times 0:500:20000 -o g10km_steady.nc

#(spinup.sh)  smoothing with SIA for 100a
mpiexec -n 8 pismr -e 3 -skip 20 -i g10km_steady.nc -atmosphere searise_greenland -surface pdd -config_override config_269.0_0.001_0.80_-0.500_9.7440.nc -y 100 -o g10km_SIA.nc -pik -eigen_calving

#(spinup.sh)  paleo-climate forcing run with full physics,
#(spinup.sh)      except bed deformation, from -40000 a to -20000a
mpiexec -n 8 pismr -e 3 -skip 20 -i g10km_SIA.nc -ssa_sliding -thk_eff -pseudo_plastic_q 0.25 -plastic_pwfrac 0.98 -topg_to_phi 5.0,20.0,-300.0,700.0,10.0 -atmosphere searise_greenland,dTforcing -surface pdd -config_override config_269.0_0.001_0.80_-0.500_9.7440.nc -paleo_precip -dTforcing pism_dT.nc -ocean constant,dSLforcing -dSLforcing pism_dSL.nc -ts_file ts_g10km_m20ka.nc -ts_times -40000:1:-20000 -extra_file ex_g10km_m20ka.nc -extra_vars thk,usurf,diffusivity,tempthk_basal,temppabase,bmelt,csurf,hardav,mask,cbase,tauc,IcebergMask -extra_times -124500:500:-20000 -ys -40000 -ye -20000 -o g10km_m20ka.nc -pik -eigen_calving

#(spinup.sh)  regrid to fine grid and do paleo-climate forcing run with full physics,
#(spinup.sh)      including bed deformation, from -20000a BPE to 0a BPE
mpiexec -n 8 pismr -e 3 -skip 20 -boot_file pism_Greenland_5km_v1.1.nc -Mx 151 -My 281 -Lz 4000 -Lbz 2000 -Mz 101 -Mbz 41 -ssa_sliding -thk_eff -pseudo_plastic_q 0.25 -plastic_pwfrac 0.98 -topg_to_phi 5.0,20.0,-300.0,700.0,10.0 -bed_def lc -atmosphere searise_greenland,dTforcing -surface pdd -config_override config_269.0_0.001_0.80_-0.500_9.7440.nc -paleo_precip -dTforcing pism_dT.nc -ocean constant,dSLforcing -dSLforcing pism_dSL.nc -regrid_file g10km_m20ka.nc -regrid_vars litho_temp,thk,enthalpy,bwat -ts_file ts_g10km_0.nc -ts_times -20000:1:0 -extra_file ex_g10km_0.nc -extra_vars thk,usurf,diffusivity,tempthk_basal,temppabase,bmelt,csurf,hardav,mask,cbase,tauc,IcebergMask -extra_times -19500:500:0 -ys -20000 -ye 0 -o g10km_0.nc -pik -eigen_calving

#(spinup.sh)  spinup done
