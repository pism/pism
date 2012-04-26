#!/bin/bash

# these

# if your results are different, try getting input files from Ed's ftp site:
#FIXME:  they were posted to dogbert but I can't find the index online
#wget -nc ???/pism_Greenland_v1.1.nc
#wget -nc ???/g20km_m5000a.nc

#  see http://www.pism-docs.org/wiki/doku.php?id=kspdiverged
#  recall that to look at the saved linear system   A v = b, do
##default options:
#mpiexec -n 4 ./ex10 -f SSAFD_kspdivergederror.petsc    # note 10000 iterations and large residual
#mpiexec -n 4 ./ex10 -f SSAFD_kspdivergederror.petsc -ksp_error_if_not_converged
##more robust iterative options:
#mpiexec -n 4 ./ex10 -f SSAFD_kspdivergederror.petsc -ksp_type gmres -ksp_norm_type unpreconditioned -ksp_pc_side right -pc_type asm -sub_pc_type lu
##direct solve only in serial:
#./ex10 -f SSAFD_kspdivergederror.petsc -ksp_type preonly -pc_type lu

echo "generates a ksp diverged that the scheme recovers from after one try:"
echo
mpiexec -n 4 pismr -config_override searise_config.nc -acab_cumulative -skip -skip_max 10 -boot_file pism_Greenland_5km_v1.1.nc -Mx 76 -My 141 -Lz 4000 -Lbz 2000 -Mz 101 -Mbz 11 -z_spacing equal -ssa_sliding -thk_eff -topg_to_phi 5.0,20.0,-300.0,700.0 -bed_def lc -atmosphere searise_greenland,delta_T -surface pdd -paleo_precip -atmosphere_delta_T_file pism_dT.nc -ocean constant,delta_SL -ocean_delta_SL_file pism_dSL.nc -ocean_kill -regrid_file g20km_m5000a.nc -regrid_vars litho_temp,thk,enthalpy,bwat,bmelt -regrid_bed_special -y 10

exit

echo "same run with better linear algebra options; no ksp diverged:"
echo

mpiexec -n 4 pismr -config_override searise_config.nc -acab_cumulative -skip -skip_max 10 -boot_file pism_Greenland_5km_v1.1.nc -Mx 76 -My 141 -Lz 4000 -Lbz 2000 -Mz 101 -Mbz 11 -z_spacing equal -ssa_sliding -thk_eff -topg_to_phi 5.0,20.0,-300.0,700.0 -bed_def lc -atmosphere searise_greenland,delta_T -surface pdd -paleo_precip -atmosphere_delta_T_file pism_dT.nc -ocean constant,delta_SL -ocean_delta_SL_file pism_dSL.nc -ocean_kill -regrid_file g20km_m5000a.nc -regrid_vars litho_temp,thk,enthalpy,bwat,bmelt -regrid_bed_special -y 10 -ksp_type gmres -ksp_norm_type unpreconditioned -ksp_pc_side right -pc_type asm -sub_pc_type lu

echo "generates a ksp diverged; scheme recovers from after 3 increases of epsilon"
echo
mpiexec -n 4 pismr -config_override searise_config.nc -acab_cumulative -skip -skip_max 10 -boot_file pism_Greenland_5km_v1.1.nc -Mx 76 -My 141 -Lz 4000 -Lbz 2000 -Mz 101 -Mbz 11 -z_spacing equal -ssa_sliding -thk_eff -topg_to_phi 5.0,20.0,-300.0,700.0 -bed_def lc -atmosphere searise_greenland,delta_T -surface pdd -ocean_kill -regrid_file g20km_m5000a.nc -regrid_vars litho_temp,thk,enthalpy,bwat,bmelt -y 10 -atmosphere_delta_T_file pism_dT.nc

