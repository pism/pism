#!/bin/bash

#FIXME:  LOOK AT Jacobshavn CASE BEFORE REVISING THIS

#WHOLEFILE=g5km_0_ftt.nc   # 1 Gb or so
BCFILE=g5km_bc.nc          # 500Mb or so

wget -nc http://www.pism-docs.org/download/peter.nc

#wget -nc http://www.pism-docs.org/download/$WHOLEFILE
#ncks -v u_ssa,v_ssa,bmelt,enthalpy $WHOLEFILE $BCFILE
#ncrename -O -v u_ssa,u_ssa_bc -v v_ssa,v_ssa_bc $BCFILE
wget -nc http://www.pism-docs.org/download/$BCFILE


# FIXME: fails because searise_greenland needs lat,lon
mpiexec -n 4 pismo -boot_file peter.nc -Mx 99 -My 210 -no_model_strip 5 -Lz 4000 -Lbz 1000 -Mz 201 -Mbz 51 -z_spacing equal -atmosphere searise_greenland -surface pdd,forcing -force_to_thk peter.nc -ocean_kill -ssa_sliding -topg_to_phi 5.0,20.0,-300.0,700.0 -diffuse_bwat -thk_eff -plastic_pwfrac 0.95 -pseudo_plastic_q 0.25 -ssa_dirichlet_bc -regrid_file g5km_bc.nc -regrid_vars bmelt,bwat,enthalpy,litho_temp,vel_ssa_bc -y 10 -skip 50 -o peter3km_y10.nc 
