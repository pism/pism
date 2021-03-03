#!/bin/bash

### This run script for a flowline ice shelf setup
#   tests the Calving Front Boundary Condition at ice shelves
#   and a modified driving stress scheme at the calving front.
#
#   contact: reese@pik-potsdam.de and albrecht@pik-potsdam.de"""


### settings ###############################

# PISM version
#pismex=/home/albrecht/pism19/pism1.1/bin/pismr
pismex=~/pism20/pism1.2.2_cfbc/bin/pismr

# initial ice thickness and velocity in m and m/yr
H0=800
v0=300

# calving front position in m
ca=1750e3

# resolution for 1800km length 
N0=601
#N0=3001

# PISM output file names
expname=${H0}_${v0}_${ca}_${N0}
of=init_${expname}.nc
op=initp_${expname}.nc

rd=result_${expname}_dev.nc
rf=result_${expname}_glfix.nc

rdp=resultp_${expname}_dev.nc
rfp=resultp_${expname}_glfix.nc

echo $of

# build PISM initial setup
python prepare_iceshelf.py -c $ca -v $v0 -H $H0 -o $of -N $N0
python prepare_iceshelf.py -c $ca -v $v0 -H $H0 -o $op -N $N0 --perturb_cf

#mkdir results
ncgen3 confover.cdl -o confover.nc

runopts="-bootstrap -Mx ${N0} -My 3 -Mz 15 -Lz 6000 -energy none -ssa_flow_law isothermal_glen -config_override confover.nc -ssa_method fd -cfbc -part_grid -ys 0 -ye 0.0001 -options_left -stress_balance ssa -ssa_dirichlet_bc -no_mass -shelf_base_melt_rate 0.0"
#runopts="-bootstrap -Mx ${N0} -My 3 -Mz 15 -Lz 6000 -energy none -ssa_flow_law isothermal_glen -yield_stress constant -tauc 7.624000e+06 -pseudo_plastic -pseudo_plastic_q 3.333333e-01 -pseudo_plastic_uthreshold 3.155693e+07 -config_override confover.nc -ssa_method fd -cfbc -part_grid -ssafd_ksp_rtol 1e-7 -ys 0 -ye 0.0001 -options_left -stress_balance ssa -ssa_dirichlet_bc -no_mass -shelf_base_melt_rate 0.0"


# PISM run options
$pismex -i $of $runopts -o $rd
$pismex -i $of $runopts -o $rf -taudcf

$pismex -i $op $runopts -o $rdp
$pismex -i $op $runopts -o $rfp -taudcf


# plot PISM output
python plot_iceshelf.py $rd $rf $of
python plot_iceshelf.py $rdp $rfp $op
