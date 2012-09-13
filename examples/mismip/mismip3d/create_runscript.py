#!/usr/bin/env python

# Copyright (C) 2012 Moritz Huetten and Torsten Albrecht

# create MISMIP-runscript

import sys
import getopt
try:
    from netCDF4 import Dataset as NC
except:
    from netCDF3 import Dataset as NC


# default values:

model=2
subgl=False
accumrate=0.5
resolutionmode=6

#### command line arguments ####
try:
  opts, args = getopt.getopt(sys.argv[1:], "m:sa:r:", ["model=","subgl","accumrate=","resolutionmode="])
  for opt, arg in opts:
    if opt in ("-m", "--model"):
      model = arg
    if opt in ("-s", "--subgl"):
      subgl = True
    if opt in ("-a", "--accumrate"):
      accumrate = arg 
    if opt in ("-r", "--resolutionmode"):
      resolutionmode = arg 

except getopt.GetoptError:
  print 'Incorrect command line arguments'
  sys.exit(2)

model=int(model)
accumrate=float(accumrate)
resolutionmode=int(resolutionmode)

print '###### MISMIP3D run Script for complete reversibility experiment ######'
print ''
print 'python="python"' #set the python path
print 'pismr="pismr"' #set the pismr path
print ''

if model==1:
	print '# model is SSA only'
	print 'modelopt="-no_sia" '
elif model==2:
	print '# model is hybrid SIA+SSA'
	print 'modelopt="-sia -sia_flow_law isothermal_glen" '
else:
 	print 'Incorrect command line arguments'
  	sys.exit(2)

print ''
print 'accumrate=%s' % accumrate
print ''


if resolutionmode==1:
	print 'resolution=0.5 # resolution in km'
	print 'Mx=3201'
	print 'My=201'
elif resolutionmode==2:
	print 'resolution=1 # resolution in km' 
	print 'Mx=1601'
	print 'My=101'
elif resolutionmode==3:
	print 'resolution=2.5 # resolution in km'
	print 'Mx=641'
	print 'My=41'
elif resolutionmode==4:
	print 'resolution=5 # resolution in km'
	print 'Mx=321'
	print 'My=21'
elif resolutionmode==5:
	print 'resolution=10 # resolution in km'
	print 'Mx=161'
	print 'My=11'
elif resolutionmode==6:
	print 'resolution=16.666 # resolution in km'
	print 'Mx=97'
	print 'My=7'
else:
 	print 'Incorrect command line arguments'
  	sys.exit(2)
print ''

print 'integration_time_stnd=3000'
print 'integration_time_PXXS=100'
print 'integration_time_PXXR=100'
print ''

# subgl usage:
if subgl==True:
	print '# subgrid grounding line interpolation is used'
	print 's="-s"'
	print 'subgl="-subgl"'
	print 'gl_mask=",gl_mask"'
else:
	print '# subgrid grounding line interpolation is not used'
	print 's=" "'
	print 'subgl=" "'
	print 'gl_mask=" "'
print ''


# create MISMIP config override file:


filename = "MISMIP3D_conf.nc"

nc = NC(filename, 'w', format="NETCDF3_CLASSIC")

var = nc.createVariable("pism_overrides", 'i')

attrs = {"is_dry_simulation" : "no",
                 "include_bmr_in_continuity" : "no",
                 "compute_surf_grad_inward_ssa" : "no",
                 "ice_softness" : 1.0e-25,
                 "ice_density" : 900.,
                 "sea_water_density" : 1000.,
                 "bootstrapping_geothermal_flux_value_no_var" : 0.0,
                 "Glen_exponent" : 3.,
                 "standard_gravity": 9.81,
                 "ocean_sub_shelf_heat_flux_into_ice" : 0.0,
                 "bed_smoother_range" : 0.0,
                 }

for name, value in attrs.iteritems():
            var.setncattr(name, value)

nc.close()


print '# Standard experiment:'
print '$python createSetup_Stnd.py -a $accumrate -r $resolution > Stnd.out'
print 'interval=$(($integration_time_stnd/50))'
print ''
print '$pismr -boot_file MISMIP3D_stnd_initialSetup.nc -Mx $Mx -My 3 -Mz 15 -Lz 6000 $subgl -cold $modelopt -ssa_flow_law isothermal_glen -no_energy -ssa_sliding  -pseudo_plastic -gradient eta -pseudo_plastic_q 0.333333333 -tauc 1.0e7 -hold_tauc -pseudo_plastic_uthreshold 3.155693e+07 -ocean_kill -config_override MISMIP3D_conf.nc -ssa_method fd -cfbc -part_grid -ksp_rtol 1e-7 -ys 0 -ye $integration_time_stnd -options_left -skip -skip_max 10 -stress_output -extra_file ex_Stnd.nc -extra_times 0:50:$integration_time_stnd -extra_vars thk,topg,cbar,cflx,mask,dHdt,usurf,hardav,velbase,velsurf,velbar,wvelbase,wvelsurf,sigma_xx,sigma_yy,sigma_xy,climatic_mass_balance$gl_mask -ts_file ts_Stnd.nc -ts_times 0:50:$integration_time_stnd -o Stnd.nc -o_order zyx -o_size big >> Stnd.out'


print ''
print '# P10S experiment:'
print 'amplitude=0.1'
print '$python createSetup_PXXS.py -a $amplitude -i ex_Stnd.nc $s > P10S.out'
print ''
print '$pismr -boot_file MISMIP3D_P10S_initialSetup.nc -Mx $Mx -My $My -Mz 15 -Lz 6000 $subgl -cold $modelopt -ssa_flow_law isothermal_glen -no_energy -ssa_sliding  -pseudo_plastic -gradient eta -pseudo_plastic_q 0.333333333 -hold_tauc -pseudo_plastic_uthreshold 3.155693e+07 -ocean_kill -config_override MISMIP3D_conf.nc -ssa_method fd -cfbc -part_grid -ksp_rtol 1e-7 -ys 0 -ye $integration_time_PXXS -options_left -skip -skip_max 10 -stress_output -o P10S.nc -o_size big -o_order zyx -extra_file ex_P10S.nc -extra_times 0:1:$integration_time_PXXS -extra_vars thk,topg,cbar,cflx,mask,dHdt,usurf,hardav,velbase,velsurf,velbar,wvelbase,wvelsurf,sigma_xx,sigma_yy,sigma_xy$gl_mask -ts_file ts_P10S.nc -ts_times 0:1:$integration_time_PXXS >>  P10S.out'

print ''
print '# P10R experiment:'
print ''
print '$pismr -i P10S.nc $subgl -cold $modelopt -ssa_flow_law isothermal_glen -no_energy -ssa_sliding  -pseudo_plastic -gradient eta -pseudo_plastic_q 0.333333333 -tauc 1.0e7 -hold_tauc -pseudo_plastic_uthreshold 3.155693e+07 -ocean_kill -config_override MISMIP3D_conf.nc -ssa_method fd -cfbc -part_grid -ksp_rtol 1e-7 -ys 0 -ye $integration_time_PXXR -options_left -skip -skip_max 10 -o_size big -o P10R.nc -o_order zyx -extra_file ex_P10R.nc -extra_times 0:1:100 -stress_output -extra_vars thk,topg,cbar,cflx,mask,dHdt,usurf,hardav,velbase,velsurf,velbar,wvelbase,wvelsurf,sigma_xx,sigma_yy,sigma_xy$gl_mask -ts_file ts_P10R.nc -ts_times 0:10:$integration_time_PXXR >>  P10R.out'

print ''
print '# P75S experiment:'
print 'amplitude=0.75'
print '$python createSetup_PXXS.py -a $amplitude -i ex_Stnd.nc $s > P75S.out'
print ''
print '$pismr -boot_file MISMIP3D_P75S_initialSetup.nc -Mx $Mx -My $My -Mz 15 -Lz 6000 $subgl -cold $modelopt -ssa_flow_law isothermal_glen -no_energy -ssa_sliding  -pseudo_plastic -gradient eta -pseudo_plastic_q 0.333333333 -hold_tauc -pseudo_plastic_uthreshold 3.155693e+07 -ocean_kill -config_override MISMIP3D_conf.nc -ssa_method fd -cfbc -part_grid -ksp_rtol 1e-7 -ys 0 -ye $integration_time_PXXS -options_left -skip -skip_max 10 -stress_output -o P75S.nc -o_size big -o_order zyx -extra_file ex_P75S.nc -extra_times 0:1:$integration_time_PXXS -extra_vars thk,topg,cbar,cflx,mask,dHdt,usurf,hardav,velbase,velsurf,velbar,wvelbase,wvelsurf,sigma_xx,sigma_yy,sigma_xy$gl_mask -ts_file ts_P75S.nc -ts_times 0:1:$integration_time_PXXS >>  P75S.out'

print ''
print '# P75R experiment:'
print ''
print '$pismr -i P75S.nc $subgl -cold $modelopt -ssa_flow_law isothermal_glen -no_energy -ssa_sliding  -pseudo_plastic -gradient eta -pseudo_plastic_q 0.333333333 -tauc 1.0e7 -hold_tauc -pseudo_plastic_uthreshold 3.155693e+07 -ocean_kill -config_override MISMIP3D_conf.nc -ssa_method fd -cfbc -part_grid -ksp_rtol 1e-7 -ys 0 -ye $integration_time_PXXR -options_left -skip -skip_max 10 -o_size big -o P10R.nc -o_order zyx -extra_file ex_P75R.nc -extra_times 0:1:100 -stress_output -extra_vars thk,topg,cbar,cflx,mask,dHdt,usurf,hardav,velbase,velsurf,velbar,wvelbase,wvelsurf,sigma_xx,sigma_yy,sigma_xy$gl_mask -ts_file ts_P75R.nc -ts_times 0:10:$integration_time_PXXR >>  P75R.out'

print ''
print '# P75D experiment based on Elmer data not done here yet'
print ''
