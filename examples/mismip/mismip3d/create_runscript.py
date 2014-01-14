#!/usr/bin/env python

# Copyright (C) 2012 Moritz Huetten and Torsten Albrecht

# create MISMIP-runscript

import sys
import getopt


# default values:

model=2
subgl=False
accumrate=0.5
resolutionmode=6
exper='Stnd'

#### command line arguments ####
try:
  opts, args = getopt.getopt(sys.argv[1:], "m:sa:r:e:",
                             ["model=","subgl","accumrate=","resolutionmode=","experiment="])
  for opt, arg in opts:
    if opt in ("-m", "--model"):
      model = arg
    if opt in ("-s", "--subgl"):
      subgl = True
    if opt in ("-a", "--accumrate"):
      accumrate = arg 
    if opt in ("-r", "--resolutionmode"):
      resolutionmode = arg 
    if opt in ("-e", "--experiment"):
      exper = arg

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


if exper=='Stnd':
    print 'integration_time=3000'
    print ''
    print '# Standard experiment:'
    print '$python createSetup_Stnd.py -a $accumrate -r $resolution'
    print 'interval=$(($integration_time/50))'
    print ''
    print '$pismr -boot_file MISMIP3D_stnd_initialSetup.nc -Mx $Mx -My 3 -Mz 15 -Lz 6000 $subgl -cold $modelopt -ssa_flow_law isothermal_glen -no_energy -ssa_sliding  -pseudo_plastic -gradient eta -pseudo_plastic_q 0.333333333 -tauc 1.0e7 -hold_tauc -pseudo_plastic_uthreshold 3.155693e+07 -ocean_kill -config_override MISMIP3D_conf.nc -ssa_method fd -cfbc -part_grid -ksp_rtol 1e-7 -ys 0 -ye $integration_time -options_left -skip -skip_max 10 -extra_file ex_Stnd.nc -extra_times 0:50:$integration_time -extra_vars thk,topg,cbar,cflx,mask,dHdt,usurf,hardav,velbase,velsurf,velbar,wvelbase,wvelsurf,deviatoric_stresses,climatic_mass_balance$gl_mask -ts_file ts_Stnd.nc -ts_times 0:50:$integration_time -o Stnd.nc -o_order zyx -o_size big'

elif exper=='P10S':
    print 'integration_time=100'
    print ''
    print '# P10S experiment:'
    print 'amplitude=0.1'
    print '$python createSetup_PXXS.py -a $amplitude -i ex_Stnd.nc $s'
    print ''
    print '$pismr -boot_file MISMIP3D_P10S_initialSetup.nc -Mx $Mx -My $My -Mz 15 -Lz 6000 $subgl -cold $modelopt -ssa_flow_law isothermal_glen -no_energy -ssa_sliding  -pseudo_plastic -gradient eta -pseudo_plastic_q 0.333333333 -hold_tauc -pseudo_plastic_uthreshold 3.155693e+07 -ocean_kill -config_override MISMIP3D_conf.nc -ssa_method fd -cfbc -part_grid -ksp_rtol 1e-7 -ys 0 -ye $integration_time -options_left -skip -skip_max 10 -stress_output -o P10S.nc -o_size big -o_order zyx -extra_file ex_P10S.nc -extra_times 0:1:$integration_time -extra_vars thk,topg,cbar,cflx,mask,dHdt,usurf,hardav,velbase,velsurf,velbar,wvelbase,wvelsurf,sigma_xx,sigma_yy,sigma_xy$gl_mask -ts_file ts_P10S.nc -ts_times 0:1:$integration_time'

elif exper=='P10R':
    print 'integration_time=100'
    print ''
    print '# P10R experiment:'
    print ''
    print '$pismr -i P10S.nc $subgl -cold $modelopt -ssa_flow_law isothermal_glen -no_energy -ssa_sliding  -pseudo_plastic -gradient eta -pseudo_plastic_q 0.333333333 -tauc 1.0e7 -hold_tauc -pseudo_plastic_uthreshold 3.155693e+07 -ocean_kill -config_override MISMIP3D_conf.nc -ssa_method fd -cfbc -part_grid -ksp_rtol 1e-7 -ys 0 -ye $integration_time -options_left -skip -skip_max 10 -o_size big -o P10R.nc -o_order zyx -extra_file ex_P10R.nc -extra_times 0:1:$integration_time -stress_output -extra_vars thk,topg,cbar,cflx,mask,dHdt,usurf,hardav,velbase,velsurf,velbar,wvelbase,wvelsurf,sigma_xx,sigma_yy,sigma_xy$gl_mask -ts_file ts_P10R.nc -ts_times 0:10:$integration_time'

elif exper=='P75S':
    print 'integration_time=100'
    print ''
    print '# P75S experiment:'
    print 'amplitude=0.75'
    print '$python createSetup_PXXS.py -a $amplitude -i ex_Stnd.nc $s'
    print ''
    print '$pismr -boot_file MISMIP3D_P75S_initialSetup.nc -Mx $Mx -My $My -Mz 15 -Lz 6000 $subgl -cold $modelopt -ssa_flow_law isothermal_glen -no_energy -ssa_sliding  -pseudo_plastic -gradient eta -pseudo_plastic_q 0.333333333 -hold_tauc -pseudo_plastic_uthreshold 3.155693e+07 -ocean_kill -config_override MISMIP3D_conf.nc -ssa_method fd -cfbc -part_grid -ksp_rtol 1e-7 -ys 0 -ye $integration_time -options_left -skip -skip_max 10 -stress_output -o P75S.nc -o_size big -o_order zyx -extra_file ex_P75S.nc -extra_times 0:1:$integration_time -extra_vars thk,topg,cbar,cflx,mask,dHdt,usurf,hardav,velbase,velsurf,velbar,wvelbase,wvelsurf,sigma_xx,sigma_yy,sigma_xy$gl_mask -ts_file ts_P75S.nc -ts_times 0:1:$integration_time'

elif exper=='P75R':
    print 'integration_time=100'
    print ''
    print '# P75R experiment:'
    print ''
    print '$pismr -i P75S.nc $subgl -cold $modelopt -ssa_flow_law isothermal_glen -no_energy -ssa_sliding  -pseudo_plastic -gradient eta -pseudo_plastic_q 0.333333333 -tauc 1.0e7 -hold_tauc -pseudo_plastic_uthreshold 3.155693e+07 -ocean_kill -config_override MISMIP3D_conf.nc -ssa_method fd -cfbc -part_grid -ksp_rtol 1e-7 -ys 0 -ye $integration_time -options_left -skip -skip_max 10 -o_size big -o P10R.nc -o_order zyx -extra_file ex_P75R.nc -extra_times 0:1:$integration_time -stress_output -extra_vars thk,topg,cbar,cflx,mask,dHdt,usurf,hardav,velbase,velsurf,velbar,wvelbase,wvelsurf,sigma_xx,sigma_yy,sigma_xy$gl_mask -ts_file ts_P75R.nc -ts_times 0:10:$integration_time'

elif exper=='P75D':
    print ''
    print '# P75D experiment based on Elmer data not done here yet'
    print 'ERROR: NO RUN'

else:
    print ''
    print '# ERROR: UNKNOWN EXPERIMENT NUMBER'
    print 'ERROR: NO RUN'
