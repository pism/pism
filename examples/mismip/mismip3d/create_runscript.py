#!/usr/bin/env python

# Copyright (C) 2012-2014 Moritz Huetten and Torsten Albrecht (and Ed Bueler)

import sys, argparse

# process command line arguments
parser = argparse.ArgumentParser(description='Create run script for a MISMIP3d experiment.',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-a', '--accumrate', metavar='A', type=float, default=0.5,
                   help='accumulation rate in meters/year')
parser.add_argument('-e', default='Stnd',
                   choices=['Stnd','P10S','P10R','P75S','P75R'],  # P75D using Elmer is
                                                                  # not implemented yet
                   help='name of experiment')
parser.add_argument('-m', choices=[1,2], type=int, default=2,
                   help='model; 1 = SSA only, 2 = hybrid SIA+SSA')
parser.add_argument('--mpiname', metavar='NAME', default='mpiexec',
                   help='name of mpi executable')
parser.add_argument('-n', metavar='N', type=int, default=2,
                   help='number of MPI processes')
parser.add_argument('--pismpath', metavar='PATH', default='pismr',
                   help='full path to PISM executable pismr')
parser.add_argument('--pythonpath', metavar='PATH', default='python',
                   help='full path to python executable')
parser.add_argument('-r', type=int, choices=[1,2,3,4,5,6], default=6,
                   help='resolution mode; 1 = finest, 6 = coarsest')
parser.add_argument('-s', '--subgl', action='store_true', # thus defaults to False
                   help='use sub-grid grounding line method')
                    
args = parser.parse_args()
#print args   # helpful for debugging

print '#!/bin/bash'
print ''
print '###### MISMIP3D run script for experiment %s, model %d, and resolution mode %d ######' \
      % (args.e, args.m, args.r)

if args.n > 1:
    print 'pismr="%s"' % (args.mpiname + (' -n %d ' % args.n) + args.pismpath)
else:
    print 'pismr="%s"' % args.pismpath

print 'python="%s"' % args.pythonpath

if args.m==1:
	print 'modelopt="-no_sia" '
elif args.m==2:
	print 'modelopt="-sia -sia_flow_law isothermal_glen" '

print 'accumrate=%s' % args.accumrate

# note My is ignored if args.e=='Stnd'
print ''
if args.r==1:
	print 'resolution=0.5 # resolution in km'
	print 'Mx=3201'
	print 'My=201'
elif args.r==2:
	print 'resolution=1 # resolution in km' 
	print 'Mx=1601'
	print 'My=101'
elif args.r==3:
	print 'resolution=2.5 # resolution in km'
	print 'Mx=641'
	print 'My=41'
elif args.r==4:
	print 'resolution=5 # resolution in km'
	print 'Mx=321'
	print 'My=21'
elif args.r==5:
	print 'resolution=10 # resolution in km'
	print 'Mx=161'
	print 'My=11'
elif args.r==6:
	print 'resolution=16.666 # resolution in km'
	print 'Mx=97'
	print 'My=7'

print ''

if args.subgl:
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

print 'listexvar="thk,topg,cbar,cflx,mask,dHdt,usurf,hardav,velbase,velsurf,velbar,wvelbase,wvelsurf,deviatoric_stresses,climatic_mass_balance$gl_mask"'

FIXME

print 'commonopts="$subgl -cold $modelopt -ssa_flow_law isothermal_glen -no_energy -ssa_sliding -pseudo_plastic -gradient eta  -ocean_kill -config_override MISMIP3D_conf.nc -ssa_method fd -cfbc -part_grid -ssafd_ksp_rtol 1e-7"'

if args.e=='Stnd':
    print 'integration_time=3000'
    print '$python createSetup_Stnd.py -a $accumrate -r $resolution'
    print 'interval=$(($integration_time/50))'
    print ''
    print '$pismr -boot_file MISMIP3D_stnd_initialSetup.nc -Mx $Mx -My 3 -Mz 15 -Lz 6000 $subgl -cold $modelopt -ssa_flow_law isothermal_glen -no_energy -ssa_sliding  -pseudo_plastic -gradient eta -pseudo_plastic_q 0.333333333 -tauc 1.0e7 -hold_tauc -pseudo_plastic_uthreshold 3.155693e+07 -ocean_kill -config_override MISMIP3D_conf.nc -ssa_method fd -cfbc -part_grid -ssafd_ksp_rtol 1e-7 -ys 0 -ye $integration_time -options_left -extra_file ex_Stnd.nc -extra_times 0:50:$integration_time -extra_vars $listexvar -ts_file ts_Stnd.nc -ts_times 0:50:$integration_time -o Stnd.nc -o_order zyx -o_size big'

elif args.e=='P10S':
    print 'integration_time=100'
    print 'amplitude=0.1'
    print '$python createSetup_PXXS.py -a $amplitude -i ex_Stnd.nc $s'
    print ''
    print '$pismr -boot_file MISMIP3D_P10S_initialSetup.nc -Mx $Mx -My $My -Mz 15 -Lz 6000 $subgl -cold $modelopt -ssa_flow_law isothermal_glen -no_energy -ssa_sliding  -pseudo_plastic -gradient eta -pseudo_plastic_q 0.333333333 -hold_tauc -pseudo_plastic_uthreshold 3.155693e+07 -ocean_kill -config_override MISMIP3D_conf.nc -ssa_method fd -cfbc -part_grid -ssafd_ksp_rtol 1e-7 -ys 0 -ye $integration_time -options_left -stress_output -o P10S.nc -o_size big -o_order zyx -extra_file ex_P10S.nc -extra_times 0:1:$integration_time -extra_vars $listexvar -ts_file ts_P10S.nc -ts_times 0:1:$integration_time'

elif args.e=='P10R':
    print 'integration_time=100'
    print ''
    print '$pismr -i P10S.nc $subgl -cold $modelopt -ssa_flow_law isothermal_glen -no_energy -ssa_sliding  -pseudo_plastic -gradient eta -pseudo_plastic_q 0.333333333 -tauc 1.0e7 -hold_tauc -pseudo_plastic_uthreshold 3.155693e+07 -ocean_kill -config_override MISMIP3D_conf.nc -ssa_method fd -cfbc -part_grid -ssafd_ksp_rtol 1e-7 -ys 0 -ye $integration_time -options_left -o_size big -o P10R.nc -o_order zyx -extra_file ex_P10R.nc -extra_times 0:1:$integration_time -stress_output -extra_vars $listexvar -ts_file ts_P10R.nc -ts_times 0:10:$integration_time'

elif args.e=='P75S':
    print 'integration_time=100'
    print 'amplitude=0.75'
    print '$python createSetup_PXXS.py -a $amplitude -i ex_Stnd.nc $s'
    print ''
    print '$pismr -boot_file MISMIP3D_P75S_initialSetup.nc -Mx $Mx -My $My -Mz 15 -Lz 6000 $subgl -cold $modelopt -ssa_flow_law isothermal_glen -no_energy -ssa_sliding  -pseudo_plastic -gradient eta -pseudo_plastic_q 0.333333333 -hold_tauc -pseudo_plastic_uthreshold 3.155693e+07 -ocean_kill -config_override MISMIP3D_conf.nc -ssa_method fd -cfbc -part_grid -ssafd_ksp_rtol 1e-7 -ys 0 -ye $integration_time -options_left -stress_output -o P75S.nc -o_size big -o_order zyx -extra_file ex_P75S.nc -extra_times 0:1:$integration_time -extra_vars $listexvar -ts_file ts_P75S.nc -ts_times 0:1:$integration_time'

elif args.e=='P75R':
    print 'integration_time=100'
    print ''
    print '$pismr -i P75S.nc $subgl -cold $modelopt -ssa_flow_law isothermal_glen -no_energy -ssa_sliding  -pseudo_plastic -gradient eta -pseudo_plastic_q 0.333333333 -tauc 1.0e7 -hold_tauc -pseudo_plastic_uthreshold 3.155693e+07 -ocean_kill -config_override MISMIP3D_conf.nc -ssa_method fd -cfbc -part_grid -ssafd_ksp_rtol 1e-7 -ys 0 -ye $integration_time -options_left -o_size big -o P10R.nc -o_order zyx -extra_file ex_P75R.nc -extra_times 0:1:$integration_time -stress_output -extra_vars $listexvar -ts_file ts_P75R.nc -ts_times 0:10:$integration_time'

