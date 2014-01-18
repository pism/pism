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
print ''

print 'accumrate=%s' % args.accumrate

# note My is ignored if args.e=='Stnd'
print ''
print '# grid'
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
if args.e in {'Stnd','P10S','P75S'}:
    print '# create bootstrap file using python ...'
else:
    print '# NOT creating bootstrap file since experiment %s starts from previously-saved state' % args.e

if args.e=='Stnd':
    print '%s createSetup_Stnd.py -a $accumrate -r $resolution' % args.pythonpath
elif args.e=='P10S':
    print 'amplitude=0.1'
    print '%s createSetup_PXXS.py -a $amplitude -i ex_Stnd.nc $s' % args.pythonpath
elif args.e=='P75S':
    print 'amplitude=0.75'
    print '%s createSetup_PXXS.py -a $amplitude -i ex_Stnd.nc $s' % args.pythonpath

print ''
print '# build the PISM command'
if args.n > 1:
    print 'pismr="%s"' % (args.mpiname + (' -n %d ' % args.n) + args.pismpath)
else:
    print 'pismr="%s"' % args.pismpath

print ''

print 'listexvar="thk,topg,cbar,cflx,mask,dHdt,usurf,hardav,velbase,velsurf,velbar,wvelbase,wvelsurf,deviatoric_stresses,climatic_mass_balance$gl_mask"'
if args.e == 'Stnd':
    print 'integration_time=3000'
    print 'extrastuff="-extra_times 0:50:$integration_time -extra_vars $listexvar"'
else:
    print 'integration_time=100'
    print 'extrastuff="-extra_times 0:1:$integration_time -extra_vars $listexvar"'

print ''
print 'stressbalance="-ssa_sliding -ssa_method fd -ssa_flow_law isothermal_glen -ssafd_ksp_rtol 1e-7"'
print 'basal="-yield_stress constant -pseudo_plastic -pseudo_plastic_q 0.333333333 -pseudo_plastic_uthreshold 3.155693e+07"'
print 'calvingfront="-cfbc -part_grid -calving ocean_kill"'
if args.m==1:
	print 'modelopt="-no_sia" '
elif args.m==2:
	print 'modelopt="-sia -sia_flow_law isothermal_glen" '

print ''
print 'opts="-config_override MISMIP3D_conf.nc $stressbalance $basal $calvingfront $subgl $modelopt -no_energy -cold -gradient eta -options_left -ts_file ts_%s.nc -ts_times 0:1:$integration_time -extra_file ex_%s.nc $extrastuff -ys 0 -ye $integration_time -o_order zyx -o_size big -o %s.nc"' % (args.e,args.e,args.e)

print ''
if args.e=='Stnd':
    print 'infile=MISMIP3D_Stnd_initialSetup.nc'
    print '$pismr -boot_file $infile -Mx $Mx -My 3 -Mz 15 -Lz 6000 -tauc 1.0e7 -ocean_kill_file $infile $opts'
elif args.e=='P10S':
    print 'infile=MISMIP3D_P10S_initialSetup.nc'
    print '$pismr -boot_file $infile -Mx $Mx -My $My -Mz 15 -Lz 6000 -ocean_kill_file $infile $opts'
elif args.e=='P10R':
    print 'infile=P10S.nc'
    print '$pismr -i $infile -tauc 1.0e7 -ocean_kill_file $infile $opts'
elif args.e=='P75S':
    print 'infile=MISMIP3D_P75S_initialSetup.nc'
    print '$pismr -boot_file $infile -Mx $Mx -My $My -Mz 15 -Lz 6000 -ocean_kill_file $infile $commonopts'
elif args.e=='P75R':
    print 'infile=P75S.nc'
    print '$pismr -i $infile -tauc 1.0e7 -ocean_kill_file $infile $opts'

