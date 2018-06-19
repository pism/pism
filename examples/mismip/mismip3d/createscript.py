#!/usr/bin/env python

# Copyright (C) 2012-2015 Moritz Huetten and Torsten Albrecht (and Ed Bueler)

import argparse

# process command line arguments
parser = argparse.ArgumentParser(description='Create run script for a MISMIP3d experiment.',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-a', '--accumrate', metavar='A', type=float, default=0.5,
                    help='accumulation rate in meters/year')
parser.add_argument('-d', '--duration', metavar='T', type=float, default=-1.0,
                    help='duration of run in years (if not set, use 3000 years for Stnd and 100 for others)')
parser.add_argument('-e', default='Stnd',
                    choices=['Stnd', 'P10S', 'P10R', 'P75S', 'P75R'],
                    help='name of experiment; note P10D and P75D are not implemented yet')
parser.add_argument('-m', choices=[1, 2], type=int, default=2,
                    help='model; 1 = SSA only, 2 = hybrid SIA+SSA')
parser.add_argument('--mpiname', metavar='NAME', default='mpiexec',
                    help='name of mpi executable')
parser.add_argument('-n', metavar='N', type=int, default=2,
                    help='number of MPI processes; if N=1 then MPI is not used')
parser.add_argument('--pismpath', metavar='PATH', default='pismr',
                    help='full path to PISM executable pismr')
parser.add_argument('--pythonpath', metavar='PATH', default='python',
                    help='full path to python executable')
parser.add_argument('-r', type=int, choices=[1, 2, 3, 4, 5, 6, 7], default=5,
                    help='resolution mode; 1 = finest, 6 = coarsest')
parser.add_argument('-s', '--subgl', action='store_true',  # thus defaults to False
                    help='use sub-grid grounding line method')

args = parser.parse_args()
# print args   # helpful for debugging

print("""#!/bin/bash
###### MISMIP3D run script for experiment %s, model %d, and resolution mode %d ######

accumrate=%s""" % (args.e, args.m, args.r, args.accumrate))

# key: resolution mode; parameters: (resolution in km, Mx, My)
grid_parameters = {1: (0.5,    3201, 201),
                   2: (1,      1601, 101),
                   3: (2,      801,  51),
                   4: (2.5,    641,  41),
                   5: (5,      321,  21),
                   6: (10,     161,  11),
                   7: (16.666, 97,   7)}

# note My is ignored if args.e=='Stnd'
print("""
# grid
resolution=%f # resolution in km
Mx=%d
My=%d""" % grid_parameters[args.r])

print('')
if args.subgl:
    print('# subgrid grounding line interpolation is used')
    print('s="-s"')
    print('subgl="-subgl"')
    print('gl_mask=",gl_mask"')
else:
    print('# subgrid grounding line interpolation is not used')
    print('s=" "')
    print('subgl=" "')
    print('gl_mask=" "')

print('')
if args.e in {'Stnd', 'P10S', 'P75S'}:
    print('# create bootstrap file using python ...')
else:
    print('# NOT creating bootstrap file since experiment %s starts from previously-saved state' % args.e)

if args.e == 'Stnd':
    print('%s setup_Stnd.py -a $accumrate -r $resolution' % args.pythonpath)
elif args.e == 'P10S':
    print('amplitude=0.1')
    print('%s setup_PXXS.py -a $amplitude -i ex_Stnd.nc $s' % args.pythonpath)
elif args.e == 'P75S':
    print('amplitude=0.75')
    print('%s setup_PXXS.py -a $amplitude -i ex_Stnd.nc $s' % args.pythonpath)

print('')
print('# build the PISM command')
if args.n > 1:
    print('pismr="%s"' % (args.mpiname + (' -n %d ' % args.n) + args.pismpath))
else:
    print('pismr="%s"' % args.pismpath)

print('')
if args.duration < 0:
    if args.e == 'Stnd':
        print('duration=3000')
    else:
        print('duration=100')
else:
    print('duration=%s' % args.duration)

print('')
print('listexvar="thk,topg,velbar_mag,flux_mag,mask,dHdt,usurf,hardav,velbase,velsurf,velbar,wvelbase,wvelsurf,tauc,deviatoric_stresses,climatic_mass_balance$gl_mask"')
if args.e == 'Stnd':
    print('extrastuff="-extra_times 0:50:$duration -extra_vars $listexvar"')
else:
    print('extrastuff="-extra_times 0:1:$duration -extra_vars $listexvar"')

print('')
print('stressbalance="-ssa_method fd -ssa_flow_law isothermal_glen -ssafd_ksp_rtol 1e-7"')
print('basal="-yield_stress constant -pseudo_plastic -pseudo_plastic_q 0.333333333 -pseudo_plastic_uthreshold 3.155693e+07"')
print('calvingfront="-cfbc -part_grid -calving ocean_kill"')
if args.m == 1:
    print('modelopt="-stress_balance ssa" ')
elif args.m == 2:
    print('modelopt="-stress_balance ssa+sia -sia_flow_law isothermal_glen" ')

print('STRONGKSP="-ssafd_ksp_type gmres -ssafd_ksp_norm_type unpreconditioned -ssafd_ksp_pc_side right -ssafd_pc_type asm -ssafd_sub_pc_type lu"')

print('')
print('opts="-config_override MISMIP3D_conf.nc $stressbalance $basal $calvingfront $subgl $modelopt -energy none -gradient eta -options_left -ts_file ts_%s.nc -ts_times 0:1:$duration -extra_file ex_%s.nc $extrastuff -ys 0 -ye $duration -o_order zyx -o_size big -o %s.nc $STRONGKSP"' % (args.e, args.e, args.e))

print('')
if args.e == 'Stnd':
    print('infile=MISMIP3D_Stnd_initialSetup.nc')
    print('cmd="$pismr -i $infile -bootstrap -Mx $Mx -My 3 -Mz 15 -Lz 6000 -tauc 1.0e7 -ocean_kill_file $infile $opts"')
elif args.e == 'P10S':
    print('infile=MISMIP3D_P10S_initialSetup.nc')
    print('cmp="$pismr -i $infile -bootstrap -Mx $Mx -My $My -Mz 15 -Lz 6000 -ocean_kill_file $infile $opts"')
elif args.e == 'P10R':
    print('infile=P10S.nc')
    print('cmd="$pismr -i $infile -tauc 1.0e7 -ocean_kill_file $infile $opts"')
elif args.e == 'P75S':
    print('infile=MISMIP3D_P75S_initialSetup.nc')
    print('cmd="$pismr -i $infile -bootstrap -Mx $Mx -My $My -Mz 15 -Lz 6000 -ocean_kill_file $infile $opts"')
elif args.e == 'P75R':
    print('infile=P75S.nc')
    print('cmd="$pismr -i $infile -tauc 1.0e7 -ocean_kill_file $infile $opts"')

print('')
print('echo "running command:"')
print('echo $cmd')
print('echo')

print('')
print('$cmd')
