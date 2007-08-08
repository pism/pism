#! /usr/bin/env python

## ANTNC.PY  is designed for modelling the current state of Antarctica on a 
## multi-processor machine.  Ice streams are determined by using balance velocities
## OR Schoof's plastic till model.  Assumes  init.nc  is present (which contains 5km 
## data mostly from the British Antarctica Service.)  Tends to require many processors 
## just to have memory for model state.  Contact ffelb@uaf.edu with questions.
## ELB 1/31/07; 5/29/07: .sh --> .py

import os
import sys
import time
import string

## default settings
NP = 8
BINITFILE = 'init.nc'

## program prog is a list; elements of prog are:
##     ['DESCRIPTION','EXECUTABLE',YEARS,INITFLAG,'OUTNAME','OTHER_OPTIONS'[,'REGRIDFILENAME','REGRIDVARS']]
## note: INITFLAG == 0   <==>   '-bif BINITFILE'
##       INITFLAG == -1   <==>   '-if [prev OUTNAME]'          {case of next}
##       INITFLAG < 0    <==>   '-if [OUTNAME from stage[current-|INITFLAG|]]'
prog = [
['input BAS data, put on 40km grid, and smooth surface',
      'pismr',    10,0,'ant10yr_40km','-Mx 141 -My 141 -Mz 101 -Mbz 21 -verbose'],
['running for 155k year to equilibriate temp (steady geometry; SIA only)',
      'pismr',154990,-1,'ant155k_40km','-no_mass'],
['regrid to 20km grid: interpolate temp and age and smooth surface',
      'pismr',    10,0,'ant155k_20km','-Mx 281 -My 281 -Mz 201 -Mbz 41 -verbose',
      'ant155k_40km','TBel'],
['running for another 3000 yrs to equilibriate temp on fine grid (steady geometry and SIA only) ',
      'pismr',  2990,-1,'ant158k_20km','-no_mass'],
['2 yrs to smooth (with MacAyeal)',
      'pismr',     2,-1,'ant158k_20kmSMMV','-mv -ksp_rtol 1e-6 -ocean_kill -verbose'],
['another 500 yrs to equilibriate temp (steady geometry with MacAyeal)',
      'pismr',   498,-1,'ant158p5k_20km','-no_mass -mv -ksp_rtol 1e-6 -ocean_kill -verbose'],
['regrid to 14km grid and smooth for 2 years',
      'pismr',     2,0,'ant158p5k_14km','-Mx 401 -My 401 -Mz 201 -Mbz 41 -gk -e 1.2 -verbose',
      'ant158p5k_20km','TBeL'],
['compute everything for 5 years on 14km grid',
      'pismr',     5,-1,'ant158p5k_p5yr_14km','-mv -ksp_rtol 1e-7 -ocean_kill -verbose'],
['pant using scalar beta on 14km grid',
      'pant',      5,-2,'antlast14km_beta',
      '-mv -ksp_rtol 1e-7 -ocean_kill -verbose -beta -balvel init.nc -obasal betamap14km'],
['pant using vector beta on 14km grid',
      'pant',      5,-3,'antlast14km_betaxy',
      '-mv -ksp_rtol 1e-7 -ocean_kill -verbose -betaxy -balvel init.nc -obasal betamapxy14km'],
['pant using plastic on 14km grid',
      'pant',      5,-4,'antlast14km_tauc',
      '-mv -ksp_rtol 1e-7 -ocean_kill -verbose -super -plastic -obasal taucmap14km']
]
always_opts = '-gk -e 1.2'

## get options: -n for number of processors
nproc = NP  ## default; will not use 'mpiexec' if equal to one
for ll in range(len(sys.argv))[1:]:
   if sys.argv[ll].find('n') > -1:
      nproc = string.atoi(sys.argv[ll+1])

## run program: stage[] above determines stages of computation
print '  ANTNC using %d processor(s)' % nproc
stage_count = 0
for stage in prog:
   print ('  stage %d.  ' % (stage_count+1)) + stage[0] + ':'
   ## use mpi only if nproc > 1
   stagedo = ''
   if nproc > 1:  stagedo += 'mpiexec -n %d' % nproc
   stagedo += ' ' + stage[1] + (' -y %d' % stage[2]) 
   ## determine infile
   if stage[3] == 0:
      stagedo += ' -bif ' + BINITFILE
   elif stage_count >= -stage[3]:
      stagedo += ' -if ' + prog[stage_count+stage[3]][4] + '.nc'
   else:
      print 'ERROR: only use previous outfile names if possible';  sys.exit(98)
   ## build rest of options
   stagedo += ' ' + stage[5] + ' ' + always_opts
   if len(stage) > 6:  stagedo += ' -regrid ' + stage[6] + '.nc -regrid_vars ' + stage[7]
   stagedo += ' -o ' + stage[4]
   ## try it with timing
   print '  trying \"' + stagedo + '\"'
   try:
      lasttime = time.time()
      status = os.system(stagedo)
      elapsetime = time.time() - lasttime      
   except KeyboardInterrupt:  sys.exit(2)
   if status:  sys.exit(status)
   print ('  stage %d ' % (stage_count+1)) + ' finished in %7.4f seconds' % elapsetime
   stage_count += 1
   


