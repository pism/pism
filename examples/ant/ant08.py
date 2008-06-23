#! /usr/bin/env python

## ant08.py is designed for modelling the current state of Antarctica on a
## multi-processor machine. Ice streams are determined by using balance
## velocities OR Schoof's plastic till model. Assumes antarctica_init.nc is
## present (which contains 5km data mostly from the British Antarctica
## Service.) Tends to require many processors just to have memory for model
## state. Contact ffelb@uaf.edu with questions. ELB 1/31/07; 5/29/07: .sh -->
## .py; 16june08: total revision

import os
import sys
import time
from getopt import getopt, GetoptError

## default settings
NP = 8
BINITFILE = 'antarctica_init.nc'

## program prog is a list; elements of prog are:
##     ['DESCRIPTION','EXECUTABLE',YEARS,INITFLAG,'OUTNAME',
##                'OTHER_OPTIONS'[,'REGRIDFILENAME','REGRIDVARS']]
## note: INITFLAG == 0   <==>   '-bif BINITFILE'
##       INITFLAG == -1   <==>   '-if [prev OUTNAME]'          {case of next}
##       INITFLAG < 0    <==>   '-if [OUTNAME from stage[current-|INITFLAG|]]'
prog = [
['input data, put on 40km grid, and smooth surface',
      'pismr',    10,0,'ant10yr_40km',
      '-Mx 141 -My 141 -Mz 51 -Mbz 11 -Lz 5000 -quadZ -ocean_kill'],
['running for 200k year to equilibriate temp (steady geometry; SIA only)',
      'pismr',200000,-1,'ant200k_40km','-no_mass'],
['regrid to 20km grid: interpolate temp and age and smooth surface',
      'pismr',    10,0,'ant200k_20km',
      '-Mx 281 -My 281 -Mz 101 -Mbz 21 -Lz 5000 -quadZ -ocean_kill',
      'ant200k_40km','BeLT'],
['running for another 10000 yrs to equilibriate temp on fine grid (steady geometry and SIA only) ',
      'pismr',  10000,-1,'ant210k_20km','-no_mass'],
['2 yrs with SSA',
      'pismr',     2,-1,'ant2yr_ssa_20km','-ssa -super -plastic -ocean_kill -ys 0 -verbose 5'],
['another 500 yrs with SSA',
      'pismr',   500,-1,'ant502yr_ssa_20km','-ssa -super -plastic -ocean_kill']]
#['regrid to 10km grid and run for 100 years',
#      'pismr',   100,0,'ant100yr_ssa_10km',
#      '-Mx 561 -My 561 -Mz 101 -Mbz 21 -Lz 5000 -quadZ -ssa -super -plastic -ocean_kill -ys 0',
#      'ant502yr_ssa_20km','BeLT']]

#always_opts = '-gk -e 1.2'
always_opts = '-e 3'

## get options: -n for number of processors, -s for the stage to start from
## (numbered from 1), -d to debug (with this option ant08.py just prints
## commands it would run in the normal mode).
try:
   opts, args = getopt(sys.argv[1:], "n:s:d")
   nproc = NP                 # default; will not use 'mpiexec' if equal to one
   start_from = 0             # start from the first stage by default
   debug = False
   for opt, arg in opts:
      if opt == "-n":
         nproc = int(arg)
      if opt == "-s":
         start_from = int(arg) - 1
      if opt == "-d":
         debug = True
except GetoptError, message:
   print "Getopt error:", message

if start_from >= len(prog):
   print "Error: requested to start from stage %d out of %d" % (start_from,
                                                                len(prog))
   print "Starting from stage #1."

sys.stdout.flush()

## run program: stage[] above determines stages of computation
print '  ANT08 using %d processor(s)' % nproc
if start_from > 0:
   print "  Requested to start from stage %d." % (start_from + 1)
for j in range(start_from, len(prog)):
   stage = prog[j]
   print ('  stage %d.  ' % (j + 1)) + stage[0] + ':'

   ## use mpi only if nproc > 1
   stagedo = ''
   if nproc > 1:  stagedo += 'mpiexec -n %d' % nproc
   stagedo += ' ' + stage[1] + (' -y %d' % stage[2]) 
   ## determine infile
   if stage[3] == 0:
      stagedo += ' -bif ' + BINITFILE
   elif j >= -stage[3]:
      stagedo += ' -if ' + prog[j + stage[3]][4] + '.nc'
   else:
      print 'ERROR: only use previous outfile names if possible';  sys.exit(98)
   ## build rest of options
   stagedo += ' ' + stage[5]
   if len(stage) > 6:  stagedo += ' -regrid ' + stage[6] + '.nc -regrid_vars ' + stage[7]
   stagedo += ' -o ' + stage[4] + ' ' + always_opts

   ## try it with timing
   if debug:
      print "  Would try\n%s" % stagedo
      continue
   print '  trying \"' + stagedo + '\"'
   sys.stdout.flush()
   try:
      lasttime = time.time()
      status = os.system(stagedo)
      elapsetime = time.time() - lasttime      
   except KeyboardInterrupt:  sys.exit(2)
   if status:  sys.exit(status)
   print ('  stage %d ' % (j + 1)) + ' finished in %7.4f seconds' % elapsetime
   

