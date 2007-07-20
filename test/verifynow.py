#! /usr/bin/env python

## VERIFYNOW.PY is a script to verify several components of PISM.  Uses tests C, I and G.  
## It is intended to do roughly the minimal amount of computation to show convergence to continuum results.

## ELB 1/31/07; 2/3/07: -ksp_rtol 1e-6 added; 5/29/07 verifynow.sh --> verifynow.py; 7/20/07 used getopt

import sys
import getopt
import time
import commands
import string

## default settings
NP = 1
LEVELS = 3
KSPRTOL = 1e-12
MVRTOL = 5e-7
PREFIX = ''
MPIDO = 'mpiexec'

## tests and additional info for verification
alltests = [
   ['C',[41,61,81,101,121],
        'isothermal SIA',0,' -Mz 31 -y 15208.0',
        ' (Mx=My=41,61,81,101,121\n       corresponds to dx=dy=50,33.3,25,20,16 km)'],
   ['I',[49,193,769,3073,12289],'plastic till ice stream',1,
        ' -Mx 5 -mv_rtol ' + str(MVRTOL) + ' -ksp_rtol ' + str(KSPRTOL),
        ' (My=49,193,769,3073,12289\n       corresponds to dy=5000,1250,312.5,78.13,19.53 m)'],
   ['G',[61,91,121,181,241],'thermocoupled SIA',2,' -y 25000.0',
        ' (Mx=My=Mz=61,91,121,181,241\n       corresponds to dx=dy=30,20,15,10,7.5 km and dz=66.7,44.4,33.3,22.2,16.7 m)']
]

## get options: -n for number of processors, -l for number of levels
nproc = NP  ## default; will not use 'mpiexec' if equal to one
levs = LEVELS
mpi = MPIDO
pref = PREFIX
try:
  opts, args = getopt.getopt(sys.argv[1:], "p:m:n:l:", ["prefix=", "mpido=", "nproc=", "levels="])
  for opt, arg in opts:
    if opt in ("-p", "--prefix"):
      pref = arg
    elif opt in ("-m", "--mpido"):
      mpi = arg
    elif opt in ("-n", "--nproc"):
      nproc = string.atoi(arg)
    elif opt in ("-l", "--levels"):
      levs = string.atoi(arg)
except getopt.GetoptError:
  print 'Incorrect command line arguments'
  sys.exit(2)

print '  VERIFYNOW using %d processor(s) and %d level(s) of refinement' % (nproc, levs)
## go through verification tests
for test in alltests:
   print '  ++++ verifying ' + test[2] + ' using test ' + test[0] + test[5] + ' ++++'
   for myMx in test[1][:levs]:
      if test[3] == 0:
     	   gridopts = ' -Mx ' + str(myMx) + ' -My ' + str(myMx)
      elif test[3] == 1:
         gridopts = ' -My ' + str(myMx)
      elif test[3] == 2:
         gridopts = ' -Mx ' + str(myMx) + ' -My ' + str(myMx) + ' -Mz ' + str(myMx)
      if nproc > 1:
         predo = mpi + ' -np ' + str(nproc) + ' '
      else:
         predo = ''
      testdo = predo + pref + 'pismv -test ' + test[0] + gridopts + test[4] + ' -verbose 1'
      print '  trying \"' + testdo + '\"'
      try:
         lasttime = time.time()
         (status,output) = commands.getstatusoutput(testdo)
         elapsetime = time.time() - lasttime      
      except KeyboardInterrupt:
         sys.exit(2)
      if status:
         sys.exit(status)
      print '  finished in %7.4f seconds; reported numerical errors as follows:' % elapsetime
      errpos = output.find('NUMERICAL ERRORS')
      if errpos:
         errreport = output[errpos:output.rfind('Writing')-1]
         print '    |' + string.replace(errreport,'\n','\n    |')
      else:
         print '  ERROR: can\'t find reported numerical error'
         sys.exit(99)


