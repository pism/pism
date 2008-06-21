#! /usr/bin/env python
## VERIFYNOW.PY is a script to verify several components of PISM.
## It specifies a refinement path for each of Tests ABCDEFGIJKL and runs pismv accordingly.

## ELB 1/31/07; 2/3/07: -ksp_rtol 1e-6 added; 5/29/07 verifynow.sh --> verifynow.py;
## 7/20/07 used getopt; 10/8/07 mods to work with vnreport.py; 3/9/08 fix to error report
## reading

import sys
import getopt
import time
import commands
import string

# run a chosen verification   
def verify(test):
   print ' ++++ verifying ' + test[2] + ' using **TEST ' + test[0] + '** ++++'
   print ' ' + test[5]
   # for myMx in test[1][:levs]:
   for lev_index in range(levs):
      myMx = test[1][lev_index]
      if test[3] == 0:
     	   gridopts = ' -Mx ' + str(myMx) + ' -My ' + str(myMx)
      elif test[3] == 1:
         gridopts = ' -My ' + str(myMx)
      elif test[3] == 2:
         gridopts = ' -Mx ' + str(myMx) + ' -My ' + str(myMx)
         if (uneq == 0):
           gridopts = gridopts + ' -Mz ' + str(myMx)
         else:
           gridopts = gridopts + ' -Mz ' + str(test[6][lev_index])
      elif test[3] == 3:
         myMz = myMx
         myMbz = (myMz - 1) / 4 + 1
         mymaxdt = 40000.0 / (float((myMbz - 1) * (myMbz - 1)))
         mymaxdt = 0.5 * int(2.0 * mymaxdt)  # floor to nearest half-year
         if (uneq != 0):
           myMz = test[6][lev_index]
         gridopts = ' -Mz ' + str(myMz) + ' -Mbz ' + str(myMbz) + ' -max_dt ' + str(mymaxdt)
      if nproc > 1:
         predo = mpi + ' -np ' + str(nproc) + ' '
      else:
         predo = ''
      testdo = predo + pref + 'pismv -test ' + test[0] + gridopts + test[4]
      if (uneq == 1):
        testdo = testdo + ' -quadZ'
      elif (uneq == 2):
        testdo = testdo + ' -chebZ'
      print ' trying \"' + testdo + '\"'
      testdo = testdo + ' -verbose 1'  # only need final errors anyway
      try:
         lasttime = time.time()
         (status,output) = commands.getstatusoutput(testdo)
         elapsetime = time.time() - lasttime      
      except KeyboardInterrupt:
         sys.exit(2)
      if status:
         sys.exit(status)
      print ' finished in %7.4f seconds; reported numerical errors as follows:' % elapsetime
      errpos = output.find('NUMERICAL ERRORS')
      if errpos >= 0:
         # errreport = output[errpos:output.rfind('Writing')-1]
         # print '  |' + string.replace(errreport,'\n','\n  |')
         errreport = output[errpos:output.find('Writing')]
         endline = errreport.find('\n')
         print '    ' + errreport[0:endline]
         errreport = errreport[endline+1:]
         while (len(errreport) > 1) and (endline > 0):
           endline = errreport.find('\n')
           if endline == -1:
             endline = len(errreport)
           print '   #' + errreport[0:endline]
           errreport = errreport[endline+1:]       
           endline = errreport.find('\n')
           if endline == -1:
             endline = len(errreport)
           print '   |' + errreport[0:endline]
           errreport = errreport[endline+1:]       
      else:
         print ' ERROR: can\'t find reported numerical error'
         sys.exit(99)

## default settings
NP = 1
LEVELS = 2
PREFIX = ''
MPIDO = 'mpiexec'
TESTS = 'CGIJ'
KSPRTOL = 1e-12 # for test I
SSARTOL = 5e-7   # ditto

## tests and additional info for verification
## order here is for convenience and speed: generally do easier (faster) tests first
alltests = [
   ['B',[31,41,61,81,121],
        'isothermal SIA w moving margin',0,' -Mz 31 -no_eta -ys 422.45 -y 25000.0',
        '(refine dx=80,60,40,30,20,km, dx=dy and Mx=My=31,41,61,81,121)'],
   ['C',[41,61,81,101,121],
        'isothermal SIA w moving margin',0,' -Mz 31 -no_eta -y 15208.0',
        '(refine dx=50,33.33,25,20,16,km, dx=dy and Mx=My=41,61,81,101,121)'],
   ['D',[41,61,81,101,121],
        'isothermal SIA w variable accumulation',0,' -Mz 31 -no_eta -y 25000.0',
        '(refine dx=50,33.33,25,20,16.67,km, dx=dy and Mx=My=41,61,81,101,121)'],
   ['A',[31,41,61,81,121],
        'isothermal SIA w marine margin',0,' -Mz 31 -no_eta -y 25000.0',
        '(refine dx=53.33,40,26.67,20,13.33,km, dx=dy and Mx=My=31,41,61,81,121)'],
   ['E',[31,41,61,81,121],
        'isothermal SIA w sliding',0,' -Mz 31 -no_eta -y 25000.0',
        '(refine dx=53.33,40,26.67,20,13.33,km, dx=dy and Mx=My=31,41,61,81,121)'],
   ['H',[31,41,61,81,121],
        'isothermal SIA w moving margin and isostatic bed deformation\n     ',
        0,' -Mz 31 -no_eta -bed_def_iso -y 60000.0',
        '(refine dx=80,60,40,30,20,km, dx=dy and Mx=My=31,41,61,81,121)'],
   ['I',[49,193,769,3073,12289],'plastic till ice stream',1,
        ' -Mx 5 -ssa_rtol ' + str(SSARTOL) + ' -ksp_rtol ' + str(KSPRTOL),
        '(refine dy=5000,1250,312.5,78.13,19.53,m, My=49,193,769,3073,12289)'],
   ['J',[30,60,120,180,240],'linearized, periodic ice shelf',0,
        ' -Mz 11 -ksp_rtol ' + str(KSPRTOL),
        '(refine dx=20,10,5,3.333,2.5, km; dx=dy and My=30,60,120,180,240)'],
   ['K',[41,81,161,321,641],
        'pure conduction problem in ice and bedrock',3,' -Mx 6 -My 6 -y 130000.0',
        '(refine dz=100,50,25,12.5,6.25,m, Mz=41,81,161,321,641)',
        [15,28,55,108,215]],
   ['L',[31,61,91,121,181],
        'isothermal SIA w non-flat bed',0,' -Mz 31 -no_eta -y 25000.0',
        '(refine dx=60,30,20,15,10,km, dx=dy and Mx=My=31,61,91,121,181)'],
   ['F',[61,91,121,181,241],'thermocoupled SIA',2,' -y 25000.0',
        '(refine dx=30,20,15,10,7.5,km, dx=dy, dz=66.67,44.44,33.33,22.22,16.67 m\n'
        + '  and Mx=My=Mz=61,91,121,181,241)',
        [21,31,41,61,81]],
   ['G',[61,91,121,181,241],'thermocoupled SIA w variable accum',2,' -y 25000.0',
        '(refine dx=30,20,15,10,7.5,km, dx=dy, dz=66.67,44.44,33.33,22.22,16.67 m\n'
        + '  and Mx=My=Mz=61,91,121,181,241)',
        [21,31,41,61,81]],
]

## get options:
##   -l for number of levels
##   -m for MPI executable (e.g. mpirun or mpiexec)
##   -n for number of processors
##   -p for prefix on pismv executable
##   -t for which tests to use
##   -u to add unequal spaced vertical (-u 0 [default] for equal spaced
##                        -u 1 for "-quadZ", -u 2 for "-chebZ")
nproc = NP  ## default; will not use 'mpiexec' if equal to one
levs = LEVELS
mpi = MPIDO
pref = PREFIX
letters = TESTS
uneq = 0
try:
  opts, args = getopt.getopt(sys.argv[1:], "p:m:n:l:t:u:",
                             ["prefix=", "mpido=", "nproc=", "levels=", "tests=", "unequal="])
  for opt, arg in opts:
    if opt in ("-p", "--prefix"):
      pref = arg
    elif opt in ("-m", "--mpido"):
      mpi = arg
    elif opt in ("-n", "--nproc"):
      nproc = string.atoi(arg)
    elif opt in ("-l", "--levels"):
      levs = string.atoi(arg)
    elif opt in ("-t", "--tests"):
      letters = arg
    elif opt in ("-u", "--uneq"):
      uneq = string.atoi(arg)
except getopt.GetoptError:
  print 'Incorrect command line arguments'
  sys.exit(2)
# should do range checking on option arguments here

print ' VERIFYNOW.PY using %d processor(s) and %d level(s) of refinement' % (nproc, levs)
## go through verification tests
for test in alltests:
   if letters.find(test[0]) > -1:
       verify(test)

