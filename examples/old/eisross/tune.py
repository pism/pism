#! /usr/bin/env python

## ELB 8/6/07; 11/17/07

# FIXME (CK, Dec 12, 2011): This script is out of date. 'ice_type custom' is
# gone and so is the -ice_custom_hardness option. We need to use
# -config_override for parameter studies like this one.

import sys
import getopt
import time
import commands
import string

## default settings
KSPRTOL = 1e-8
MVRTOL = 1e-5

pref = ''
nproc = 1
try:
  opts, args = getopt.getopt(sys.argv[1:], "p:n:",
                             ["prefix=", "nproc="])
  for opt, arg in opts:
    if opt in ("-p", "--prefix"):
      pref = arg
    elif opt in ("-n", "--nproc"):
      nproc = string.atoi(arg)
except getopt.GetoptError:
  print 'Incorrect command line arguments'
  sys.exit(2)

print 'TUNE.PY (for EISMINT Ross; compare Chi^2 results to table 1 in MacAyeal et al 1996)'
for hard in [1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2]:
   rossdo = ''
   if nproc > 1:
     rossdo += 'mpiexec -n ' + str(nproc) + ' '
   rossdo += 'pross -boot_file ross.nc -ssaBC ross.nc -riggs riggs.nc'
   rossdo += ' -ksp_rtol ' + str(KSPRTOL) + ' -ssa_rtol ' + str(MVRTOL) 
   rossdo += ' -Mx 147 -My 147 -Lz 1000 -Mz 3 -ssa -flow_law isothermal_glen -ice_custom_hardness ' + str(hard) + 'e8'
   print 'trying \"' + rossdo + '\"'
   try:
      lasttime = time.time()
      (status,output) = commands.getstatusoutput(rossdo)
      elapsetime = time.time() - lasttime      
   except KeyboardInterrupt:
      sys.exit(2)
   if status:
      sys.exit(status)
   print '  finished in %7.4f seconds; max computed speed and Chi^2 as follows:' % elapsetime
   maxvelpos = output.find('maximum computed speed')
   if maxvelpos:
      report = output[maxvelpos:output.rfind('ERRORS')-1]
      print '  |' + report
   else:
      print ' ERROR: can\'t find maximum computed speed in output'
   chisqrpos = output.find('Chi^2 statistic')
   if chisqrpos:
      report = output[chisqrpos:output.rfind('is')+11]
      print '  |' + report
   else:
      print ' ERROR: can\'t find Chi^2 statistic in output'

