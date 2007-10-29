#! /usr/bin/env python
# This script will run the EISMINT-Greenland SSL2 or SSL3 experiment. It runs pgrn
# until there is less than a .01% change in volume in 10,000 years (standard given by
# EISMINT-Greenland)

import sys
import getopt
import time
import commands
import string

# command line arguments
nproc = 1
stdout_file = "out.txt"
IN_FILE = "green20km_Tsteady.nc"  # this input file should contain credible temp field
steadyState = "-ssl2"
ending = "SSL2"
criterion = 0.0001  # 0.01% is default criterion
interval = 10 # 10k model years; option gives multiple of 1000 years
prev_year = 0

try:
  opts, args = getopt.getopt(sys.argv[1:], "i:n:s:t:c:s:",
         ["infile","nproc", "ssl3", "timeinterval", "criterion", "startyeark"])
  for opt, arg in opts:
    if opt in ("-i", "--infile"):
      IN_FILE = arg
    if opt in ("-n", "--nproc"):
      nproc = int(arg)
    if opt in ("-s", "--ssl3"):
      steadyState = "-ssl3 -gk"
      ending = "SSL3"
      print "  [running SSL3 experiment]"
    if opt in ("-t", "--timeinterval"):
      interval = int(arg)
    if opt in ("-c", "--criterion"):
      criterion = float(arg)
    if opt in ("-s", "--startyeark"):
      prev_year = int(arg)
except getopt.GetoptError:
  print 'INCORRECT COMMAND LINE ARGUMENTS; EXITING'
  sys.exit(2)

print "  [using " + str(nproc) + " processors]"
print "  [using interval of " + str(interval) + "k years between volume measurements]"
print "  [using criterion of " + str(criterion*100) + " percent change between volume measurements]"

# set up iteration and test
curr_year = prev_year + interval

# iterate
inname = IN_FILE
count = 0
while True:
  print '  running from ' + str(prev_year) + 'k years until ' + str(curr_year) + 'k years ...'
  outname = 'green_' + ending + '_' + str(curr_year) + 'k'
  cmd = 'mpiexec -n ' + str(nproc) + ' pgrn ' + steadyState + ' -if ' + inname
  cmd += ' -y ' + str(interval) + '000'
  if (count == 0):
    cmd += ' -ys 0'
  cmd += ' -ocean_kill -tempskip 3 -o ' + outname + ' >> ' + stdout_file
  print '  trying:  ' + cmd
  try:
    (status, output)=commands.getstatusoutput(cmd)
  except KeyboardInterrupt:
    sys.exit(2)
  except:
    print 'ERROR OUTPUT: ' + output
    sys.exit(status)

  # Compute the change in "volume".  Actually use the sum of the thicknesses, but since the grid is
  # equally spaced, the result will be the same.  Use the NCO.

  areaPrevName = 'area_' + str(prev_year) + 'k.nc'
  areaCurrName = 'area_' + str(curr_year) + 'k.nc'
  outnamefull = outname + '.nc'
  cmd1='ncwa -O -N -v H -a x,y ' + inname + ' -o ' + areaPrevName
  cmd2='ncwa -O -N -v H -a x,y ' + outnamefull + ' -o ' + areaCurrName
  cmd3='ncks -O -s \'%f\n\' -C -v H ' + areaPrevName + ' | tail -n 1'
  cmd4='ncks -O -s \'%f\n\' -C -v H ' + areaCurrName + ' | tail -n 1'
  try:
    (status, output)=commands.getstatusoutput(cmd1)
    (status, output)=commands.getstatusoutput(cmd2)
  except KeyboardInterrupt:
    sys.exit(2)
  except:
    print 'ERROR OUTPUT: ' + output
    sys.exit(status)
  try:  
    (status, area1)=commands.getstatusoutput(cmd3)
    print '  sum of H for ' + str(prev_year) + 'k years: ' + str(area1)
    (status, area2)=commands.getstatusoutput(cmd4)
    print '  sum of H for ' + str(curr_year) + 'k years: ' + str(area2)
  except KeyboardInterrupt:
    sys.exit(2)
  except:
    print 'ERROR: PROBLEM WITH AREAS'
    sys.exit(status)
  result = abs((float(area1)-float(area2))/float(area2))
  print '  percent difference is: ' + str(result*100) + '%'

  if (result <= criterion):
    print "  resulting volume difference less than criterion; ending"
    break
    
  prev_year = curr_year
  curr_year = curr_year + interval
  inname = 'green_' + ending + '_' + str(prev_year) + 'k.nc'
  count += 1
  #end while True:

