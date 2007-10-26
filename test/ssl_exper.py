#! /usr/bin/env python
# This script will run the EISMINT-Greenland SSL2 or SSL3 experiment. It runs pgrn
# until there is less than a .01% change in volume in 10,000 years

import sys
import getopt
import time
import commands
import string

# command line arguments
num_proc=1
stdout_file = "out.txt"
IN_FILE = "green20km_Tsteady.nc"
steadyState = "-ssl2"
out_end = "SSL2"

try:
  opts, args = getopt.getopt(sys.argv[1:], "i:n:s", ["infile","num_proc", "ssl3"])
  for opt, arg in opts:
    if opt in ("-i", "--infile"):
      IN_FILE = arg
    if opt in ("-n", "--num_proc"):
      num_proc = int(arg)
    if opt in ("-s", "--ssl3"):
      steadyState = "-ssl3 -gk"
      out_end = "SSL3"
      print "Running SSL3 experiment"
except getopt.GetoptError:
  print 'Incorrect command line arguments'
  sys.exit(2)

print 'Running with ' + str(num_proc) + ' processors'

# run
try:
  print '1. run for 100k years to head toward geometry and thermocoupled steady state'
  cmd = 'mpiexec -n '+str(num_proc)+' pgrn ' + steadyState + ' -if ' + IN_FILE 
  cmd += ' -y 1e5 -ys 0 -ocean_kill -tempskip 3 -o green20km_' + out_end + '100k >> ' + stdout_file
  print cmd
  (status, output)=commands.getstatusoutput(cmd)

except KeyboardInterrupt:
  sys.exit(2)
except:
  print 'output: '+output
  sys.exit(status)

# set up the final section
result = 1
curr_year=110
prev_year=100

# build command
while result > .0001:
  #run pism
  print 'running until '+str(curr_year)+'k years'
  run_pism = 'mpiexec -n ' + str(num_proc) + ' pgrn ' + steadyState
  run_pism += ' -if green20km_' + out_end + str(prev_year) + 'k.nc -y 1e4 -ocean_kill -tempskip 3 '
  run_pism += ' -o green20km_' + out_end + str(curr_year) + 'k >> ' + stdout_file
  print run_pism
  try:
    (status, output)=commands.getstatusoutput(run_pism)
  except KeyboardInterrupt:
    sys.exit(2)

  # compute the change in "volume" it's actually the
  # sum of the thicknesses, but since the grid is
  # equally spaced, the result will be the same

  cmd1='ncwa -O -N -v H -a x,y green20km_'+out_end+str(prev_year)+'k.nc -o area_'+str(prev_year)+'k.nc'
  cmd2='ncwa -O -N -v H -a x,y green20km_'+out_end+str(curr_year)+'k.nc -o area_'+str(curr_year)+'k.nc'
  cmd3='ncks -O -s \'%f\n\' -C -v H area_'+str(prev_year)+'k.nc | tail -n 1'
  cmd4='ncks -O -s \'%f\n\' -C -v H area_'+str(curr_year)+'k.nc | tail -n 1'
  try:
    (status, output)=commands.getstatusoutput(cmd1)
    (status, output)=commands.getstatusoutput(cmd2)
  except KeyboardInterrupt:
    sys.exit(2)
  except:
    print output
    sys.exit(status)
  try:  
    (status, area1)=commands.getstatusoutput(cmd3)
    print 'Sum of H for '+str(prev_year)+'k: '+str(area1)
    (status, area2)=commands.getstatusoutput(cmd4)
    print 'Sum of H for '+str(curr_year)+'k: '+str(area2)
  except KeyboardInterrupt:
    sys.exit(2)
  except:
    print 'problem with areas'
    sys.exit(status)
  result=abs((float(area1)-float(area2))/float(area2))
  print 'percent difference is: '+str(result*100)+'%'
  prev_year=curr_year
  curr_year=curr_year+10

