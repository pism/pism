#!/usr/bin/env python

# Produces a PNG figure based on final steady-state MISMIP output.
# For example, if you have just run  examples/mismip/mismip.sh
# and gotten files 'ABC1_1a_M1_A1_ss' and 'ABC1_1a_M1_A1_f'
# (among others), then do
#   ./figsMISMIP.py -p ABC1_1a_M1_A1

from pylab import *
import os
import sys
from getopt import getopt, GetoptError

MAXTHICK = 4500.0  # so that all plotting has save vertical scale

haveExtras = False
mprefix='ABC1_1a_M1_A1'
outputname='foo.png'
try:
  opts, args = getopt(sys.argv[1:], "e:o:p:",["extras=:outfile=:prefix="])
  for opt, arg in opts:
    if opt in ("-e", "--extras"):
      haveExtras = True
    if opt in ("-p", "--prefix"):
      mprefix = arg
    if opt in ("-o", "--outfile"):
      outputname = arg
    else:
      outputname = 'profile_' + mprefix + '.png'
except GetoptError:
  print "Incorrect command line arguments. Exiting..."
  sys.exit(-1)
  
try:
  name = mprefix + '_ss'
  A = load(name)
except IOError:
  print "file '%s' not found" % name
  sys.exit(2)

try:
  name = mprefix + '_f'
  B = load(name)
except IOError:
  print "file '%s' not found" % name
  sys.exit(2)

figure(1)

if haveExtras:
  try:
    name = mprefix + '_extras'
    E = load(name)
  except IOError:
    print "can't find _extras, so showing just thickness from _ss"
    haveExtras = False

if haveExtras:
  plot(A[:,0],E[:,0]-A[:,1],'c.-',markersize=9)
  hold(True)
  plot(A[:,0],zeros(size(A[:,0])),'r:',linewidth=2)
  plot(A[:,0],E[:,0],'b.-',markersize=9)
  plot(A[:,0],E[:,1],'k.-',markersize=9)
  plot(array([B[0]]),array([-1000.0]),'kd',markersize=14)
  hold(False)
  title('profile ($t_f$ = %.2f a, $x_g$ = %.3f km)' % (B[1],B[0]))
  ylabel(r'elevation (m)',size=14)
  myax = axis()
  axis((myax[0],myax[1],-1000.0,MAXTHICK+700.0))
else:
  plot(A[:,0],A[:,1],'b.-',markersize=12)
  hold(True)
  plot(array([B[0]]),array([0.0]),'kd',markersize=16)
  hold(False)
  title('thickness ($t_f$ = %.2f a, $x_g$ = %.3f km)' % (B[1],B[0]))
  ylabel(r'$H(x,t_f)$ (m)',size=14)
  myax = axis()
  axis((myax[0],myax[1],0.0,MAXTHICK))

text(1200,4500,mprefix,size=14)
xlabel(r'$y$ (km)',size=14)

print "  saving figure %s" % outputname
savefig(outputname, dpi=300, facecolor='w', edgecolor='w')
  

