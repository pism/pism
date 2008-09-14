#!/usr/bin/env python

# produces figures based on final (steady-state) MISMIP output

from pylab import *
import os
import sys
from getopt import getopt, GetoptError

mprefix='EBU1_1a_M1_A1'
try:
  opts, args = getopt(sys.argv[1:], "p:",["prefix="])
  for opt, arg in opts:
    if opt in ("-p", "--prefix"):
      mprefix = arg
except GetoptError:
  print "Incorrect command line arguments. Exiting..."
  sys.exit(-1)
  
try:
  A=load(mprefix + '_ss')
except IOError:
  print "file '%s' not found" % mprefix + '_ss'
  sys.exit(2)

B=load(mprefix + '_f')
figure(1)
plot(A[:,0],A[:,1],'b.-',array([B[0]]),array([0.0]),'kd')
title('profile at steady state ($t_f$ = %.2f a; $x_g$ = %.3f km.)' % (B[1],B[0]))
xlabel(r'$x$')
ylabel(r'$H(x,t_f)$')

figfilename='steady_' + mprefix + '.png'
print "  saving figure %s" % figfilename
savefig(figfilename, dpi=300, facecolor='w', edgecolor='w')

