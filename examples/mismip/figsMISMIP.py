#!/usr/bin/env python

# produces figures based on final (steady-state) MISMIP output

from pylab import *

mprefix='EBU1_1a_M1_A1'

A=load(mprefix + '_ss')
B=load(mprefix + '_f')
figure(1)
plot(A[:,0],A[:,1],'b.-',array([B[0]]),array([0.0]),'kd')
title('profile at steady state ($t_f$ = %.2f a)' % B[1])
xlabel(r'$x$')
ylabel(r'$H(x,t_f)$')

figfilename='steady_' + mprefix + '.png'
print "  saving figure %s" % figfilename
savefig(figfilename, dpi=300, facecolor='w', edgecolor='w')

