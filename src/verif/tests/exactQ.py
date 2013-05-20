#!/usr/bin/env python

# Computes and plots exact similarity solution "test Q".  See equation
# (3.19) in [\ref PeglerListerWorster2012] = PLW2012:
#   S. Pegler, J. Lister, and M. G. Worster, 2012.  Release of a viscous
#   power-law fluid over an inviscid ocean", J. Fluid Mech. 700, 63--76.

from pylab import *

SperA= 31556926.0

g    = 9.81
rho  = 910.0    # density of ice; kg/m^3
rhow = 1028.0   # density of ocean water; kg/m^3

barB = 1.9e8    # strength of shelf; Pa s^(1/3); from MacAyeal et al 1996;
                # is this equal to \tilde mu or not?

n = 3.0
m = (1.0/n) - 1.0

gprime = (rhow - rho) * g / rhow  # see just after (2.1) in PLW2012

nurescale = 3.0**(m/2) * barB / rho # see just after (3.19) in PLW2012

# FIXME: just guessing 

V = 1.0e14   # try this in m^3 = 10^5 km^3

t = 10.0 * SperA

tnt = 2.0 * n * t

H0 = (12.0 * nurescale / gprime) * tnt**(-1.0/n)

cR = (V * gprime / (12.0 * pi * nurescale))**(1.0/2.0)
R0 = cR * tnt**(1.0/(2.0*n))

print 'volume    V  = %f m^3 = %f km^3' % (V, V/1.0e9)
print 'thickness H0 = %f m' % H0
print 'radius    R0 = %f m = %f km' % (R0, R0/1.0e3)

