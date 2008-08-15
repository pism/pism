#!/usr/bin/env python

# uses NCO to create full length record for ivol from P1cont run;
# then computes exponential time constant

from pylab import *
import os
import sys
from netCDF3 import Dataset as NC

## a successful experiment with polyfit:
#t=arange(0,1,.1)
#v=3.2*exp(-0.8*t)
#vnoise=v+0.1*squeeze(randn(1,10))
#p=polyfit(t,log(vnoise),1)
#plot(t,v,'s-',t,vnoise,'o:',t,exp(p[1]+p[0]*t),'--')

def exponentialFit(t,v,vfinal):
  vshift = v - vfinal
  p = polyfit(t,log(vshift),1)
  plot(t,v,'k',t,exp(p[1] + p[0]*t)+vfinal,'k:')
  show()
  return -1./p[0]

fullname = "P1cont_ivol_full.ser.nc"
runtogether = "ncrcat -O -v ivol P1_10km.ser.nc P1_10km_20k.ser.nc P1_10km_40k.ser.nc P1_10km_60k.ser.nc P1_10km_80k.ser.nc P1_10km_100k.ser.nc %s" % fullname

print "creating %s using NCO (ncrcat)" % fullname
try:
  status = os.system(runtogether)
except KeyboardInterrupt:  sys.exit(2)
if status:  sys.exit(status)

print "opening %s to get ivol time series" % fullname
nc = NC(fullname, 'r')
time = nc.variables["t"][:]
ivol = nc.variables["ivol"][:]
nc.close()

estFinalVol = 2.106  # from direct inspection of .ser.nc
print "graphing exponential fit; close figure to finish"
ivolConst = exponentialFit(time,ivol,estFinalVol)
print "exponential constant for P1cont ivol is %f years" % ivolConst

