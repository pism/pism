#!/usr/bin/env python

'''Check if variable is stationary

'''

try:
    from netCDF3 import Dataset as CDF
except:
    from netCDF4 import Dataset as CDF

from scipy.stats.stats import nanmean
import numpy as np
import pylab as plt
from optparse import OptionParser

__author__ = "Andy Aschwanden"

# default values
THRESHOLD = 1e-5
VAR = 'tempbase'

parser = OptionParser()
parser.usage = "usage: %prog [options] FILE"
parser.description = "Check stationarity of a variable in FILE."
parser.add_option("-t", "--threshold",dest="threshold",
                  help="draws a line horizontal line at THRESHOLD",
                  metavar="THRESHOLD",default=THRESHOLD)
parser.add_option("-v", "--variable",dest="varname",type='string',
                  help="use VAR (default=tempbase)",
                  metavar="VAR",default=VAR)

(options, args) = parser.parse_args()

threshold = options.threshold
varname = options.varname

if len(args) == 1:
    infile = args[0]
else:
    print('wrong number arguments, must be 1')

try:
    nc = CDF(infile,'r')
    try:
        t = nc.variables['t']
    except:
        print('Variable t not found in file %s' % varname)
except IOError:
    pass

        
if varname in nc.variables.keys():
    var = nc.variables[varname]
    dim = var.ndim
else:
    print('error: %s not found in %s' % (varname,infile))
    exit(0)


# differenetiate along time axis
dt = np.diff(t[:],axis=0)

if dim == 1:
    Var = var[:]
    dVar = np.diff(Var)
    dVardt = dVar/dt

elif dim == 3:
    if 'mask' in nc.variables.keys():
        mask     = np.array(nc.variables['mask'][:]) # (t,y,x)
        k        = np.nonzero((mask==1) ^ (mask==2) ^ (mask==3))
        mask2    = np.ones_like(mask)
        mask2[k] = 0
        # use masked values (i.e. ignore ocean and ice-free land)
        Var = np.ma.array(data=var[:],mask=mask2)
    else:
        Var = var[:]

    Varmean = nanmean(Var,axis=1)
    Varmean = nanmean(Varmean,axis=1)
    dVar = np.diff(Varmean)
    dVardt = dVar/dt
    
else:
    print('error: dim n = %i of variable %s not supported, must be 1 or 3' % (dim, varname))


# plot with log-scale y-axis
plt.figure()
plt.hold(True)
plt.semilogy(t[1:],dVardt, 'b', lw = 2)
plt.plot([t[1:][0], t[1:][-1]],[threshold, threshold], 'r',lw = 1)
plt.xlabel(('time [%s]' % t.units))
plt.ylabel(('d %s / dt [%s/a]' % (varname,var.units)))
plt.title(('rate of change of %s' % varname), fontsize=12)
plt.show()

# eventuall, close the file
nc.close()
