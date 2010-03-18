#!/usr/bin/env python

'''Check if variable is stationary by using a p-norm

'''

try:
    from netCDF3 import Dataset as CDF
except:
    from netCDF4 import Dataset as CDF

import numpy as np
import pylab as plt
from optparse import OptionParser

__author__ = "Andy Aschwanden"

# default values
THRESHOLD = 1e-3
VAR = 'enthalpybase'
PNORM = float(2) # euclidian norm, 2-norm

parser = OptionParser()
parser.usage = "usage: %prog [options] FILE"
parser.description = "Check stationarity of a variable in FILE by calculating the rate of change of its p-norm. That is \
d/dt || X ||_{p} = (\sum_{i}^{m} ( E_{i}^{n+1}-E_{i}^{n} )^p)^(1/p)/(t^{n+1}-t^{n+1}), \
where E_{i}^{n} is the value at time n and coordinate i."
parser.add_option("-p", "--pnorm",dest="pnorm",type='float',
                  help="use P norm (default p = 2)",
                  metavar="P",default=PNORM)
parser.add_option("-t", "--threshold",dest="threshold",
                  help="draws a line horizontal line at THRESHOLD",
                  metavar="THRESHOLD",default=THRESHOLD)
parser.add_option("-v", "--variable",dest="varname",type='string',
                  help="calculate from from variable VAR (default=enthalpybase)",
                  metavar="VAR",default=VAR)

(options, args) = parser.parse_args()

threshold = options.threshold
varname = options.varname
p = options.pnorm

if len(args) == 1:
    infile = args[0]
else:
    print('wrong number arguments, must be 1')
    parser.print_help()
    exit(0)

def load_file(infile):
    
    try:
        nc = CDF(infile,'r')
        try:
            't' in nc.variables.keys()
        except:
            print('Variable t not found in file %s' % infile)
    except IOError:
        pass

    return nc


def getRateOfChange(t,var,p):

    '''
       Calculate rate of change of p-norm of variable
    '''

    # differenetiate time along time axis
    dt = np.diff(t[:],axis=0)

    dim = var.ndim
    if dim == 1:
        Var = var[:]
        dVp = np.diff(Var)
        dVpdt = dVar/dt

    elif dim == 3:
        if 'mask' in nc.variables.keys():
            mask     = np.array(nc.variables['mask'][:]) # (t,y,x)
            k        = np.nonzero((mask==1) ^ (mask==2) ^ (mask==3))
            mask2    = np.ones_like(mask)
            mask2[k] = 0
            # use masked values (i.e. ignore ocean and ice-free land)
            Var = np.ma.array(data=var[:],mask=mask2)
        else:
            Var = np.array(var[:])

        dVar = np.diff(Var,axis=0)
        nt,ny,nx = dVar.shape
        dV = dVar.reshape(nt,nx*ny) # convert (t,y,x) -> (t,y*x) so that np.nansum needs to be applied only once
        dVp = (np.nansum(dV**p,axis=1))**(1./p)
        dVpdt = dVp/dt
    else:
        print('error: dim n = %i of variable %s not supported, must be 1 or 3' % (dim, varname))

    return dVpdt


if __name__ == "__main__":

    '''
       Load file, make sure that time dimension exists, return pointer to file
    '''

    nc = load_file(infile)

    t = nc.variables['t']
    if varname in nc.variables.keys():
        var = nc.variables[varname]
    else:
        print("error: variable '%s' not found in %s" % (varname,infile))
        exit(0)

    dVpdt = getRateOfChange(t,var,p)


    # plot with log-scale y-axis
    plt.figure()
    plt.hold(True)
    plt.semilogy(t[1:],dVpdt, 'b', lw = 2)
    plt.plot([t[1:][0], t[1:][-1]],[threshold, threshold], 'r',lw = 1)
    plt.xlabel(('time [%s]' % t.units))
    plt.ylabel(('d %s / dt [%s/a]' % (varname,var.units)))
    plt.title(('rate of change of %s' % varname), fontsize=12)
    plt.show()

    # eventually, close the file
    nc.close()
