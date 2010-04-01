#!/usr/bin/env python

## @package check_stationarity
# \author Andy Aschwanden, University of Alaska Fairbanks, USA
# \brief Script to evaluate stationarity of a variable.
# \details Given a time series of a variabale \f$ X \f$, it computes the time rate of change of the p-norm of \f$ X \f$.
# \f[\frac{d}{dt} || X ||_{p} = \frac{ \left [\sum_{i}^{m} \left( X_{i}^{n+1} - X_{i}^{n} \right)^{p} \right] ^{1/p} }{ t^{n+1} - t^{n+1} } \f] where \f$ X_{i}^{n} \f$ is the value at time \f$ n \f$ and coordinate \f$ i \f$.


try:
    from netCDF3 import Dataset as CDF
except:
    from netCDF4 import Dataset as CDF

import numpy as np
import pylab as plt
from optparse import OptionParser

__author__ = "Andy Aschwanden"

## If no threshold is given, set it to 1e1
THRESHOLD = 1e1
## This is the default variable to calculate norm from
X = 'enthalpybase'
## Default norm is the euclidian norm, 2-norm
PNORM = float(2) 

## Set up the option parser
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
                  help="calculate from from variable X (default=enthalpybase)",
                  metavar="VAR",default=X)
## Run the option parser
(options, args) = parser.parse_args()

## Assign threshold
threshold = options.threshold
## Assign variable name
varname = options.varname
## Assign norm
p = options.pnorm

if len(args) == 1:
    ## Retrieve command line argument (name of input file)
    infile = args[0]
else:
    print('wrong number arguments, must be 1')
    ## If number of arguments is wrong, print help to give user some clues.
    parser.print_help()
    exit(0)


## Opens a netCDF file.
#
# Open netCDF file and check that time dimension exists. On success, a pointer to the file is returned,
# otherwise an error is issued.
def open_ncfile(infile):

    try:
        nc = CDF(infile,'r')
        try:
            't' in nc.variables.keys()
        except:
            print('Variable t not found in file %s' % infile)
    except IOError:
        pass

    return nc

## Calculate rate of change.
#
# Calculate rate of change of p-norm of variable \f$ X \f$.
def getRateOfChange(t,X,p,varname):

    # differenetiate time along time axis
    dt = np.diff(t[:],axis=0)

    Xdim = X.ndim
    Xunits = X.units
    
    if Xdim == 1:
        X = X[:]
        dXp = np.diff(X)
        dXpdt = dX/dt
    elif Xdim == 3:
        if 'mask' in nc.variables.keys():
            mask     = np.array(nc.variables['mask'][:]) # (t,y,x)
            k        = np.nonzero((mask==1) ^ (mask==2) ^ (mask==3))
            mask2    = np.ones_like(mask)
            mask2[k] = 0
            # use masked values (i.e. ignore ocean and ice-free land)
            X = np.ma.array(data=X[:],mask=mask2)
        else:
            X = np.array(X[:])

        dX = np.diff(X,axis=0)
        nt,ny,nx = dX.shape
        dX = dX.reshape(nt,nx*ny) # convert (t,y,x) -> (t,y*x) so that np.nansum needs to be applied only once
        dXp = (np.nansum(dX**p,axis=1))**(1./p)
        dXpdt = dXp/dt
    else:
        print('error: dim n = %i of variable %s not supported, must be 1 or 3' % (Xdim, varname))

    return dXpdt


if __name__ == "__main__":

    ## Open netCDF file
    nc = open_ncfile(infile)

    ## time variable t
    t = nc.variables['t']
    if varname in nc.variables.keys():
        var = nc.variables[varname]
    else:
        print("error: variable '%s' not found in %s" % (varname,infile))
        exit(0)

    ## Calculate rate of change from time t, variable var, norm p and variable name varname 
    dVpdt = getRateOfChange(t,var,p,varname)


    ## Make plot with log-scale y-axis
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
