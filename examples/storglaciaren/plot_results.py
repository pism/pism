#!/usr/bin/env python

try:
    from netCDF3 import Dataset
except:
    from netCDF4 import Dataset
import numpy as np
import pylab as plt

from optparse import OptionParser

parser = OptionParser()
parser.usage = "usage: %prog [options] FILE"
parser.description = "A script to compare PISM flowline velocities with full Stokes solution."

(options, args) = parser.parse_args()

if len(args) != 1:
    print('wrong number of arguments, 1 expected')
    exit(1)
    
try:
    nc = Dataset(args[0], 'r')
except:
    print(("file %s not found ... ending ..." % args[0]))
    exit(2)

def permute(variable, output_order = ('t', 'z', 'zb', 'y', 'x')):
    """Permute dimensions of a NetCDF variable to match the output storage order."""
    input_dimensions = variable.dimensions

    # filter out irrelevant dimensions
    dimensions = filter(lambda(x): x in input_dimensions,
                        output_order)

    # create the mapping
    mapping = map(lambda(x): dimensions.index(x),
                  input_dimensions)

    if mapping:
        return np.transpose(variable[:], mapping)
    else:
        return variable[:]              # so that it does not break processing "mapping"


x = nc.variables["x"][:]
b = np.squeeze(nc.variables["topg"][:])
s = np.squeeze(nc.variables["usurf"][:])
h = np.squeeze(nc.variables["thk"][:])
z = nc.variables["z"][:]

mask = np.zeros_like(h)
mask[h<=1] = 1

us = np.ma.array(data=np.squeeze(nc.variables["uvelsurf"][:]),mask=mask)
ub = np.ma.array(data=np.squeeze(nc.variables["uvelbase"][:]),mask=mask)
## stuff needed for contour plots
xx = np.tile(x,[len(z),1]) 
zz = (np.tile(z,[len(x),1])).transpose() + b

cts = np.squeeze(permute(nc.variables["cts"]))
liqfrac = np.squeeze(permute(nc.variables["liqfrac"]))
temppa = np.squeeze(permute(nc.variables["temp_pa"]))
mask2 = np.zeros_like(cts)
mask2[zz>s] = 1
cts = np.ma.array(data=cts,mask=mask2)
liqfrac = np.ma.array(data=liqfrac,mask=mask2)
temppa = np.ma.array(data=temppa,mask=mask2)


fig = plt.figure(figsize=(6.4,7.4))

axUpper = plt.axes([0.1,0.6,0.8,0.25])
axLower = plt.axes([0.1,0.05,0.8,0.5])

axUpper.plot(x,us,color='#377EB8', lw = 1.5)
axUpper.plot(x,ub,'--',color='#377EB8', lw = 1.5)
axUpper.axis([-250, 3500, 0, 50])
plt.ylabel("velocity [m a$^{-1}$]")
plt.setp(axUpper, xticks=[])
## Contour level of the CTS
cts_level = [1,1]
liqfrac_levels = np.arange(0,1.5,.25)
#temppa_levels = np.arange(-6,0,1)
temppa_levels = [-6,-5,-4,-3,-2,-1,-.0001]

axLower.plot(x,b,color='black', lw = 1.5)
axLower.plot(x,s,color='black', lw = 1.5)
c1=axLower.pcolor(xx,zz,liqfrac*100,cmap=plt.cm.Reds,vmin=0.5)
plt.colorbar(mappable=c1,ax=axLower,orientation='horizontal',pad=0.05,shrink=0.75)
#c2=axLower.contourf(xx,zz,temppa,temppa_levels,cmap=plt.cm.Blues_r,lw = 1)
c2=axLower.pcolor(xx,zz,temppa,cmap=plt.cm.Blues_r,vmin=-6,vmax=0)
c3=axLower.pcolor(xx,zz,liqfrac*100,cmap=plt.cm.Reds,vmin=0.01)
plt.colorbar(mappable=c2,ax=axLower,orientation='horizontal',ticks=[-6,-5,-4,-3,-2,-1,0],pad=0.20,shrink=0.75)
axLower.contour(xx,zz,cts,cts_level,colors='black',linestyles='dashed',lw = 1)
axLower.axis([-250, 3500, 1100, 1800])
plt.xlabel("distance from bergschrund [m]")
plt.ylabel("elevation [m a.s.l.]")


plt.savefig('sg_results.pdf',bbox_inches='tight',pad_inches=0.35)

nc.close()
