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


x = nc.variables["x"][:]
b = np.squeeze(nc.variables["topg"][:])
s = np.squeeze(nc.variables["usurf"][:])
z = nc.variables["z"][:]
us = np.squeeze(nc.variables["uvelsurf"][:])
ub = np.squeeze(nc.variables["uvelbase"][:])
## stuff needed for contour plots
xx = np.tile(x,[len(z),1]) 
zz = (np.tile(z,[len(x),1])).transpose() + b

cts = np.squeeze(nc.variables["cts"][:])
liqfrac = np.squeeze(nc.variables["liqfrac"][:])
temppa = np.squeeze(nc.variables["temp_pa"][:])
mask = np.zeros_like(cts)
mask[zz>s] = 1
cts = np.ma.array(data=cts,mask=mask)
liqfrac = np.ma.array(data=liqfrac,mask=mask)
temppa = np.ma.array(data=temppa,mask=mask)


fig = plt.figure()
ax1=fig.add_subplot(211)
plt.plot(x,us,color='#377EB8', lw = 1.5)
plt.plot(x,ub,'--',color='#377EB8', lw = 1.5)
plt.axis([-250, 3750, 0, 40])
plt.grid()
plt.ylabel("velocity [m a$^{-1}$]")
plt.setp(ax1, xticks=[])
## Contour level of the CTS
cts_level = [1,1]
liqfrac_levels = np.arange(0,2.1,.25)
temppa_levels = np.arange(-6,0,1)

fig.add_subplot(212)
plt.plot(x,b,color='black', lw = 1.5)
plt.plot(x,s,color='black', lw = 1.5)
plt.contourf(xx,zz,liqfrac*100,liqfrac_levels,cmap=plt.cm.Reds,lw = 1)
plt.colorbar(orientation='horizontal',pad=0.10,shrink=0.75)
plt.contourf(xx,zz,temppa,temppa_levels,cmap=plt.cm.Blues_r,lw = 1)
plt.colorbar(orientation='horizontal',pad=0.35,shrink=0.75)
plt.contour(xx,zz,cts,cts_level,colors='black',linestyles='dashed',lw = 1)
plt.axis([-250, 3750, 1000, 2000])
plt.xlabel("distance from bergschrund [m]")
plt.ylabel("elevation [m a.s.l.]")
plt.savefig('sg_results.pdf',bbox_inches='tight',pad_inches=0.35)

nc.close()
