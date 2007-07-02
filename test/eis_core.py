#! /usr/bin/env python

from numpy import *
from pycdf import *

SUM8_FILE = 'sum89-92-ss09-50yr.stp'
CORE_FILE = 'eis_grn_core.nc'

##### sum89-92-ss09-50yr.stp #####

dim=[]

input = open(SUM8_FILE, 'r')

# remove headers 
for n in zeros(8):
  input.readline()

for num in input.readline().split():
  dim.append(num)

y_bp = zeros(int(dim[0])/2, float32)
d18O = zeros(int(dim[0])/2, float32)
dT = zeros(int(dim[0])/2, float32)
dSea = zeros(int(dim[0])/2, float32)
count = 0
y_count= int(dim[0])/2-1

for num in input.read().split():
  if count%4 == 0:
    y_bp[y_count] = -float(num)
  if count%4 == 1:
    d18O[y_count] = float(num)
    dT[y_count] = 1.5*(d18O[y_count]+35.27)
    dSea[y_count] = -34.83*(d18O[y_count]+1.93)
    y_count = y_count-1
  count = count + 1

print str(count) + ' numbers read.'

# open the nc file to write to
ncfile = CDF(CORE_FILE, NC.WRITE|NC.CREATE|NC.TRUNC)
ncfile.automode()


# define the one dimension
tdim = ncfile.def_dim('t', int(dim[0])/2)

# define variables
tvar = ncfile.def_var('t', NC.FLOAT, (tdim,))
d18Ovar = ncfile.def_var('delta_O_18', NC.FLOAT, (tdim,))
dTvar = ncfile.def_var('delta_T', NC.FLOAT, (tdim,))
dSeavar = ncfile.def_var('delta_Sea', NC.FLOAT, (tdim,))

# define attributes
setattr(tvar, 'units', 'years since 1989')

setattr(dTvar, 'units', 'm')
setattr(dTvar, 'long_name', 'change in surface temperature from the previous time step')
setattr(dTvar, 'interpolation', 'constant_piecewise_forward')

setattr(dSeavar, 'units', 'm')
setattr(dSeavar, 'long_name', 'change in sea level from the previous time step')
setattr(dSeacar, 'interpolation', 'constant_piecewise_forward')

# write data
tvar[:] = y_bp
d18Ovar[:] = d18O
dTvar[:] = dT
dSeavar[:] = dSea

ncfile.close()
