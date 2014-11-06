#!/usr/bin/env python

# creates flowline data set "fevor-slab.nc" which can be used with
# flowline.py to setup a flow-line fevor run. 

# see the "Using PISM for flow-line modeling" subsection of the Users Manual

from numpy import linspace, minimum, maximum, abs
import netCDF4

Lx         = 10e3                       # 10 km
Mx         = 201
# topg_slope = -0.08748866352             # 5 degrees
topg_slope = 0                          # zero degrees (we prescribe the surface gradient)
thk_0      = 1e3                        # 1000 meters
climatic_mass_balance_0 = 0             # kg m-2 s-1
ice_surface_temp_0      = -60           # Celsius

nc = netCDF4.Dataset("fevor-slab.nc", 'w')
nc.createDimension('x', Mx)

x    = nc.createVariable('x',    'f4', ('x',))
topg = nc.createVariable('topg', 'f4', ('x',))
thk  = nc.createVariable('thk',  'f4', ('x',))
climatic_mass_balance = nc.createVariable('climatic_mass_balance', 'f4', ('x',))
ice_surface_temp = nc.createVariable('ice_surface_temp', 'f4', ('x',))

x[:]    = linspace(-Lx, Lx, Mx);    x.units = "m"
topg[:] = topg_slope * (x[:] - Lx); topg.units = "m"

thk[:] =  x[:]*0 + thk_0
thk.units = "m"

climatic_mass_balance[:] = climatic_mass_balance_0;
climatic_mass_balance.units = "kg m-2 s-1"
ice_surface_temp[:] = ice_surface_temp_0;
ice_surface_temp.units = "Celsius"

nc.close()
