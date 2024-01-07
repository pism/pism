#!/usr/bin/env python3

# creates an example flowline data set "slab.nc" which can be used with
# flowline.py to show how to run PISM in flow-line mode

# see the "Using PISM for flow-line modeling" subsection of the Users Manual

# We recommend creating more metadata-rich datasets than is done here.
# Here we generate the bare minimum to work with flowline.py and PISM.

from numpy import linspace, minimum, maximum, abs
import xarray as xr
import sys
import time


Lx = 1000e3                     # 1000 km
Mx = 501
topg_slope = -1e-4
thk_0 = 1e3                        # meters
climatic_mass_balance_0 = 0             # kg m-2 s-1
ice_surface_temp_0 = -10           # Celsius

nc = netCDF4.Dataset("slab.nc", 'w')
nc.createDimension('x', Mx)

x = nc.createVariable('x',    'f4', ('x',))
topg = nc.createVariable('topg', 'f4', ('x',))
thk = nc.createVariable('thk',  'f4', ('x',))
climatic_mass_balance = nc.createVariable('climatic_mass_balance', 'f4', ('x',))
ice_surface_temp = nc.createVariable('ice_surface_temp', 'f4', ('x',))

x[:] = linspace(-Lx, Lx, Mx)
x.units = "m"
topg[:] = topg_slope * (x[:] - Lx)
topg.units = "m"

thk[:] = maximum(minimum(5e3 - abs(x[:]) * 0.01, thk_0), 0)
thk.units = "m"

climatic_mass_balance[:] = climatic_mass_balance_0
climatic_mass_balance.units = "kg m-2 s-1"
ice_surface_temp[:] = ice_surface_temp_0
ice_surface_temp.units = "Celsius"

nc.close()

historysep = ' '
historystr = time.asctime() + ': ' + historysep.join(sys.argv) + '\n'

coords = {
    "x": (
        ["x"],
        x,
        {
            "units": "m",
            "axis": "X",
            "standard_name": "projection_x_coordinate",
            "long_name": "x-coordinate in projected coordinate system",
        },
    ),
    "y": (
        ["y"],
        y,
        {
            "units": "m",
            "axis": "Y",
            "standard_name": "projection_y_coordinate",
            "long_name": "y-coordinate in projected coordinate system",
        },
    ),
}

ds = xr.Dataset(
    {
        "topg": xr.DataArray(
            data=topg,
            dims=["y", "x"],
            coords=coords,
            attrs={
                "standard_name": "bedrock_altitude",
                "units": "m",
            },
        ),
        "thk": xr.DataArray(
            data=thk,
            dims=["y", "x"],
            coords=coords,
            attrs={
                "standard_name": "land_ice_thickness",
                "units": "m",
            },
        ),
        "ice_surface_temp": xr.DataArray(
            data=ice_surface_temp,
            dims=["y", "x"],
            coords=coords,
            attrs={
                "units": "Celsius",
            },
        ),
        "climatic_mass_balance": xr.DataArray(
            data=climatic_mass_balance,
            dims=["y", "x"],
            coords=coords,
            attrs={
                "standard_name": "land_ice_surface_specific_mass_balance",
                "units": "kg m-2 s-1",
            },
        ),
    },
    attrs={"Conventions": "CF-1.7", "history": historystr},
)
ds.to_netcdf("slab.nc")
