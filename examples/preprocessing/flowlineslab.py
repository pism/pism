#!/usr/bin/env python3

# creates an example flowline data set "slab.nc" which can be used with
# pism_flowline to show how to run PISM in flow-line mode

# see the "Using PISM for flow-line modeling" subsection of the Users Manual

# We recommend creating more metadata-rich datasets than is done here.
# Here we generate the bare minimum to work with pism_flowline and PISM.

from numpy import linspace, minimum, maximum, abs
import xarray as xr

Lx = 1000e3                     # 1000 km
Mx = 501
topg_slope = -1e-4
thk_0 = 1e3                        # meters
climatic_mass_balance_0 = 0             # kg m-2 s-1
ice_surface_temp_0 = -10           # Celsius

x = linspace(-Lx, Lx, Mx)
topg = topg_slope * (x - Lx)
thk = maximum(minimum(5e3 - abs(x) * 0.01, thk_0), 0)

ds = xr.Dataset(
    coords={"x": ("x", x.astype("f4"), {"units": "m"})},
    data_vars={
        "topg": (("x",), topg.astype("f4"), {"units": "m"}),
        "thk":  (("x",), thk.astype("f4"),  {"units": "m"}),
        "climatic_mass_balance": (("x",),
                                  (climatic_mass_balance_0 * abs(x) * 0).astype("f4"),
                                  {"units": "kg m^-2 s^-1"}),
        "ice_surface_temp": (("x",),
                             (ice_surface_temp_0 + abs(x) * 0).astype("f4"),
                             {"units": "degree_Celsius"}),
    },
)
ds.to_netcdf("slab.nc", mode="w")
