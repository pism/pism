import numpy as np
import xarray as xr
import cf_xarray.units
import pint_xarray

Lx = 640e3
Ly = 80e3
dx = 1e3
dy = 1e3
x = np.arange(0, Lx + 160e3 + dx, dx)
y = np.arange(0, Ly + dy, dy)

X, Y = np.meshgrid(x, y)

bed_deep = -720

def Bx(x, B0 = -150, B2 = -728.8, B4 = 343.91, B6 = -50.57, x_bar = 300e3):
    """
    Equation 3-4
    """
    
    return B0 + B2 * (x / x_bar)**2  + B4 * (x / x_bar)**4 + B6 * (x / x_bar)**6

def By(y, Ly = 80e3, dc = 500, fc = 4e3, wc = 24e3):
    """
    Equation 2
    """
    return dc / (1 + np.exp(-2*(y -Ly/2 - wc) /fc)) + dc / (1 + np.exp(2*(y -Ly/2 + wc) /fc))

bed = np.maximum(Bx(X) + By(Y), bed_deep)  # Equation 1
bed[X>Lx] = bed_deep
bed[X>Lx+140e3] = -800

# Initial ice thickness
thickness = np.zeros_like(bed) + 100

# Mask to calve off all ice where x>Lx
liafr = np.zeros_like(bed)
liafr[X<=Lx] = 1

x_dim = "x"
y_dim = "y"

coords = {
    x_dim: (
        x_dim,
        x,
        {
            "units": "m",
            "axis": x_dim.upper(),
            "standard_name": "projection_x_coordinate",
            "long_name": f"{x_dim}-coordinate in projected coordinate system",
        },
    ),
    y_dim: (
        y_dim,
        y,
        {
            "units": "m",
            "axis": y_dim.upper(),
            "standard_name": "projection_y_coordinate",
            "long_name": f"{y_dim}-coordinate in projected coordinate system",
        },
    ),
}

ds = xr.Dataset(
    {
        "bed": xr.DataArray(
            data=bed,
            dims=[y_dim, x_dim],
            coords={y_dim: coords[y_dim], x_dim: coords[x_dim]},
            attrs={
                "standard_name": "bedrock_altitude",
                "units": "m"
            },
        ),
        "thickness": xr.DataArray(
            data=thickness,
            dims=[y_dim, x_dim],
            coords={y_dim: coords[y_dim], x_dim: coords[x_dim]},
            attrs={
                "standard_name": "land_ice_thickness",
                "units": "m"
            },
        ),
        "land_ice_area_fraction_retreat": xr.DataArray(
            data=liafr,
            dims=[y_dim, x_dim],
            coords={y_dim: coords[y_dim], x_dim: coords[x_dim]},
            attrs={
                "units": "1",
            },
        ),            
    },
    attrs={"Conventions": "CF-1.8"}
)
ds.to_netcdf("mismip+.nc")


rho_i = xr.DataArray(918.).pint.quantify("kg m^-3")
rho_sw = xr.DataArray(1028.).pint.quantify("kg m^-3")
accum = xr.DataArray(0.3).pint.quantify("m yr^-1") 
cmb = (accum * rho_i)
ice_surface_temp = xr.DataArray(-2).pint.quantify("degC")
A = xr.DataArray(6.338e-25).pint.quantify("Pa^−3 s^−1")
so = xr.DataArray(35).pint.quantify("g/kg")
to = xr.DataArray(1).pint.quantify("degC")

basins = np.zeros_like(bed)
basins[X<Lx+140e3] = 1

to = np.zeros_like(bed)
to[X<=Lx+140e3] = 1

theta_ocean = xr.zeros_like(ds["bed"]) + to
theta_ocean.name = "theta_ocean"
theta_ocean.attrs.update({"units": "degC"})

salinity_ocean =  xr.zeros_like(ds["bed"]) + so
salinity_ocean.name = "salinity_ocean"
salinity_ocean.attrs.update({"units": "g/kg"})

basins =  xr.zeros_like(ds["bed"]) + basins
basins.name = "basins"

ocean_ds = xr.merge([theta_ocean, salinity_ocean, basins])
ocean_ds.to_netcdf("ocean.nc")

ice_surface_temp = xr.zeros_like(ds["bed"]) + to
ice_surface_temp.name = "ice_surface_temp"
ice_surface_temp.attrs.update({"units": "degC"})

climatic_mass_balance = xr.zeros_like(ds["bed"]) + cmb
climatic_mass_balance.name = "climatic_mass_balance"
climatic_mass_balance.attrs.update({"units": "kg m^-2 yr^-1"})

climate_ds = xr.merge([climatic_mass_balance, ice_surface_temp])
climate_ds.to_netcdf("climate.nc")
