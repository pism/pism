#!/usr/bin/env python
from numpy import arange, linspace, squeeze, meshgrid, shape
from pyproj import Proj
from griddata import griddata
from sys import stderr

try:
    from netCDF4 import Dataset as NC
except:
    from netCDF3 import Dataset as NC

def read_usurf(input):
    nc = NC(input, 'r')
    usurf = squeeze(nc.variables['ht'][:])
    nc.close()
    return usurf

def read_input(input):
    nc = NC(input, 'r')

    time_hours = nc.variables['time'][:];
    time_years = time_hours / 24.0 / 360.0;
    # ALTERNATE: time_months = arange(0.0, 95.0 * 12.0, 1.0)
    # COMMENT:  conversion to months could/should be based on udunits2:
    #  You have: 1 month
    #  You want: hours
    #      1 month = 730.484 hours

    latitude = nc.variables['latitude'][:]
    longitude = nc.variables['longitude'][:]
    usurf = nc.variables['topg'][:]
    precipitation = nc.variables['precip'][:]
    temp = nc.variables['temp'][:]
    nc.close()

    return (time_years, latitude, longitude, usurf, precipitation, temp)

def prepare_file(fname, x, y):
    write("Writing data to '%s'...\n" % fname)
    nc = NC(fname, "w",format="NETCDF3_CLASSIC")
    nc.set_fill_off()
    nc.createDimension("time", None)
    nc.createDimension("x", size=x.size)
    nc.createDimension("y", size=y.size)

    t_var = nc.createVariable("time", 'f', dimensions=("time",))

    x_var = nc.createVariable("x", 'f', dimensions=("x",))
    x_var[:] = x

    y_var = nc.createVariable("y", 'f', dimensions=("y",))
    y_var[:] = y

    # put dimension metadata in a dictionary:
    #    name : (unit, long_name, standard_name)
    attributes = {"x" : ("m", "X-coordinate in Cartesian system", "projection_x_coordinate"),
                  "y" : ("m", "Y-coordinate in Cartesian system", "projection_y_coordinate"),
                  "time" : ("years since 2004-1-1",  "time", "time")
                  }

    for each in list(attributes.keys()):
        var = nc.variables[each]
        var.units = attributes[each][0]
        var.long_name = attributes[each][1]
        var.standard_name = attributes[each][2]

    nc.variables['time'].calendar = "360_day"

    ps = nc.createVariable("mapping", 'b')
    ps.grid_mapping_name = "polar_stereographic"
    ps.straight_vertical_longitude_from_pole = -39.0
    ps.latitude_of_projection_origin = 90.0
    ps.standard_parallel = 71.0
    nc.Conventions = "CF-1.3"

    return (nc, t_var)

def prepare_precipitation_file(output, x, y):
    (nc, t) = prepare_file(output, x, y)

    precipitation = nc.createVariable("precipitation", 'f', dimensions=("time", "y", "x"))
    precipitation.units = "m/a"
    precipitation.long_name = "ice-equivalent precipitation rate"
    precipitation.grid_mapping = "mapping"
    precipitation.description = "AR4 precipitation anomaly"
    return (nc, t, precipitation)

def prepare_temp_file(output, x, y):
    (nc, t) = prepare_file(output, x, y)

    temp = nc.createVariable("air_temp", 'f', dimensions=("time", "y", "x"))
    temp.units = "Kelvin"
    temp.long_name = "air temperature at 2m above the surface"
    temp.mapping = "mapping"
    temp.description = "AR4 temperature anomaly"

    usurf = nc.createVariable("usurf", 'f', dimensions=("y", "x"))
    usurf.units = "m"
    usurf.long_name = "surface elevation"
    usurf.mapping = "mapping"
    usurf.standard_name = "surface_altitude"

    return (nc, t, temp, usurf)

write   = stderr.write

rho_ice = 910.0                  # ice density, kg/m3
rho_w = 1000.0                   # pure water density, kg/m3; Charles Jackson confirms this density

input = "Climate_forcing_2004_2098_v3.nc"

# output projection:
projection = "+proj=stere +ellps=WGS84 +datum=WGS84 +lon_0=-39 +lat_0=90 +lat_ts=71 +units=m"

# describe the 5km SeaRISE grid, even though we use coarser
e0 =  -800e3 # m
n0 = -3400e3 # m
de =     5e3 # m
dn =     5e3 # m
M  =     301
N  =     561
e1 = e0 + (M-1)*de
n1 = n0 + (N-1)*dn

# pick a coarse (50km) grid; finer than GCM but coarser than PISM runs
M  =     31
N  =     57

x = linspace(e0, e1, M)
y = linspace(n0, n1, N)

time_years, latitude, longitude, usurf, precipitation, temp = read_input(input)

# Set up the interpolation code and extract Greenland:
lon, lat = meshgrid(longitude, latitude)
mask = (lat >= 57.0) & (lat <= 88) & ((lon >= 360.0 - 94.0) | (lon <= 12.0))

h = usurf[:, mask]
T = temp[:, mask]
P = precipitation[:, mask]
lon = lon[mask]
lat = lat[mask]

ee, nn = (Proj(projection))(lon, lat)
ee_out,nn_out = meshgrid(x,y)
usurf_out = griddata(ee,nn,h,ee_out,nn_out)

output_temp = "air_temp.nc"
output_precip = "precipitation.nc"

nc_temp, t_temp_var, temp_var, usurf_var_temp = prepare_temp_file(output_temp, x, y)
usurf_var_temp[:] = usurf_out

nc_precip, t_precip_var, precip_var = prepare_precipitation_file(output_precip, x, y)

write("Interpolating")
# The main interpolation loop:
for j in range(time_years.size):

  f_temp = griddata(ee,nn,T[j],ee_out,nn_out)
  f_precip = griddata(ee,nn,P[j],ee_out,nn_out)
  t_temp_var[j] = time_years[j]
  t_precip_var[j] = time_years[j]
  temp_var[j,:,:] = f_temp
  # convert to m/year ice equivalent:
  precip_var[j,:,:] = f_precip * (rho_w / rho_ice)

  write(".")

write("done.\n")
nc_temp.close()
nc_precip.close()

