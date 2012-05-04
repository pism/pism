#!/usr/bin/env python
from numpy import arange, squeeze, zeros, shape
from sys import exit

try:
    from netCDF4 import Dataset as NC
except:
    from netCDF3 import Dataset as NC


## prepare the output file:
def prepare_file(fname, x, y):
    print "  preparing file '%s'..." % fname
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

    # put dimensions metadata in a dictionary:
    #   name : (unit, long_name, standard_name)
    attributes = {"x" : ("m", "X-coordinate in Cartesian system", "projection_x_coordinate"),
                  "y" : ("m", "Y-coordinate in Cartesian system", "projection_y_coordinate"),
                  "time" : ("years since 2004-1-1",  "time", "time")
                  }

    for each in list(attributes.keys()):
        var = nc.variables[each]
        var.units = attributes[each][0]
        var.long_name = attributes[each][1]
        var.standard_name = attributes[each][2]

    ps = nc.createVariable("mapping", 'b')
    ps.grid_mapping_name = "polar_stereographic"
    ps.straight_vertical_longitude_from_pole = -39.0
    ps.latitude_of_projection_origin = 90.0
    ps.standard_parallel = 71.0
    nc.Conventions = "CF-1.3"

    return (nc, t_var)

input_temp = "air_temp.nc"
input_precip = "precipitation.nc"

## read air temperatures, and get space/time grid info from this file:
nc_temp = NC(input_temp, 'r')
temp_in = nc_temp.variables['air_temp']
x = nc_temp.variables['x'][:]
y = nc_temp.variables['y'][:]
N = len(nc_temp.dimensions['time'])
years = N/12                    # number of years covered
print "  found N = %d frames covering %d years in file %s" % (N, years, input_temp)

## read precipitation:
nc_precip = NC(input_precip, 'r')
precip_in = nc_precip.variables['precipitation']
if len(nc_precip.dimensions['time']) != N:
    print "ERROR: number of frames in precipitation file '%s' does not match that in air temperature file '%s'" \
         % (input_temp,input_precip)
    exit(1)
else:
    print "  found same number of frames in file %s" % input_precip

output_temp = "ar4_temp_anomaly.nc"
output_precip = "ar4_precip_anomaly.nc"

nc_temp_out, t_temp = prepare_file(output_temp, x, y)
nc_precip_out, t_precip = prepare_file(output_precip, x, y)

temp = nc_temp_out.createVariable("air_temp_anomaly", 'f', dimensions=("time", "y", "x"))
temp.units = "Kelvin"
temp.long_name = \
  "mean annual air temperature at 2m above the surface (as anomaly from first year)"
temp.mapping = "mapping"
temp.description = "AR4 temperature anomaly"

precip = nc_precip_out.createVariable("precipitation_anomaly", 'f', dimensions=("time", "y", "x"))
precip.units = "m year-1"
precip.long_name = \
  "mean annual ice-equivalent precipitation rate (as anomaly from first year)"
precip.mapping = "mapping"
precip.description = "AR4 precipitation anomaly"

print "  averaging monthly temperature data to give annual mean"
for year in arange(years):
    t_temp[year] = year + 0.5
    temp[year,:,:] = zeros((y.size, x.size))

    for month in arange(12):
        j = 12 * year + month
        #print "    [year=%d,j=%d]" % (year,j)
        temp[year,:,:] += temp_in[j,:,:]

    temp[year,:,:] /= 12.0

print "  averaging monthly precipitation data to give annual mean"
for year in arange(years):
    t_precip[year] = year + 0.5
    precip[year,:,:] = zeros((y.size, x.size))

    for month in arange(12):
       j = 12 * year + month
       precip[year,:,:] += precip_in[j,:,:]

    precip[year,:,:] /= 12.0

# convert to anomalies by subtracting-off first year averages:
print "  converting annual mean temperature to 'air_temp_anomaly'"
temp_0 = temp[0,:,:].copy()
for year in arange(years):
    temp[year,:,:] -= temp_0

print "  converting annual mean precipitation to 'precipitation_anomaly'"
precip_0 = precip[0,:,:].copy()
for year in arange(years):
    precip[year,:,:] -= precip_0

## warning: earlier closure of these input files causes problems; why?
nc_temp.close()
nc_precip.close()

nc_temp_out.close()
nc_precip_out.close()

print "done"

