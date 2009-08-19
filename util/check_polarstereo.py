#!/usr/bin/env python
# This script depends on the following non-standard Python packages:
# 1) pyproj (http://code.google.com/p/pyproj/)
# 2) netcdf4-python (http://code.google.com/p/netcdf4-python/)

try:
    from netCDF4 import Dataset as NC
except:
    from netCDF3 import Dataset as NC
from pyproj import Proj
from sys import argv, exit
from numpy import *
from getopt import getopt, GetoptError

def read_vars(filename, var1, var2):
    nc = NC(filename, "r")
    try:
        x = nc.variables[var1][:]
        y = nc.variables[var2][:]
    except:
        print "Can't read %s and %s variables. Exiting..." % (var1, var2)
        exit()
    nc.close()

    return (squeeze(x), squeeze(y))

def get_ps_info(filename):
    nc = NC(filename, "r")
    
    # Get the projection info:
    try:
        ps = nc.variables["polar_stereographic"]
        straight_vertical_longitude_from_pole = ps.straight_vertical_longitude_from_pole
        latitude_of_projection_origin = ps.latitude_of_projection_origin
        standard_parallel = ps.standard_parallel
    except:
        print "Can't read polar stereographic parameters. Exiting..."
        exit()

    proj_string = "+proj=stere +ellps=WGS84 +datum=WGS84 +lon_0=%f +lat_0=%f +lat_ts=%f +units=m" % (
        straight_vertical_longitude_from_pole,
        latitude_of_projection_origin,
        standard_parallel)
    
    nc.close()
    return Proj(proj_string), proj_string

def fix_latlon(filename):
    print "Fixing latitude and longitude..."
    nc = NC(filename, "a")

    # Note the "fundamental transpose".
    x, y = read_vars(filename, "y", "x")
    p,_ = get_ps_info(filename)

    xx, yy = meshgrid(x, y)
    lon, lat = p(xx, yy, inverse=True)

    try:
        lon_var = nc.variables["lon"]
        lat_var = nc.variables["lat"]
    except:
        lon_var = nc.createVariable("lon", 'f', ("x", "y"))
        lat_var = nc.createVariable("lat", 'f', ("x", "y"))

    try:
        lon_var[:] = lon
        lat_var[:] = lat
    except:
        print "Can't write longitude and latitude. Exiting..."
        exit()

    nc.close()

def write_polarstereo(filename, lon_0, lat_0, lat_ts):
    nc = NC(filename, "a")

    name = "polar_stereographic"

    try:
        ps = nc.variables[name]
    except:
        ps = nc.createVariable(name, "i")
    
    ps.grid_mapping_name = name
    ps.straight_vertical_longitude_from_pole = lon_0
    ps.latitude_of_projection_origin = lat_0
    ps.standard_parallel = lat_ts
    
    nc.close()

def fix_xy(filename):    
    nc = NC(filename, "a")

    lon, lat = read_vars(filename, "lon", "lat")
    p,_ = get_ps_info(filename)

    x0, y0 = p(lon[0,0], lat[0,0])
    x1, y1 = p(lon[-1,-1], lat[-1,-1])

    # Note the "fundamental transpose".
    try:
        x = nc.variables["y"]
        y = nc.variables["x"]
    except:
        print "Can't read x and y variables. Exiting..."
        exit()

    x[:] = linspace(x0, x1, x.shape[0])
    y[:] = linspace(y0, y1, y.shape[0])

    nc.close()

def check_polarstereo(filename):
    print "Checking variables..."
    # Set up the projection:
    p,_ = get_ps_info(filename)

    # Read the lat and lon variables:
    lon_in, lat_in = read_vars(filename, "lon", "lat")

    # Read x and y. Note the "fundamental transpose":
    x, y = read_vars(filename, "y", "x")

    # Compute lon and lat corresponding to x and y:
    xx, yy = meshgrid(x,y)
    lon, lat = p(xx, yy, inverse=True)

    # Compare:
    dlon = abs(lon_in - lon).max()
    dlat = abs(lat_in - lat).max()

    # Report:
    if maximum(dlat, dlon) < 1e-5:
        print "%s has x and y variables consistent with lon and lat." % filename
        print "Maximum difference: %e degrees." % maximum(dlon, dlat)
    else:
        print "%s has x and y variables inconsistent with lon and lat." % filename
        print "Maximum difference: %e degrees." % maximum(dlon, dlat)

def print_options_and_exit():
    print """Options:
With arguments:
-i       <input file>
--lon_0  <straight vertical longitude from pole, degrees>
--lat_0  <latitude of projection origin, degrees>
--lat ts <standard parallel, degrees>
Without arguments:
--check      to check if x and y variables are consistent
             with the polar stereographic projection information
             and lon, lat variables
--fix_latlon to fix lat and lon variables using the x and y variables
             plus the projection information form the file
--fix_xy     to fix x and y using lon and lat variables
             plus the projection information form the file
--write_ps   to write the projection information specified using
             command-line arguments
"""
    exit()

try:
    opts, args = getopt(argv[1:], "i:", ["check", "fix_latlon", "fix_xy", "write_ps", "lon_0=", "lat_0=", "lat_ts="])
    filename   = ""
    lon_0_set  = False
    lat_0_set  = False
    lat_ts_set = False
    check      = False
    latlon_set = False
    xy_set     = False
    ps_set     = False
    lon_0      = 0.0
    lat_0      = 90.0
    lat_ts     = 70.0
    for (opt, arg) in opts:
       if opt == "--check":
           check = True
       if opt == "--fix_latlon":
           latlon_set = True
           xy_set     = False
       if opt == "--fix_xy":
           xy_set     = True
           latlon_set = False
       if opt == "--write_ps":
           ps_set = True
       if opt == "--lon_0":
           lon_0_set = True
           lon_0 = double(arg)
       if opt == "--lat_0":
           lat_0_set = True
           lat_0 = double(arg)
       if opt == "--lat_ts":
           lat_ts_set = True
           lat_ts = double(arg)
       if opt == "-i":
           filename = arg

    if filename == "":
        print_options_and_exit()

    if (ps_set & lon_0_set & lat_0_set & lat_ts_set):
        write_ps(filename, lon_0, lat_0, lat_ts)

    if (latlon_set):
        fix_latlon(filename)

    if (xy_set):
        fix_xy(filename)

    if (check):
        check_polarstereo(filename)
except GetoptError:
    print_options_and_exit()
