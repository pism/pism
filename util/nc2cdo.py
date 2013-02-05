#!/usr/bin/env python

## @package nc2cdo
# \author Andy Aschwanden, University of Alaska Fairbanks, USA
# \brief Script makes netCDF file ready for Climate Data Operators (CDO).
# \details Script adds attributes and variables to a netCDF which are required for
# post-processing with <a href="http://www.mad.zmaw.de/Pingo/post/post.cdo.home.html">CDO</a>.
# That is, the attribute 'coordinates = "lon lat"' is added to any gridded variable
# found in the file, except for lat and lon. For remapping methods other than bilinear,
# CDO additionally needs the corners of each grid cell (the lat/lon variables define the
# cell centers). We compute the corners of each grid cell on the x-y grid because this is
# relatively straightforward (compared to the lat-lon world), and then project the corners
# onto the lat-lon grid. The script attemps to retrieve the required mapping information from,
# if present, a global attribute called 'projection' which (hopefully) contains a valid
# Proj4 projection string. If this string is not available (this is true for standard PISM
# output), the script searches variables for the 'grid_mapping' attribute and translates information
# from the corresponding mapping variable into a Proj4 projection string. If neither is found,
# the script exists with an error message.
#
# Usage:
#
# \verbatim $ nc2cdo.py foo.nc \endverbatim

import numpy as np
from argparse import ArgumentParser

from pyproj import Proj

# try different netCDF modules
try:
    from netCDF4 import Dataset as CDF
except:
    from netCDF3 import Dataset as CDF


## Set up the option parser
parser = ArgumentParser()
parser.description = '''Script makes netCDF file ready for Climate Data Operators (CDO). Either a global attribute "projection", a mapping variable, or a command-line proj4 string or a EPSG code must be given.'''
parser.add_argument("FILE", nargs=1)
parser.add_argument("--srs",dest="srs",
                  help='''
                  a valid proj4 string describing describing the projection
                  ''', default=None)
options = parser.parse_args()
args = options.FILE
srs = options.srs

if len(args) == 1:
    nc_outfile = args[0]
else:
    print('wrong number arguments, 1 expected')
    parser.print_help()
    exit(0)

## Get projection information from netCDF file
def get_projection_from_file(nc):

    ## First, check if we have a global attribute 'projection'
    ## which contains a Proj4 string:
    try:
        p = Proj(str(nc.projection))
        print('Found projection information in global attribute, using it')
    except:
        try:
            ## go through variables and look for 'grid_mapping' attribute
            for var in nc.variables.keys():
                if hasattr(nc.variables[var], 'grid_mapping'):
                    mappingvarname = nc.variables[var].grid_mapping
                    exit
            print('Found projection information in variable "%s", using it' % mappingvarname)
            var_mapping = nc.variables[mappingvarname]
            p = Proj(proj   = "stere",
                     ellps  = var_mapping.ellipsoid,
                     datum  = var_mapping.ellipsoid,
                     units = "m",
                     lat_ts = var_mapping.standard_parallel,
                     lat_0  = var_mapping.latitude_of_projection_origin,
                     lon_0  = var_mapping.straight_vertical_longitude_from_pole,
                     x_0    = var_mapping.false_easting,
                     y_0    = var_mapping.false_northing)
        except:
            print('No mapping information found, exiting.')
            exit(1)

    return p


if __name__ == "__main__": 

    ## open netCDF file in 'append' mode
    nc = CDF(nc_outfile,'a')

    ## a list of possible x-dimensions names
    xdims = ['x', 'x1']
    ## a list of possible y-dimensions names
    ydims = ['y', 'y1']

    ## assign x dimension
    for dim in xdims:
        if dim in nc.dimensions.keys():
            xdim = dim
    ## assign y dimension
    for dim in ydims:
        if dim in nc.dimensions.keys():
            ydim = dim

    ## coordinate variable in x-direction
    x_var = nc.variables[xdim]
    ## coordinate variable in y-direction
    y_var = nc.variables[ydim]

    ## values of the coordinate variable in x-direction
    easting  = x_var[:]
    ## values of the coordinate variable in y-direction
    northing = y_var[:]

    ## grid spacing in x-direction
    de = np.diff(easting)[0]
    ## grid spacing in y-direction
    dn = np.diff(northing)[0]

    ## number of grid points in x-direction
    M = easting.shape[0]
    ## number of grid points in y-direction
    N = northing.shape[0]

    ## number of grid corners
    grid_corners = 4

    ## array holding x-component of grid corners
    gc_easting = np.zeros((M,grid_corners))
    ## array holding y-component of grid corners
    gc_northing = np.zeros((N,grid_corners))
    ## array holding the offsets from the cell centers in x-direction (counter-clockwise)
    de_vec = np.array([-de/2, de/2, de/2, -de/2])
    ## array holding the offsets from the cell centers in y-direction (counter-clockwise)
    dn_vec = np.array([-dn/2, -dn/2, dn/2, dn/2])
    ## array holding lat-component of grid corners
    gc_lat = np.zeros((N, M, grid_corners))
    ## array holding lon-component of grid corners
    gc_lon = np.zeros((N, M, grid_corners))

    if srs:
        ## Use projection from command line
        try:
            proj = Proj(init=srs)
        except:
            proj = Proj(srs)
    else:
        ## Get projection from file
        proj = get_projection_from_file(nc)
    
    ## If it does not yet exist, create dimension 'grid_corners'    
    if 'grid_corners' not in nc.dimensions.keys():
        for corner in range(0,grid_corners):
            ## grid_corners in x-direction
            gc_easting[:,corner] = easting + de_vec[corner]
            ## grid corners in y-direction
            gc_northing[:,corner] = northing + dn_vec[corner]
            ## meshgrid of grid corners in x-y space
            gc_ee, gc_nn = np.meshgrid(gc_easting[:,corner],gc_northing[:,corner])
            ## project grid corners from x-y to lat-lon space
            gc_lon[:,:,corner], gc_lat[:,:,corner] = proj(gc_ee,gc_nn,inverse=True)

        nc.createDimension("grid_corners",size=grid_corners)
    
        var = 'lon_bounds'
        ## Create variable 'lon_bounds'
        var_out = nc.createVariable(var, 'f', dimensions=(ydim, xdim, "grid_corners"))
        ## Assign units to variable 'lon_bounds'
        var_out.units = "degreesE";
        ## Assign values to variable 'lon_bounds'
        var_out[:] = gc_lon

        var = 'lat_bounds'
        ## Create variable 'lat_bounds'
        var_out = nc.createVariable(var, 'f', dimensions=(ydim, xdim, "grid_corners"))
        ## Assign units to variable 'lat_bounds'
        var_out.units = "degreesN";
        ## Assign values to variable 'lat_bounds'
        var_out[:] = gc_lat


    if (not 'lon' in nc.variables.keys()) or ( not 'lat' in nc.variables.keys()):
        print("No lat/lon coordinates found, creating them")
        ee, nn = np.meshgrid(easting,northing)
        lon, lat = proj(ee, nn, inverse=True)

    var = 'lon'
    ## If it does not yet exist, create variable 'lon'
    if not var in nc.variables.keys():
        var_out = nc.createVariable(var, 'f', dimensions=(ydim, xdim))
        ## Assign values to variable 'lon'
        var_out[:] = lon
    else:
        var_out = nc.variables[var]
    ## Assign units to variable 'lon'
    var_out.units = "degreesE"
    ## Assign long name to variable 'lon'
    var_out.long_name = "Longitude"
    ## Assign standard name to variable 'lon'
    var_out.standard_name = "longitude"
    ## Assign bounds to variable 'lon'
    var_out.bounds = "lon_bounds"
    
    var = 'lat'
    ## If it does not yet exist, create variable 'lat'
    if not var in nc.variables.keys():
        var_out = nc.createVariable(var, 'f', dimensions=(ydim, xdim))
        var_out[:] = lat
    else:
        var_out = nc.variables[var]
    ## Assign units to variable 'lat'
    var_out.units = "degreesN";
    ## Assign long name to variable 'lat'
    var_out.long_name = "Latitude"
    ## Assign standard name to variable 'lat'
    var_out.standard_name = "latitude"
    ## Assign bounds to variable 'lat'
    var_out.bounds = "lat_bounds"
    ## Assign values to variable 'lat'

    ## Make sure variables have 'coordinates' attribute
    for var in nc.variables.keys():
        if (nc.variables[var].ndim >= 2):
            nc.variables[var].coordinates = "lon lat"

    ## lat/lon coordinates must not have mapping and coordinate attributes
    ## if they exist, delete them
    for var in ['lat', 'lon', 'lat_bounds', 'lon_bounds']:
        if hasattr(nc.variables[var], 'grid_mapping'):
            delattr(nc.variables[var], 'grid_mapping')
        if hasattr(nc.variables[var], 'coordinates'):
            delattr(nc.variables[var], 'coordinates')

    ## If present prepend history history attribute, otherwise create it
    from time import asctime
    histstr = asctime() + ' : grid infos for CDO added by nc2cdo.py (PISM util) '
    if 'History' in nc.ncattrs():
        nc.History = histstr + nc.History
    elif 'history' in nc.ncattrs():
        nc.history = histstr + nc.history
    else:
        nc.history = histstr
    ## Write projection attribute
    nc.projection = proj.srs
    ## Close file    
    nc.close()
