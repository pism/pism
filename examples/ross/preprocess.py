#!/usr/bin/env python

try:
    import netCDF4 as NC
except:
    import netCDF3 as NC

import subprocess
import numpy as np
import os

# used to compute latitude and longitude fields
import pyproj

# used to identify the "continental" part of the grounded ice (locations of
# Dirichlet boundary conditions)
from PIL import Image, ImageDraw

def preprocess_ice_velocity():
    """
    Download and preprocess the antarctic ice velocity dataset from
    http://nsidc.org/data/nsidc-0484.html
    """
    url = "ftp://anonymous@sidads.colorado.edu/pub/DATASETS/nsidc0484_MEASURES_antarc_vel_V01/"
    input_filename = "Antarctica_ice_velocity.nc"
    output_filename = os.path.splitext(input_filename)[0] + "_cutout.nc"

    commands = ["wget -nc %s%s.gz" % (url, input_filename), # NSIDC supports compression on demand!
                "gunzip %s.gz" % input_filename,
                "ncrename -d nx,x -d ny,y -O %s %s" % (input_filename, input_filename)
                ]

    try:
        os.stat(input_filename)
    except:
        for cmd in commands:
            print "Running '%s'..." % cmd
            subprocess.call(cmd.split(' '))

    nc = NC.Dataset(input_filename, 'a')

    # Create x and y coordinate variables and set projection parameters; cut
    # out the Ross area.

    # Metadata provided with the dataset describes the *full* grid, so it is a
    # lot easier to modify this file instead of adding grid information to the
    # "cutout" file.
    if 'x' not in nc.variables and 'y' not in nc.variables:
        nx = nc.nx
        ny = nc.ny
        x_min = float(nc.xmin.strip().split(' ')[0])
        y_max = float(nc.ymax.strip().split(' ')[0])
        x_max = y_max
        y_min = x_min

        x = np.linspace(x_min, x_max, nx)
        y = np.linspace(y_max, y_min, ny)

        nc.projection = "+proj=stere +ellps=WGS84 +datum=WGS84 +lon_0=0 +lat_0=-90 +lat_ts=-71 +units=m"

        try:
            x_var = nc.createVariable('x', 'f8', ('x',))
            y_var = nc.createVariable('y', 'f8', ('y',))
        except:
            x_var = nc.variables['x']
            y_var = nc.variables['y']

        x_var[:] = x
        y_var[:] = y

        x_var.units = "meters"
        x_var.standard_name = "projection_x_coordinate"

        y_var.units = "meters"
        y_var.standard_name = "projection_y_coordinate"

        nc.close()

        cmd = "ncks -d x,2200,3700 -d y,3500,4700 -O %s %s" % (input_filename, output_filename)
        subprocess.call(cmd.split(' '))

    # Compute and save the velocity magnitude, latitude and longitude
    nc = NC.Dataset(output_filename, 'a')

    # fix units of 'vx' and 'vy'
    nc.variables['vx'].units = "m / year"
    nc.variables['vy'].units = "m / year"

    if 'v_magnitude' not in nc.variables:
        vx = nc.variables['vx'][:]
        vy = nc.variables['vy'][:]

        v_magnitude = np.zeros_like(vx)

        v_magnitude = np.sqrt(vx**2 + vy**2)

        magnitude = nc.createVariable('v_magnitude', 'f8', ('y', 'x'))
        magnitude.units = "m / year"

        magnitude[:] = v_magnitude

    if 'lon' not in nc.variables:
        grid_shape = nc.variables['vx'].shape

        lon = np.zeros(grid_shape, dtype='f8')
        lat = np.zeros_like(lon)

        p = pyproj.Proj(nc.projection.encode("ASCII"))

        x = nc.variables['x'][:]
        y = nc.variables['y'][:]

        xx,yy = np.meshgrid(x, y)

        lon,lat = p(xx, yy, inverse=True)

        lat_var = nc.createVariable('lat', 'f8', ('y', 'x'))
        lon_var = nc.createVariable('lon', 'f8', ('y', 'x'))

        lon_var[:] = lon
        lat_var[:] = lat

        lat_var.units = "degrees"
        lon_var.units = "degrees"

    nc.close()
    return output_filename

def preprocess_albmap():
    """
    Download and preprocess the ALBMAP dataset from http://doi.pangaea.de/10.1594/PANGAEA.734145
    """
    url = "http://www.pangaea.de/Publications/LeBrocq_et_al_2010/ALBMAPv1.nc.zip"
    input_filename = "ALBMAPv1.nc"
    output_filename = os.path.splitext(input_filename)[0] + "_cutout.nc"

    commands = ["wget -nc %s" % url,
                "unzip %s" % input_filename,
                "ncks -O -d x1,435,645 -d y1,250,460 %s %s" % (input_filename, output_filename),
                "ncks -O -v usrf,lsrf,topg,temp,acca,mask %s %s" % (output_filename, output_filename),
                "ncrename -O -d x1,x -d y1,y -v x1,x -v y1,y %s %s" % (output_filename, output_filename),
                "ncrename -O -v temp,ice_surface_temp -v acca,climatic_mass_balance %s %s" % (output_filename, output_filename)]

    for cmd in commands:
        print "Running '%s'..." % cmd
        subprocess.call(cmd.split(' '))

    nc = NC.Dataset(output_filename, 'a')

    # fix acab
    acab = nc.variables['climatic_mass_balance']
    acab.units = "m / year"
    SMB = acab[:]
    SMB[SMB == -9999] = 0
    acab[:] = SMB

    # fix artm
    artm = nc.variables['ice_surface_temp']
    artm.units = "Celsius"

    # compute ice thickness
    if 'thk' not in nc.variables:
        usrf = nc.variables['usrf'][:]
        lsrf = nc.variables['lsrf'][:]

        thk = nc.createVariable('thk', 'f8', ('y', 'x'))
        thk.units = "meters"

        thk[:] = usrf - lsrf

    # compute locations of Dirichlet boundary conditions:
    if 'bcflag' not in nc.variables:
        mask = nc.variables['mask'][:]
        grid_shape = mask.shape

        bcflag = np.zeros(grid_shape, dtype='i')

        # preprocess the mask using the Python Imaging Library
        land = 5

        try:
            img = Image.fromarray(mask)
            ImageDraw.floodfill(img, (200, 200), land)
            mask = np.asarray(img)
            bcflag = mask == land
        except:
            pass

        bcflag_var = nc.createVariable('bcflag', 'i', ('y', 'x'))
        bcflag_var[:] = bcflag

    nc.projection = "+proj=stere +ellps=WGS84 +datum=WGS84 +lon_0=0 +lat_0=-90 +lat_ts=-71 +units=m"
    nc.close()

    # Remove usrf and lsrf variables:
    command = "ncks -x -v usrf,lsrf -O %s %s" % (output_filename, output_filename)
    subprocess.call(command.split(' '))

    return output_filename

if __name__ == "__main__":
    velocity = preprocess_ice_velocity()
    albmap = preprocess_albmap()
    albmap_velocity = os.path.splitext(albmap)[0] + "_velocity.nc" # ice velocity on the ALBMAP grid
    output = "Ross_combined.nc"

    commands = ["nc2cdo.py %s" % velocity,
                "nc2cdo.py %s" % albmap,
                "cdo remapbil,%s %s %s" % (albmap, velocity, albmap_velocity),
                "ncks -x -v mask -O %s %s" % (albmap, output),
                "ncks -v vx,vy,v_magnitude -A %s %s" % (albmap_velocity, output),
                "ncrename -v vx,u_ssa_bc -v vy,v_ssa_bc -O %s %s" % (output, output)]

    for cmd in commands:
        print "Running '%s'..." % cmd
        subprocess.call(cmd.split(' '))
