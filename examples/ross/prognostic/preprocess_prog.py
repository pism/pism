#!/usr/bin/env python

# Import all necessary modules here so that if it fails, it fails early.
try:
    import netCDF4 as NC
except:
    import netCDF3 as NC

import subprocess
import numpy as np
import os

smb_name = "climatic_mass_balance"
temp_name = "ice_surface_temp"


def run(commands):
    """Run a list of commands (or one command given as a string)."""
    if isinstance(commands, (list, tuple)):
        for cmd in commands:
            print "Running '%s'..." % cmd
            subprocess.call(cmd.split(' '))
    else:
        run([commands])

def preprocess_ice_velocity():
    """
    Download and preprocess the ~95Mb Antarctic ice velocity dataset from NASA MEASURES project
    http://nsidc.org/data/nsidc-0484.html
    """
    url = "ftp://anonymous@sidads.colorado.edu/pub/DATASETS/nsidc0484_MEASURES_antarc_vel_V01/"
    input_filename = "Antarctica_ice_velocity.nc"
    output_filename = os.path.splitext(input_filename)[0] + "_cutout.nc"

    commands = ["wget -nc %s%s.gz" % (url, input_filename), # NSIDC supports compression on demand!
                "gunzip %s.gz" % input_filename,
                "ncrename -d nx,x -d ny,y -O %s %s" % (input_filename, input_filename)
                ]

    if not os.path.exists(input_filename):
        run(commands)

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

    if not os.path.exists(output_filename):
        cmd = "ncks -d x,2200,3700 -d y,3500,4700 -O %s %s" % (input_filename, output_filename)
        run(cmd)

        nc = NC.Dataset(output_filename, 'a')

        # fix units of 'vx' and 'vy'
        nc.variables['vx'].units = "m / year"
        nc.variables['vy'].units = "m / year"

        # Compute and save the velocity magnitude
        if 'magnitude' not in nc.variables:
            vx = nc.variables['vx'][:]
            vy = nc.variables['vy'][:]

            v_magnitude = np.zeros_like(vx)

            v_magnitude = np.sqrt(vx**2 + vy**2)

            magnitude = nc.createVariable('v_magnitude', 'f8', ('y', 'x'))
            magnitude.units = "m / year"

            magnitude[:] = v_magnitude

        nc.close()

    return output_filename

def preprocess_albmap():
    """
    Download and preprocess the ~16Mb ALBMAP dataset from http://doi.pangaea.de/10.1594/PANGAEA.734145
    """
    url = "http://store.pangaea.de/Publications/LeBrocq_et_al_2010/ALBMAPv1.nc.zip"
    input_filename = "ALBMAPv1.nc"
    output_filename = os.path.splitext(input_filename)[0] + "_cutout.nc"

    commands = ["wget -nc %s" % url,                # download
                "unzip -n %s.zip" % input_filename, # unpack
                "ncks -O -d x1,439,649 -d y1,250,460 %s %s" % (input_filename, output_filename), # cut out
                "ncks -O -v usrf,lsrf,topg,temp,acca %s %s" % (output_filename, output_filename), # trim
                "ncrename -O -d x1,x -d y1,y -v x1,x -v y1,y %s" % output_filename, # fix metadata
                "ncrename -O -v temp,%s -v acca,%s %s" % (temp_name, smb_name, output_filename)]

    run(commands)

    nc = NC.Dataset(output_filename, 'a')

    # fix acab
    acab = nc.variables[smb_name]
    acab.units = "m / year"
    acab.standard_name = "land_ice_surface_specific_mass_balance"
    SMB = acab[:]
    SMB[SMB == -9999] = 0
    acab[:] = SMB

    # fix artm and topg
    nc.variables[temp_name].units = "Celsius"
    nc.variables["topg"].standard_name = "bedrock_altitude"

    # compute ice thickness
    if 'thk' not in nc.variables:
        usrf = nc.variables['usrf'][:]
        lsrf = nc.variables['lsrf'][:]

        thk = nc.createVariable('thk', 'f8', ('y', 'x'))
        thk.units = "meters"
        thk.standard_name = "land_ice_thickness"

        thk[:] = usrf - lsrf


    nc.projection = "+proj=stere +ellps=WGS84 +datum=WGS84 +lon_0=0 +lat_0=-90 +lat_ts=-71 +units=m"
    nc.close()

    # Remove usrf and lsrf variables:
    command = "ncks -x -v usrf,lsrf -O %s %s" % (output_filename, output_filename)
    run(command)

    return output_filename

def final_corrections(filename):
    """
    * replaces missing values with zeros
    * computes Dirichlet B.C. locations
    """
    nc = NC.Dataset(filename, 'a')

    # replace missing values with zeros
    for var in ['u_ssa_bc', 'v_ssa_bc', 'magnitude']:
        tmp = nc.variables[var][:]
        tmp[tmp.mask == True] = 0
        nc.variables[var][:] = tmp

    thk = nc.variables['thk'][:]
    topg = nc.variables['topg'][:]

    # compute the grounded/floating mask:
    mask = np.zeros(thk.shape, dtype='i')
    rho_ice = 910.0
    rho_seawater = 1028.0

    ice_free = 0
    grounded = 1
    floating = 2

    My, Mx = thk.shape
    for j in xrange(My):
        for i in xrange(Mx):
            if topg[j,i] + thk[j,i] > 0 + (1 - rho_ice/rho_seawater) * thk[j,i]:
                mask[j,i] = grounded
            else:
                if thk[j,i] < 1:
                    mask[j,i] = ice_free
                else:
                    mask[j,i] = floating

    # compute the B.C. locations:
    bcflag_var = nc.createVariable('bcflag', 'i', ('y', 'x'))
    bcflag_var[:] = mask == grounded

    # mark floating cells next to grounded ones too:
    row = np.array([-1,  0,  1, -1, 1, -1, 0, 1])
    col = np.array([-1, -1, -1,  0, 0,  1, 1, 1])
    for j in xrange(1, My-1):
        for i in xrange(1, Mx-1):
            nearest = mask[j + row, i + col]

            if mask[j,i] == floating and np.any(nearest == grounded):
                bcflag_var[j,i] = 1
                topg[j,i]=-2000

    #modifications for prognostic run
    tempma = nc.variables[temp_name][:]

    for j in xrange(My):
        for i in xrange(Mx):
            if bcflag_var[j,i] == 0:
                topg[j,i]=-2000 # to avoid grounding
            if tempma[j,i] > -20.0:
                tempma[j,i]=-20.0 # to adjust open ocean temperatures

    nc.variables[temp_name][:] = tempma
    nc.variables['topg'][:] = topg

    nc.close()

if __name__ == "__main__":
    velocity = preprocess_ice_velocity()
    albmap = preprocess_albmap()
    albmap_velocity = os.path.splitext(albmap)[0] + "_velocity.nc" # ice velocity on the ALBMAP grid
    output = "Ross_combined_prog.nc"

    commands = ["nc2cdo.py %s" % velocity,
                "nc2cdo.py %s" % albmap,
                "cdo remapbil,%s %s %s" % (albmap, velocity, albmap_velocity),
                "ncks -x -v mask -O %s %s" % (albmap, output),
                "ncks -v vx,vy,v_magnitude -A %s %s" % (albmap_velocity, output),
                "ncrename -v vx,u_ssa_bc -v vy,v_ssa_bc -v v_magnitude,magnitude -O %s" % output]
    run(commands)

    final_corrections(output)
