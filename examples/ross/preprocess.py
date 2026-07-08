#!/usr/bin/env python3
import os
import uuid
import subprocess
from sys import exit

try:
    import xarray as xr
except ImportError:
    print("xarray is not installed!")
    exit(1)

import numpy as np

# This seems to be needed by NCO:
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"

smb_name = "climatic_mass_balance"
temp_name = "ice_surface_temp"


def run(commands):
    """Run a list of commands (or one command given as a string)."""
    if isinstance(commands, (list, tuple)):
        for cmd in commands:
            print("Running '%s'..." % cmd)
            subprocess.call(cmd.split(' '))
    else:
        run([commands])


def _open(path):
    """Load a NetCDF file fully into memory and detach the file handle so
    we can rewrite the same path safely."""
    ds = xr.open_dataset(path, decode_times=False, decode_cf=False).load()
    ds.close()
    return ds


def _save_inplace(ds, path):
    """Crash-safe in-place rewrite via tmp + os.replace."""
    tmp = path + ".tmp"
    ds.to_netcdf(tmp, mode="w")
    os.replace(tmp, path)


def preprocess_ice_velocity():
    """
    Download and preprocess the latest Antarctic ice velocity dataset from NASA MEASURES project
    https://nsidc.org/data/NSIDC-0484/versions/2
    """
    input_filename = "antarctica_ice_velocity_450m_v2.nc"

    output_filename = os.path.splitext(input_filename)[0] + "_cutout.nc"

    if not os.path.exists(input_filename):
        print("Please download the InSAR velocity dataset from https://nsidc.org/data/NSIDC-0484/versions/2")
        print("See overview of MEaSUREs data products at https://nsidc.org/data/measures/aiv")
        exit(1)


    ds = _open(input_filename)

    # Create x and y coordinate variables and set projection parameters; cut
    # out the Ross area.

    # Metadata provided with the dataset describes the *full* grid, so it is a
    # lot easier to modify this file instead of adding grid information to the
    # "cutout" file.
    if 'x' not in ds.variables and 'y' not in ds.variables:
        nx = ds.attrs['nx']
        ny = ds.attrs['ny']
        x_min = float(ds.attrs['xmin'].strip().split(' ')[0])
        y_max = float(ds.attrs['ymax'].strip().split(' ')[0])
        x_max = y_max
        y_min = x_min

        x = np.linspace(x_min, x_max, nx)
        y = np.linspace(y_max, y_min, ny)

        ds.attrs["projection"] = "+proj=stere +ellps=WGS84 +datum=WGS84 +lon_0=0 +lat_0=-90 +lat_ts=-71 +units=m"

        ds = ds.assign_coords({
            "x": ("x", x, {"units": "meters",
                            "standard_name": "projection_x_coordinate"}),
            "y": ("y", y, {"units": "meters",
                            "standard_name": "projection_y_coordinate"}),
        })

        _save_inplace(ds, input_filename)

    if not os.path.exists(output_filename):

        #cmd = "ncatted -O -a grid_mapping,VX,d,, -a grid_mapping,VY,d,, %s" % (input_filename)
        #run(cmd)

        # modify this command to cut-out a different region
        cmd = "ncks -d x,3500,7500 -d y,6000,10000 -O %s %s" % (input_filename, output_filename)
        run(cmd)

        tmp_filename = "tmp-{}.nc".format(uuid.uuid4())
        run("cdo setmisstoc,0 -setattribute,VX@_FillValue=-1.0f -setattribute,VY@_FillValue=-1.0f {} {}".format(output_filename, tmp_filename))
        run("cdo chname,VX,vx -chname,VY,vy {} {}".format(tmp_filename, output_filename))

        ds = _open(output_filename)

        # fix units of 'vx' and 'vy'
        ds["vx"].attrs["units"] = "m / year"
        ds["vy"].attrs["units"] = "m / year"

        # Compute and save the velocity magnitude
        if 'magnitude' not in ds.variables:
            vx = np.asarray(ds["vx"].values)
            vy = np.asarray(ds["vy"].values)
            mag = np.sqrt(vx ** 2 + vy ** 2)
            ds["v_magnitude"] = xr.DataArray(
                mag, dims=("y", "x"), attrs={"units": "m / year"})

        _save_inplace(ds, output_filename)

    return output_filename


def preprocess_albmap():
    """
    Download and preprocess the ~16Mb ALBMAP dataset from https://doi.pangaea.de/10.1594/PANGAEA.734145
    """
    url = "https://store.pangaea.de/Publications/LeBrocq_et_al_2010/ALBMAPv1.nc.zip"
    input_filename = "ALBMAPv1.nc"
    output_filename = os.path.splitext(input_filename)[0] + "_cutout.nc"

    commands = ["wget -nc %s" % url,                # download
                "unzip -n %s.zip" % input_filename,  # unpack
                # modify this command to cut out a different region
                "ncks -O -d x1,439,649 -d y1,250,460 %s %s" % (input_filename, output_filename),  # cut out
                "ncks -O -v usrf,lsrf,topg,temp,acca,mask %s %s" % (output_filename, output_filename),  # trim
                "ncrename -O -d x1,x -d y1,y -v x1,x -v y1,y %s" % output_filename,  # fix metadata
                "ncrename -O -v temp,%s -v acca,%s %s" % (temp_name, smb_name, output_filename)]

    run(commands)

    ds = _open(output_filename)

    # fix acab
    rho_ice = 910.0             # kg m-3
    SMB = np.asarray(ds[smb_name].values)
    SMB[SMB == -9999] = 0
    # convert from m/year to kg m-2 / year:
    ds[smb_name].attrs["standard_name"] = "land_ice_surface_specific_mass_balance_flux"
    ds[smb_name].attrs["units"] = "kg m-2 / year"
    ds[smb_name].values[...] = SMB * rho_ice

    # fix artm and topg
    ds[temp_name].attrs["units"] = "degree_Celsius"
    ds["topg"].attrs["standard_name"] = "bedrock_altitude"

    # compute ice thickness
    if 'thk' not in ds.variables:
        usrf = np.asarray(ds["usrf"].values)
        lsrf = np.asarray(ds["lsrf"].values)
        ds["thk"] = xr.DataArray(
            usrf - lsrf, dims=("y", "x"),
            attrs={"units": "meters",
                   "standard_name": "land_ice_thickness"})

    ds.attrs["projection"] = (
        "+proj=stere +ellps=WGS84 +datum=WGS84 +lon_0=0 +lat_0=-90 "
        "+lat_ts=-71 +units=m"
    )
    _save_inplace(ds, output_filename)

    # Remove usrf and lsrf variables:
    command = "ncks -x -v usrf,lsrf -O %s %s" % (output_filename, output_filename)
    run(command)

    return output_filename


def final_corrections(filename):
    """
    * replaces missing values with zeros
    * computes Dirichlet B.C. locations
    """
    ds = _open(filename)

    # replace missing values with zeros
    for var in ['u_bc', 'v_bc', 'magnitude']:
        tmp = np.asarray(ds[var].values)
        fill = ds[var].attrs.get("_FillValue", ds[var].encoding.get("_FillValue"))
        if fill is not None:
            tmp[tmp == fill] = 0
        else:
            tmp[np.isnan(tmp)] = 0
        ds[var].values[...] = tmp

    thk = np.asarray(ds["thk"].values).copy()
    topg = np.asarray(ds["topg"].values).copy()

    # compute the grounded/floating mask:
    mask = np.zeros(thk.shape, dtype='i')

    def is_grounded(thickness, bed):
        rho_ice = 910.0
        rho_seawater = 1028.0
        return bed + thickness > 0 + (1 - rho_ice / rho_seawater) * thickness

    grounded_icy = 0
    grounded_ice_free = 1
    ocean_icy = 2
    ocean_ice_free = 3

    My, Mx = thk.shape
    for j in range(My):
        for i in range(Mx):
            if is_grounded(thk[j, i], topg[j, i]):
                if thk[j, i] > 1.0:
                    mask[j, i] = grounded_icy
                else:
                    mask[j, i] = grounded_ice_free
            else:
                if thk[j, i] > 1.0:
                    mask[j, i] = ocean_icy
                else:
                    mask[j, i] = ocean_ice_free

    # compute the B.C. locations:
    vel_bc_mask = np.logical_or(mask == grounded_icy, mask == grounded_ice_free)

    # mark ocean_icy cells next to grounded_icy ones too:
    row = np.array([0, -1, 1,  0])
    col = np.array([-1,  0, 0,  1])
    for j in range(1, My - 1):
        for i in range(1, Mx - 1):
            nearest = mask[j + row, i + col]

            if mask[j, i] == ocean_icy and np.any(nearest == grounded_icy):
                vel_bc_mask[j, i] = 1

    # Do not prescribe SSA Dirichlet B.C. in ice-free ocean areas:
    vel_bc_mask[thk < 1.0] = 0

    # modifications for the prognostic run
    # this is to avoid grounding in the ice-shelf interior to make the results comparable to the diagnostic flow field
    topg[np.logical_or(mask == ocean_icy, mask == ocean_ice_free)] = -2000.0

    # cap temperature out in the ocean:
    temperature = np.asarray(ds[temp_name].values).copy()
    temperature[temperature > -20.0] = -20.0

    ds[temp_name].values[...] = temperature
    ds["topg"].values[...] = topg
    ds["vel_bc_mask"] = xr.DataArray(vel_bc_mask.astype("i4"),
                                     dims=("y", "x"))
    bad_bc_mask = np.logical_and(thk < 1.0, vel_bc_mask == 1)
    ds["bad_bc_mask"] = xr.DataArray(bad_bc_mask.astype("i4"),
                                     dims=("y", "x"))
    ds["mask"] = xr.DataArray(mask.astype("i4"), dims=("y", "x"))

    _save_inplace(ds, filename)


if __name__ == "__main__":

    velocity = preprocess_ice_velocity()
    albmap = preprocess_albmap()
    albmap_velocity = os.path.splitext(albmap)[0] + "_velocity.nc"  # ice velocity on the ALBMAP grid
    output = "Ross_combined.nc"

    commands = ["./nc2cdo.py %s" % velocity,
                "./nc2cdo.py %s" % albmap,
                "cdo remapbil,%s %s %s" % (albmap, velocity, albmap_velocity),
                "ncks -x -v mask -O %s %s" % (albmap, output),
                "ncks -v vx,vy,v_magnitude -A %s %s" % (albmap_velocity, output),
                "ncrename -v vx,u_bc -v vy,v_bc -v v_magnitude,magnitude -O %s" % output]
    run(commands)

    final_corrections(output)
