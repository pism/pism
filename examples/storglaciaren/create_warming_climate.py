#!/usr/bin/env python
# Copyright (C) 2017 Andy Aschwanden

import numpy as np
import time
from netCDF4 import Dataset as NC
from argparse import ArgumentParser


# Set up the option parser
parser = ArgumentParser()
parser.description = "Create climate forcing for a warming climate"
parser.add_argument("FILE", nargs='*')
parser.add_argument("-T_max", dest="T_max", type=float,
                    help="Maximum temperature", default=1)
parser.add_argument("-t_max", dest="t_max", type=float,
                    help="lower time bound for maximum temperature", default=100)
parser.add_argument("-amplitude", dest="amplitude", type=float,
                    help="Amplitde of seasonal cycle.", default=12)


options = parser.parse_args()
args = options.FILE
start = 0
end = 1000
step = 1./12.
amplitude = options.amplitude
t_max = options.t_max
T_max = options.T_max
bnds_interval_since_refdate = np.linspace(start, end, end * 12 + 1)
time_interval_since_refdate = (bnds_interval_since_refdate[0:-1] +
                               np.diff(bnds_interval_since_refdate) / 2)

infile = args[0]

nc = NC(infile, 'w')


def def_var(nc, name, units):
    var = nc.createVariable(name, 'f', dimensions=('time'))
    var.units = units
    return var


# create a new dimension for bounds only if it does not yet exist
time_dim = "time"
if time_dim not in list(nc.dimensions.keys()):
    nc.createDimension(time_dim)

# create a new dimension for bounds only if it does not yet exist
bnds_dim = "nb2"
if bnds_dim not in list(nc.dimensions.keys()):
    nc.createDimension(bnds_dim, 2)

# variable names consistent with PISM
time_var_name = "time"
bnds_var_name = "time_bnds"

# create time variable
time_var = nc.createVariable(time_var_name, 'd', dimensions=(time_dim))
time_var[:] = time_interval_since_refdate
time_var.bounds = bnds_var_name
time_var.units = 'years since 1-1-1'
time_var.calendar = '365_day'
time_var.standard_name = time_var_name
time_var.axis = "T"

# create time bounds variable
time_bnds_var = nc.createVariable(bnds_var_name, 'd', dimensions=(time_dim, bnds_dim))
time_bnds_var[:, 0] = bnds_interval_since_refdate[0:-1]
time_bnds_var[:, 1] = bnds_interval_since_refdate[1::]

var = 'delta_T'
dT_var = def_var(nc, var, "K")
T_0 = 0.

temp = np.zeros_like(time_interval_since_refdate) + T_max
temp[0:int(t_max/step)] = np.linspace(T_0, T_max, t_max / step)
temp[:] += -np.cos(time_interval_since_refdate * 2 * np.pi) * amplitude
dT_var[:] = temp

# writing global attributes
script_command = ' '.join([time.ctime(), ':', __file__.split('/')[-1],
                           ' '.join([str(x) for x in args])])
nc.history = script_command
nc.Conventions = "CF 1.6"
nc.close()
