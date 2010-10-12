#!/usr/bin/env python

from time import sleep
from sys import exit
from fcntl import flock, LOCK_SH, LOCK_EX, LOCK_UN
from optparse import OptionParser
from numpy import zeros_like, double

try:
    from netCDF4 import Dataset as NC
except:
    from netCDF3 import Dataset as NC

def lock(name, op):
    f = None

    while (not f):
        try:
            f = open(name)
        except:
            print "File '%s' does not exist; will try again in 5 seconds..." % name
            sleep(5)

    flock(f, op)

    return f

def unlock(f):
    flock(f, LOCK_UN)
    f.close()

def read_data(name, t_old):
    print "Reading data from %s..." % name

    file_is_ready = False

    while not file_is_ready:
        f = lock(name, LOCK_SH)

        nc = NC(name, 'r')
        t = nc.variables['t'][0]

        if t < t_old:
            print "File %s is not ready (time found: %f, time expected: %f); waiting 5 seconds..." % (
                name, t, t_old)
            sleep(5)
        else:
            file_is_ready = True

        nc.close
        unlock(f)

    f = lock(name, LOCK_SH)

    nc = NC(name, 'r')

    x = nc.variables['x'][:]
    y = nc.variables['y'][:]
    t = nc.variables['t'][:]

    usurf = nc.variables['usurf'][:]
    topg  = nc.variables['topg'][:]

    unlock(f)

    return (x, y, t, usurf, topg)

def write_data(data, name):
    print "Writing data to %s..." % name

    (x, y, t, acab) = data

    # create the file if it does not exist
    nc = NC(name, 'w')
    f = lock(name, LOCK_EX)

    # pretent to be taking a very long time
    sleep(5)


    # try to create variables; if it fails (i.e. if they exist already), do nothing
    try:
        t_dim = nc.createDimension('t', None)
        x_dim = nc.createDimension('x', len(x))
        y_dim = nc.createDimension('y', len(y))
        t_var = nc.createVariable('t', 'f8', ('t',))
        x_var = nc.createVariable('x', 'f8', ('x',))
        y_var = nc.createVariable('y', 'f8', ('y',))
        acab_var = nc.createVariable('acab', 'f8', ('t', 'y','x'))

        t_var.units = "years since 1-1-1"
        x_var.units = "m"
        y_var.units = "m"

        acab_var.long_name     = "ice-equivalent surface mass balance (accumulation/ablation) rate"
        acab_var.standard_name = "land_ice_surface_specific_mass_balance"
        acab_var.units         = "m year-1"
    except:
        pass

    t_var[0] = t
    x_var[:] = x
    y_var[:] = y
    acab_var[:] = acab

    nc.close()

    unlock(f)

def compute_acab(usurf, topg, dt):
    return zeros_like(usurf)
    
parser = OptionParser()

parser.usage = "%prog [options]"

parser.description = "Pretends to be a surface model; use with '-surface external'."
parser.add_option("-i", "--input", dest="input",
                  help="The file to read usurf and topg from")
parser.add_option("-o", "--output", dest="output",
                  help="The file to read write acab to")
parser.add_option("--ys", dest="ys",
                  help="Start year")
parser.add_option("--dt", dest="dt",
                  help="Time step")

(options, args) = parser.parse_args()


if not options.input:
    print "Please specify the -i file"
    exit()

if not options.output:
    print "Please specify the -o file"
    exit()

if not options.dt:
    print "Please specify the -dt file"
    exit()

# read data, do stuff, write data...
dt = double(options.dt)
ys = double(options.ys)
t_new = ys
while True:
    x,y,t,usurf,topg = read_data(options.input, t_new)

    acab = compute_acab(usurf, topg, dt)

    write_data((x,y,t_new,acab), options.output)

    t_new = t + dt
