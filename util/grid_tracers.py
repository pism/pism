#!/usr/bin/env python
from matplotlib.mlab import griddata
from argparse import ArgumentParser
from netCDF4 import Dataset
import numpy as np
import matplotlib.pylab as pl

def split (data, tracer_ids, timeslices):
    subsets = [ np.array([i in ids_t for i in tracer_ids], dtype=bool) for ids_t in timeslices]
    data = [data[i] for i in subsets]
    return data

def extract_sets(filename, created_file):
    cf = Dataset(created_file, "r")

    # split up tracer creations by creation time
    tracer_times = cf.variables["time"][:]
    tracer_ids   = cf.variables["tracer_id"][:]
    timeslices = np.array(sorted(set(tracer_times)))
    print timeslices
    ts_lookup = {t: n for (n, t) in enumerate (timeslices)}
    time_lookup = {tid: time for (tid, time) in zip (tracer_ids, tracer_times)}
    print 1
    ids = [tracer_ids [tracer_times == t ] for t in timeslices]
    print 2 
    # get all tracers from file that were created at a specified time
    df = Dataset(filename, "r")
    t_id = df.variables["tracer_id"][:]
    t_x = df.variables["tracer_x"][:]
    t_y = df.variables["tracer_y"][:]
    t_z = df.variables["tracer_z"][:]
    print 3
    data = np.zeros ([5,len(t_z)])
    data[0] = np.array([ts_lookup[time_lookup[td]] for td in t_id])
    data[1] = t_id
    data[2] = t_x
    data[3] = t_y
    data[4] = t_z
    print 4 
    xi = df.variables["x"][:]
    yi = df.variables["y"][:]

    useful =  np.array([(data[0] == x).any() for x in xrange(len(timeslices)) ], dtype=bool)
    timeslices = np.array(timeslices[useful])
    
    print 5
    arrays = [data[:, data[0] == n ] for n in xrange (len(timeslices))]
    print 5.2
    arrays = [ x for x in arrays if len(x[0]) != 0 ]
    print 5.5
    zd = np.array([ griddata(a[3], a[2], a[4], yi, xi) for a in arrays  ])
    print 6 
    of = Dataset(filename[:-3] + "_layers.nc" , 'w')
    for v in "xy":
        print v
        varin = df.variables[v]
        of.createDimension(v, len(varin[:]) )
        print " dt = " , varin[:].dtype
        print "dims = ", varin.dimensions
        outvar = of.createVariable(v, varin[:].dtype, varin.dimensions)
        outvar.setncatts({k: varin.getncattr(k) for k in varin.ncattrs()})

    of.createDimension("creation_time", len(timeslices) )
    ct = of.createVariable("creation_time", timeslices.dtype, ("creation_time",))
    ct[:] = timeslices
    ct.setncatts({k: df.variables["time"].getncattr(k) for k in df.variables["time"].ncattrs()})
    
    of.createVariable("z", zd[0].dtype, ("creation_time", "x", "y"))
    print of.variables["z"].shape
    print zd.shape
    of.variables["z"][:] = zd
    of.close()
    return (timeslices, ids, xs, ys, zs, zd)



def parse_args():
  parser = ArgumentParser()
  parser.description = "compare slopes of two variables from two files"
  parser.add_argument("FILES", nargs='*')
  parser.add_argument("-v", "--verbose",
                    help='''Be verbose''', action="store_true")
  parser.add_argument("-c", "--created_file",
                      help='''File with particle creation date and time''',
                      default = None)
  # parser.add_argument("-s", "--state",
  #                      help='''file with reference values''', required = True)
  # parser.add_argument("-b", "--var_b",
  #                    help='''variable b''', default="data")
  options = parser.parse_args()
  return options


def main():
  options = parse_args()
  if options.verbose:
    print (dir(options))
    print options.FILES
  extract_sets(options.FILES[0], options.created_file)

if __name__ == "__main__":
    main()
