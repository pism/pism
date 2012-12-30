#!/usr/bin/env python

# FIXME:  for now need
#   $ ln -s ../../util/PISMNC.py
#   $ (cd ../../src/verif/tests/ && make convertP)
#   $ ln -s ../../src/verif/tests/convertP

#PLAN:  construction of a PISM input file (-boot_file)
#  with fields x,y,thk,bwat,bwp.  Then run
#    pismr -boot_file inputforP.nc -Mx 51 -My 51 -Mz 11 -Lz 4000 -hydrology distributed -y 20 -no_mass -no_energy -o end.nc  
#  Then compare
#    nccmp.py -v bwat,bwp end.nc inputforP.nc difftestP.nc
#  and collect norms etc. for difftestP.nc.

import numpy as np
from sys import exit, argv, stderr
from os import system
from PISMNC import PISMDataset

if len(argv)<2:
  pism_path="."
else:
  pism_path=argv[1]
if len(argv)<3:
  mpiexec=""   # or: mpiexec = "mpiexec -n 1"
else:
  mpiexec=argv[2]

stderr.write("Testing: Test P verification of '-hydrology distributed'.\n")

Lx = 25.0e3  # outside L = 22.5 km
Mx = 51      # 1 km grid

x = np.linspace(-Lx, Lx, Mx)
xx, yy = np.meshgrid(x, x)
h = np.zeros(np.shape(xx))
W = np.zeros(np.shape(xx))

# create 1D array of tuples (r,j,k), sorted by r-value
dtype = [('r', float), ('j', int), ('k', int)]
rr = np.empty((Mx,Mx),dtype=dtype)
for j in range(Mx):
  for k in range(Mx):
    rr[j,k] = (np.sqrt(xx[j,k]*xx[j,k]+yy[j,k]*yy[j,k]), j, k)
r = np.sort(rr.flatten(),order='r')
r = np.flipud(r)

# write a file of r-values which can be read by convertP
system('rm -f rvaluesforP.txt')
fl = open('rvaluesforP.txt', 'w')
N = len(r)
fl.write('%d\n' % N)
for j in range(len(r)):
  fl.write('%.14f\n' % r[j]['r'])
fl.close()

# use convertP to get exact solution on in radial list
system('./convertP rvaluesforP.txt exactforP.txt')
system('rm rvaluesforP.txt')

# read back and put on grid
fl = open('exactforP.txt', 'r')
fl.readline()
for n in range(len(r)):
  tmp = fl.readline().split()
  j = r[n]['j']
  k = r[n]['k']
  h[j,k] = float(tmp[1])
  W[j,k] = float(tmp[4])
fl.close()
system('rm exactforP.txt')

system('rm -f inputforP.nc')
nc = PISMDataset("inputforP.nc", 'w')
nc.create_dimensions(x, x, time_dependent = True, use_time_bounds = True)

nc.define_2d_field("thk", time_dependent = False,
                   attrs = {"long_name"   : "ice thickness",
                            "units"       : "m",
                            "valid_range" : (0.0, 1.0e6)})
nc.define_2d_field("topg", time_dependent = False,
                   attrs = {"long_name"   : "bedrock topography",
                            "units"       : "m",
                            "valid_range" : (-1.0e6, 1.0e6)})
nc.define_2d_field("climatic_mass_balance", time_dependent = False,
                   attrs = {"long_name"   : "climatic mass balance for -surface given",
                            "units"       : "m year-1",
                            "valid_range" : (-1.0e6, 1.0e6)})
nc.define_2d_field("ice_surface_temp", time_dependent = False,
                   attrs = {"long_name"   : "ice surface temp (K) for -surface given",
                            "units"       : "Kelvin",
                            "valid_range" : (0.0, 1.0e6)})
nc.define_2d_field("bwat", time_dependent = False,
                   attrs = {"long_name"   : "thickness of basal water layer",
                            "units"       : "m",
                            "valid_range" : (0.0, 1.0e6)})

# fill with constants
nc.write("topg", 0.0*xx, time_dependent = False)
nc.write("climatic_mass_balance", 0.0*xx, time_dependent = False)
nc.write("ice_surface_temp", 260.0*np.ones(np.shape(xx)), time_dependent = False)

nc.write("thk", h, time_dependent = False)
nc.write("bwat", W, time_dependent = False)

#FIXME: need to write BWP

nc.close()

print "NetCDF file %s written" % "inputforP.nc"

cmd = "%s %s/pismr -boot_file inputforP.nc -Mx 51 -My 51 -Mz 11 -Lz 4000 -hydrology distributed -y 20 -no_mass -no_energy -o end.nc" % (mpiexec, pism_path)
stderr.write(cmd + '\n')
e = system(cmd)
if e != 0:
  exit(1)

#system("rm inputforP.nc")

exit(0)

