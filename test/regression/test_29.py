#!/usr/bin/env python

# FIXME:  for now need
#   $ ln -s ../../util/PISMNC.py
#   $ (cd ../../src/verif/tests/ && make convertP)
#   $ ln -s ../../src/verif/tests/convertP

import numpy as np
from sys import exit, argv, stderr
from os import system
from PISMNC import PISMDataset

# example:  ./test_29.py ../../build "mpiexec -n 4" 101

# generate config file
system("ncgen -o test29.nc test29.cdl")

if len(argv)<2:
  pism_path="."
else:
  pism_path=argv[1]

if len(argv)<3:
  mpiexec=""
else:
  mpiexec=argv[2]

if len(argv)<4:
  Mx = 51     # generates 1 km grid
else:
  Mx = int(argv[3])

stderr.write("Testing: Test P verification of '-hydrology distributed'.\n")

Lx = 25.0e3  # outside L = 22.5 km
Phi0 = 0.20  # 20 cm a-1 basal melt rate

x = np.linspace(-Lx, Lx, Mx)
xx, yy = np.meshgrid(x, x)

h = np.zeros(np.shape(xx))
W = np.zeros(np.shape(xx))
P = np.zeros(np.shape(xx))

def radially_outward(mag, x, y):
  """return components of a vector field  V(x,y)  which is radially-outward from
the origin and has magnitude mag"""
  r = np.sqrt(x*x + y*y)
  if r == 0.0:
    return (0.0, 0.0)
  vx = mag * x / r
  vy = mag * y / r
  return (vx, vy)

bcflag = np.ones(np.shape(xx))
magvb = np.zeros(np.shape(xx))
ussa = np.zeros(np.shape(xx))
vssa = np.zeros(np.shape(xx))

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
  h[j,k] = float(tmp[1])     # ice thickness in m
  magvb[j,k] = float(tmp[2]) # sliding speed in m s-1
  ussa[j,k], vssa[j,k] = radially_outward(magvb[j,k],xx[j,k],yy[j,k])
  W[j,k] = float(tmp[4])     # water thickness in m
  P[j,k] = float(tmp[5])     # water pressure in Pa
fl.close()
system('rm exactforP.txt')

system('rm -f inputforP.nc')
nc = PISMDataset("inputforP.nc", 'w')
nc.create_dimensions(x, x, time_dependent = True, use_time_bounds = True)

nc.define_2d_field("thk", time_dependent = False,
                   attrs = {"long_name"   : "ice thickness",
                            "units"       : "m",
                            "valid_min"   : 0.0,
                            "standard_name" : "land_ice_thickness"})
nc.define_2d_field("topg", time_dependent = False,
                   attrs = {"long_name"   : "bedrock topography",
                            "units"       : "m",
                            "standard_name" : "bedrock_altitude"})
nc.define_2d_field("climatic_mass_balance", time_dependent = False,
                   attrs = {"long_name"   : "climatic mass balance for -surface given",
                            "units"       : "m year-1",
                            "standard_name" : "land_ice_surface_specific_mass_balance"})
nc.define_2d_field("ice_surface_temp", time_dependent = False,
                   attrs = {"long_name"   : "ice surface temp (K) for -surface given",
                            "units"       : "Kelvin",
                            "valid_min"   : 0.0})
nc.define_2d_field("bmelt", time_dependent = False,
                   attrs = {"long_name"   : "basal melt rate",
                            "units"       : "m year-1",
                            "standard_name" : "land_ice_basal_melt_rate"})

nc.define_2d_field("bwat", time_dependent = False,
                   attrs = {"long_name"   : "thickness of basal water layer",
                            "units"       : "m",
                            "valid_min"   : 0.0})
nc.define_2d_field("bwp", time_dependent = False,
                   attrs = {"long_name"   : "water pressure in basal water layer",
                            "units"       : "Pa",
                            "valid_min"   : 0.0})

nc.define_2d_field("bcflag", time_dependent = False,
                   attrs = {"long_name"   : "if =1, apply u_ssa_bc and v_ssa_bc as sliding velocity"})
nc.define_2d_field("u_ssa_bc", time_dependent = False,
                   attrs = {"long_name"   : "x-component of prescribed sliding velocity",
                            "units"       : "m s-1"})
nc.define_2d_field("v_ssa_bc", time_dependent = False,
                   attrs = {"long_name"   : "y-component of prescribed sliding velocity",
                            "units"       : "m s-1"})

# fill with constants
nc.write("topg", 0.0*xx, time_dependent = False)
nc.write("climatic_mass_balance", 0.0*xx, time_dependent = False)
nc.write("ice_surface_temp", 260.0 * np.ones(np.shape(xx)), time_dependent = False)
nc.write("bmelt", Phi0 * np.ones(np.shape(xx)), time_dependent = False)

nc.write("thk", h, time_dependent = False)
nc.write("bwat", W, time_dependent = False)
nc.write("bwp", P, time_dependent = False)

nc.write("bcflag", bcflag, time_dependent = False)
nc.write("u_ssa_bc", ussa, time_dependent = False)
nc.write("v_ssa_bc", vssa, time_dependent = False)

nc.close()

print "NetCDF file %s written" % "inputforP.nc"

cmd = "%s %s/pismr -config_override test29.nc -boot_file inputforP.nc -Mx %d -My %d -Mz 11 -Lz 4000 -hydrology distributed -report_mass_accounting -y 0.08333333333333 -max_dt 0.01 -no_mass -no_energy -ssa_sliding -ssa_dirichlet_bc -o end.nc" % (mpiexec, pism_path, Mx, Mx)

stderr.write(cmd + '\n')
e = system(cmd)
if e != 0:
  exit(1)

system("rm -f diffP.nc")

#FIXME this ain't elegant
system("ncdiff -v bwat,bwp end.nc inputforP.nc diffP.nc")
system(r"ncap2 -O -s'errbwat=avg(abs(bwat));errbwp=avg(abs(bwp))' diffP.nc diffP.nc")
system(r"ncdump -v errbwat diffP.nc |grep 'errbwat ='")
system(r"ncdump -v errbwp  diffP.nc |grep 'errbwp ='")

#cleanup:
#system("rm test29.nc inputforP.nc end.nc foo.txt diffP.nc")

exit(0)

