#!/usr/bin/env python

# FIXME:  for now need
#   $ ln -s ../../util/PISMNC.py
#   $ (cd ../../src/verif/tests/ && make convertP)
#   $ ln -s ../../src/verif/tests/convertP

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
P = np.zeros(np.shape(xx))

magvb = np.zeros(np.shape(xx))
bcflag = np.ones(np.shape(xx))
u_ssa_bc = np.zeros(np.shape(xx))
v_ssa_bc = np.zeros(np.shape(xx))

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
                            "valid_min"   : 0.0})
nc.define_2d_field("topg", time_dependent = False,
                   attrs = {"long_name"   : "bedrock topography",
                            "units"       : "m"})
nc.define_2d_field("climatic_mass_balance", time_dependent = False,
                   attrs = {"long_name"   : "climatic mass balance for -surface given",
                            "units"       : "m year-1"})
nc.define_2d_field("ice_surface_temp", time_dependent = False,
                   attrs = {"long_name"   : "ice surface temp (K) for -surface given",
                            "units"       : "Kelvin",
                            "valid_min"   : 0.0})
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
nc.write("ice_surface_temp", 260.0*np.ones(np.shape(xx)), time_dependent = False)

nc.write("thk", h, time_dependent = False)
nc.write("bwat", W, time_dependent = False)
nc.write("bwp", P, time_dependent = False)

nc.write("bcflag", bcflag, time_dependent = False)
FIXME: u_ssa_bc, v_ssa_bc

#FIXME  need to write velocity components for SSA boundary conditions
#       following examples/ross/

nc.close()

print "NetCDF file %s written" % "inputforP.nc"

cmd = "%s %s/pismr -boot_file inputforP.nc -Mx 51 -My 51 -Mz 11 -Lz 4000 -hydrology distributed -y 1.0 -max_dt 0.1 -no_mass -no_energy -ssa_sliding -ssa_dirichlet_bc -o end.nc" % (mpiexec, pism_path)
# FIXME probably need ssabc ish option

stderr.write(cmd + '\n')
e = system(cmd)
if e != 0:
  exit(1)

#system("rm inputforP.nc")

#FIXME:  need to compare
#    nccmp.py -v bwat,bwp end.nc inputforP.nc difftestP.nc
#  and collect norms etc. for difftestP.nc.

exit(0)

