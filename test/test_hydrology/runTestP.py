#!/usr/bin/env python

# high-res and parallel example:
#    ./runTestP.py ../../buildbwp "mpiexec -n 4" 201
# example which should suffice for regression:
#    ./runTestP.py ../../buildbwp "" 21

from os import system
system("ln -sf ../../util/PISMNC.py")

import numpy as np
from sys import exit, argv, stderr
from PISMNC import PISMDataset
from exactP import exactP_list

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

# generate config file
cdlcontent = """netcdf pism_overrides {
    variables:
    byte pism_overrides;
    pism_overrides:ice_softness = 3.1689e-24;
    pism_overrides:ice_softness_doc = "Pa-3 s-1; ice softness; NOT DEFAULT";
    pism_overrides:hydrology_hydraulic_conductivity = 1.0e-2;
    pism_overrides:hydrology_hydraulic_conductivity_doc = "m s-1; = K; NOT DEFAULT";
    pism_overrides:hydrology_hydraulic_conductivity_at_large_W = 1.0e-2;
    pism_overrides:hydrology_hydraulic_conductivity_doc = "m s-1; = K";
}"""
cdlf = file("testP.cdl","w")
cdlf.write(cdlcontent)
cdlf.close()
system("ncgen -o testP.nc testP.cdl")

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

EPS_ABS = 1.0e-12
EPS_REL = 1.0e-15
h_r, magvb_r, Wcrit_r, W_r, P_r = exactP_list(np.array(r[:]['r']), EPS_ABS, EPS_REL, 1)

# put on grid
for n in range(len(r)):
  j = r[n]['j']
  k = r[n]['k']
  h[j,k] = h_r[n]         # ice thickness in m
  magvb[j,k] = magvb_r[n] # sliding speed in m s-1
  ussa[j,k], vssa[j,k] = radially_outward(magvb[j,k],xx[j,k],yy[j,k])
  W[j,k] = W_r[n]         # water thickness in m
  P[j,k] = P_r[n]         # water pressure in Pa

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

cmd = "%s %s/pismr -config_override testP.nc -boot_file inputforP.nc -Mx %d -My %d -Mz 11 -Lz 4000 -hydrology distributed -report_mass_accounting -y 0.08333333333333 -max_dt 0.01 -no_mass -no_energy -ssa_sliding -ssa_dirichlet_bc -o end.nc" % (mpiexec, pism_path, Mx, Mx)

stderr.write(cmd + '\n')
e = system(cmd)
if e != 0:
  exit(1)

#evaluate error using NCO
system("ncdiff -O -v bwat,bwp end.nc inputforP.nc diffP.nc")
system(r"ncap2 -O -s'averrbwat=avg(abs(bwat));averrbwp=avg(abs(bwp));maxerrbwat=max(abs(bwat));maxerrbwp=max(abs(bwp));' diffP.nc diffP.nc")
system(r"ncdump -v averrbwat,maxerrbwat diffP.nc |grep 'errbwat ='")
system(r"ncdump -v averrbwp,maxerrbwp diffP.nc |grep 'errbwp ='")

#cleanup:
system("rm testP.cdl testP.nc inputforP.nc end.nc diffP.nc PISMNC.py*")

exit(0)

