#!/usr/bin/env python

# FIXME:  for now you need   $ ln -s ../../util/PISMNC.py


#PLAN:  From this script generate list of r values for PISM map-plane grid of
#  given dims Mx x My.  Save this as ascii foo.txt.  Write short C program to make new
#  bar.txt containing exact solution to Test P at these r values; this C program
#  will call exactP_list() from src/verif/tests/exactTestP.[h|c].  Read this
#  bar.txt back in and complete the construction of a PISM input file start.nc
#  with fields x,y,thk,bwat,bwp.  Then run
#    pismr -i start.nc -hydrology distributed -y 20 -no_mass -no_energy -o end.nc  
#  Then compare
#    nccmp.py -v bwat,bwp end.nc start.nc difftestP.nc
#  and collect norms etc. for difftestP.nc.

from numpy import *
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

#stderr.write("Testing: temperature continuity at ice-bed interface (polythermal case).\n")
#cmd = "%s %s/pismv -test F -y 10 -verbose 1 -o bar-temp-continuity.nc" % (mpiexec, pism_path)
#stderr.write(cmd + '\n')
#e = system(cmd)
#if e != 0:
#  exit(1)

Lx = 25.0e3  # outside L = 22.5 km
Mx = 51      # 1 km grid

x = linspace(-Lx, Lx, Mx)
xx, yy = meshgrid(x, x)

nc = PISMDataset("foo.nc", 'w')
nc.create_dimensions(x, x, time_dependent = True, use_time_bounds = True)

nc.define_2d_field("bwat", time_dependent = False,
                   attrs = {"long_name"   : "thickness of basal water layer",
                            "comment"     : "for now a test variable",
                            "valid_range" : (0.0, 1.0e6)})

nc.write("bwat", (xx + Lx) / (2.0*Lx), time_dependent = False)  # silly test

nc.close()

#plan: get list rr which is sorted version of sqrt(xx*xx + yy*yy)
#      write this list as ascii file rvalues.txt
#      call ../../src/verif/tests/convertP to convert this to testPresults.txt
#      put this back onto cartesian grid and save a new .nc file which PISM can run
#      do "pismr -hydrology distributed" run
#      compare input and output .nc from this run using nccmp.py because values should be steady

#system("rm foo-temp-continuity.nc foo-temp-continuity.nc~ bar-temp-continuity.nc temp-temp-continuity.nc litho_temp-temp-continuity.nc")

exit(0)

