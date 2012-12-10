#! /usr/bin/env python
#
# Copyright (C) 2011, 2012 Ed Bueler and Constantine Khroulev and David Maxwell
# 
# This file is part of PISM.
# 
# PISM is free software; you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation; either version 2 of the License, or (at your option) any later
# version.
# 
# PISM is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.
# 
# You should have received a copy of the GNU General Public License
# along with PISM; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA


import PISM

help = \
"""
SSA_TESTCFBC
  Testing program for PISM's implementations of the SSA.
  Does a time-independent calculation.  Does not run IceModel or a derived
  class thereof. Uses the van der Veen flow-line shelf geometry. Also may be
  used in a PISM software (regression) test.
"""

usage = \
"""
usage of SSA_TEST_CFBC:
  run ssa_test_cfbc -Mx <number> -My <number>
"""

H0 = 600.;          # meters
V0 = 300./PISM.secpera;  # 300 meters/year
C  = 2.45e-18;     # "typical constant ice parameter"
T  = 400;          # time used to compute the calving front location

Q0 = V0*H0;
Hc1 = 4.*C/Q0
Hc2 = 1./(H0**4)
def H_exact(x):
  return (Hc1*x+Hc2)**(-1/4.)
  
def u_exact(x):
  return Q0 / H_exact(x)

class test_cfbc(PISM.ssa.SSAExactTestCase):

  def _initGrid(self):
    self.grid = PISM.Context().newgrid()
    halfWidth = 250.0e3;  # 500.0 km length
    Lx = halfWidth; Ly = halfWidth;
    PISM.util.init_shallow_grid(self.grid,Lx,Ly,self.Mx,self.My,PISM.Y_PERIODIC)

  def _initPhysics(self):
    config = self.config
    config.set_flag("compute_surf_grad_inward_ssa", True);
    config.set_flag("calving_front_stress_boundary_condition", True);

    basal = PISM.IceBasalResistancePlasticLaw(config)

    enthalpyconverter = PISM.EnthalpyConverter(config);

    config.set_string("ssa_flow_law", "isothermal_glen")
    config.set("ice_softness", pow(1.9e8, -config.get("Glen_exponent")))
    
    self.modeldata.setPhysics(ice,basal,enthalpyconverter)

  def _initSSACoefficients(self):
    self._allocStdSSACoefficients()
    self._allocateBCs()
    vecs = self.modeldata.vecs

    vecs.tauc.set(0.0)     # irrelevant
    vecs.bed.set(-1000.0); # assures shelf is floating
    vecs.enthalpy.set(528668.35); # arbitrary; corresponds to 263.15 Kelvin at depth=0.

    grid = self.grid
    thickness=vecs.thickness;
    surface = vecs.surface
    bc_mask = vecs.bc_mask
    vel_bc = vecs.vel_bc
    ice_mask = vecs.ice_mask

    ocean_rho = self.config.get("sea_water_density");
    ice_rho   = self.modeldata.ice.rho
    
    with PISM.util.Access(comm=[thickness,surface,bc_mask,vel_bc,ice_mask]):
      for (i,j) in grid.points():
        x = grid.x[i]
        if x <= 0:
          thickness[i,j] = H_exact( x + self.grid.Lx )
          ice_mask[i,j]  = PISM.MASK_FLOATING
        else:
          thickness[i,j] = 0
          ice_mask[i,j]  = PISM.MASK_ICE_FREE_OCEAN
      
        surface[i,j] = (1.0 - ice_rho/ocean_rho)*thickness[i,j]

        if i==0:
          bc_mask[i,j] = 1
          vel_bc[i,j].u = V0
          vel_bc[i,j].v = 0.
        else:
          bc_mask[i,j] = 0
          vel_bc[i,j].u = 0.
          vel_bc[i,j].v = 0.


  def exactSolution(self,i,j,x,y):
    if x<= 0:
      u = u_exact(x+self.grid.Lx)
    else:
      u = 0
    return [u,0]

if __name__ == '__main__':
  context = PISM.Context()

  # if PISM.optionsSet('-usage') or PISM.optionsSet('-help'):
  #   PISM.verbPrintf(1,context.com,help)
  #   PISM.verbPrintf(1,context.com,usage)

  for o in PISM.OptionsGroup(context.com,"","Test CFBC"):
    Mx = PISM.optionsInt("-Mx","Number of grid points in x-direction",default=61)
    My = PISM.optionsInt("-My","Number of grid points in y-direction",default=61)
    output_file = PISM.optionsString("-o","output file",default="ssa_test_cfbc.nc")
    verbosity = PISM.optionsInt("-verbose","verbosity level",default=3)

  PISM.setVerbosityLevel(verbosity)

  context.config.set_string('ssa_method','fd')

  tc = test_cfbc(Mx,My)
  tc.run(output_file)
