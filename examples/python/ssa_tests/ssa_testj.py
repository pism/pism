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

class testj(PISM.ssa.SSAExactTestCase):

  def _initGrid(self):
    halfWidth = 300.0e3
    Lx = halfWidth; Ly = halfWidth
    PISM.util.init_shallow_grid(self.grid,Lx,Ly,self.Mx,self.My,PISM.XY_PERIODIC);

  def _initPhysics(self):
    config = self.modeldata.config
    basal = PISM.IceBasalResistancePlasticLaw(config)

    enthalpyconverter = PISM.EnthalpyConverter(config)

    config.set_string("ssa_flow_law", "isothermal_glen")

    self.modeldata.setPhysics(basal,enthalpyconverter)


  def _initSSACoefficients(self):
    self._allocStdSSACoefficients()
    self._allocateBCs()

    vecs = self.modeldata.vecs
    vecs.tauc.set(0.0) # irrelevant for test J
    vecs.bed.set(0.0) 
    vecs.ice_mask.set(PISM.MASK_FLOATING)
    vecs.bc_mask.set(0) # No dirichlet data.

    vecs.enthalpy.set(528668.35); # arbitrary; corresponds to 263.15 Kelvin at depth=0.

    ocean_rho = self.config.get("sea_water_density");
    ice_rho = self.config.get("ice_density");
    
    # The PISM.utils.Access object ensures that we call beginAccess for each
    # variable in 'vars', and that endAccess is called for each one on exiting
    # the 'with' block.
    
    with PISM.util.Access(comm=[vecs.thickness, vecs.surface, vecs.bc_mask, vecs.vel_bc]):
      grid = self.grid
      for (i,j) in grid.points():
        x = grid.x[i]; y = grid.y[j]
        (H,junk,u,v) = PISM.exactJ(x,y);
        vecs.thickness[i,j] = H;
        vecs.surface[i,j] = (1.0 - ice_rho / ocean_rho) * H; #// FIXME task #7297
  
        # // special case at center point (Dirichlet BC)
        if (i == (grid.Mx)/2) and (j == (grid.My)/2):
          vecs.bc_mask[i,j] = 1;
          vecs.vel_bc[i,j] = [u,v]

  def _initSSA(self):
    # Test J has a viscosity that is independent of velocity.  So we force a 
    # constant viscosity by settting the strength_extension
    # thickness larger than the given ice thickness. (max = 770m).

    nu0 = 30.0 * 1.0e6 * PISM.secpera; # 9.45e14 Pa s 
    H0 = 500.0;                        # 500 m typical thickness

    ssa = self.ssa
    ssa.strength_extension.set_notional_strength(nu0 * H0);
    ssa.strength_extension.set_min_thickness(800.);

  def exactSolution(self,i,j,x,y):
    (j1,j2,u,v) = PISM.exactJ(x,y);
    return [u,v]


# The main code for a run follows:
if __name__ == '__main__':
  context = PISM.Context()

  for o in PISM.OptionsGroup(context.com,"","Test J"):
    Mx = PISM.optionsInt("-Mx","Number of grid points in x-direction",default=61)
    My = PISM.optionsInt("-My","Number of grid points in y-direction",default=61)
    output_file = PISM.optionsString("-o","output file",default="testj.nc")
    verbosity = PISM.optionsInt("-verbose","verbosity level",default=3)

  PISM.setVerbosityLevel(verbosity)
  tc = testj(Mx,My)
  tc.run(output_file)
