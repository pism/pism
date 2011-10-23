#! /usr/bin/env python
#
# Copyright (C) 2011 Ed Bueler and Constantine Khroulev and David Maxwell
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


import sys, petsc4py
petsc4py.init(sys.argv)
from petsc4py import PETSc

import PISM, math

L=50.e3; #// 50km half-width
H0=500; #// m
dhdx = 0.005; #// pure number, slope of surface & bed
nu0 = 30.0 * 1.0e6 * PISM.secpera; #/* = 9.45e14 Pa s */
tauc0 = 1.e4; #// 1kPa

class test_linear(PISM.ssa.SSAExactTestCase):
  def _initGrid(self):
    PISM.util.init_shallow_grid(self.grid,L,L,self.Mx,self.My,PISM.NONE);

  def _initPhysics(self):
    config = self.config
    linear_q = 1.
    basal = PISM.IceBasalResistancePlasticLaw(
           config.get("plastic_regularization","1/year","1/second"),
           True, # do not force pure plastic
           linear_q,
           config.get("pseudo_plastic_uthreshold","1/year","1/second"));

    enthalpyconverter = PISM.EnthalpyConverter(config)
    ice = PISM.CustomGlenIce(self.grid.com, "", config, enthalpyconverter)
    self.solver.setPhysics(ice,basal,enthalpyconverter)

  def _initSSACoefficients(self):
    solver = self.solver
    solver.allocateCoeffs()
    
    solver.thickness.set(H0)
    solver.surface.set(H0)
    solver.bed.set(0.)
    solver.tauc.set(tauc0)

    # Set boundary conditions (Dirichlet all the way around).
    self.solver.allocateBCs()
    
    bc_mask = self.solver.bc_mask
    bc_mask.set(PISM.MASK_GROUNDED)

    vel_bc  = self.solver.vel_bc
    
    grid = self.grid
    vars = [bc_mask, vel_bc]
    with PISM.util.Access(vars):
      for (i,j) in grid.points():
        edge = ( (j == 0) or (j == grid.My - 1) ) or ( (i==0) or (i==grid.Mx-1) );
        if edge:
          bc_mask[i,j] = 1;
          x = grid.x[i]; y=grid.y[j];
          [u,v] = self.exactSolution(i,j,x,y);
          vel_bc[i,j].u = u;
          vel_bc[i,j].v = v;
      for v in vars:
        v.beginGhostComm(); v.endGhostComm();

  def _initSSA(self):
    # The following ensure that the strength extension is used everywhere
    se = self.solver.ssa.strength_extension
    se.set_notional_strength(nu0 * H0);
    se.set_min_thickness(4000*10);

    # For the benefit of SSAFD on a non-periodic grid
    self.config.set_flag("compute_surf_grad_inward_ssa", True);


  def exactSolution(self, i, j, x, y ):
    tauc_threshold_velocity = self.config.get("pseudo_plastic_uthreshold") / PISM.secpera;
    v0 = 100./PISM.secpera ; #// 100 m/s.
    alpha = math.sqrt( (tauc0/tauc_threshold_velocity) / (4*nu0*H0) );
    return [v0*math.exp( -alpha*(x-L)), 0]

# The main code for a run follows:
if __name__ == '__main__':
  context = PISM.Context()

  PISM.set_abort_on_sigint(True)

  for o in PISM.OptionsGroup(context.com,"","Test Lineaer"):
    Mx = PISM.optionsInt("-Mx","Number of grid points in x-direction",default=61)
    My = PISM.optionsInt("-My","Number of grid points in y-direction",default=61)
    output_file = PISM.optionsString("-o","output file",default="test_linear.nc")
    verbosity = PISM.optionsInt("-verbose","verbosity level",default=3)

  PISM.setVerbosityLevel(verbosity)
  tc = test_linear(Mx,My)
  tc.run(output_file)
