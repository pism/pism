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
    config.set_flag("do_pseudo_plastic_till", True)
    config.set("pseudo_plastic_q", 1.0)
    basal = PISM.IceBasalResistancePseudoPlasticLaw(config)

    enthalpyconverter = PISM.EnthalpyConverter(config)

    config.set_string("ssa_flow_law", "isothermal_glen")

    self.modeldata.setPhysics(basal,enthalpyconverter)

  def _initSSACoefficients(self):
    self._allocStdSSACoefficients()
    self._allocateBCs()

    vecs = self.modeldata.vecs
    
    vecs.thickness.set(H0)
    vecs.surface.set(H0)
    vecs.bed.set(0.)
    vecs.tauc.set(tauc0)

    bc_mask = vecs.bc_mask
    bc_mask.set(0)

    vel_bc  = vecs.vel_bc
    
    grid = self.grid
    with PISM.util.Access(comm=[bc_mask, vel_bc]):
      for (i,j) in grid.points():
        edge = ( (j == 0) or (j == grid.My - 1) ) or ( (i==0) or (i==grid.Mx-1) );
        if edge:
          bc_mask[i,j] = 1;
          x = grid.x[i]; y=grid.y[j];
          [u,v] = self.exactSolution(i,j,x,y);
          vel_bc[i,j].u = u;
          vel_bc[i,j].v = v;

  def _initSSA(self):
    # The following ensure that the strength extension is used everywhere
    se = self.ssa.strength_extension
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
