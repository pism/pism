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

m_schoof = 10;        # (pure number)
L_schoof = 40e3;      # meters
aspect_schoof = 0.05; # (pure)
H0_schoof = aspect_schoof * L_schoof; 
                      # = 2000 m THICKNESS
B_schoof = 3.7e8;     # Pa s^{1/3}; hardness 
                      # given on p. 239 of Schoof; why so big?
p_schoof = 4.0/3.0;   # = 1 + 1/n

class testi(PISM.ssa.SSAExactTestCase):

  def _initGrid(self):
    Mx = self.Mx; My = self.My;
    Ly = 3*L_schoof   # 300.0 km half-width (L=40.0km in Schoof's choice of variables)
    Lx = max(60.0e3, ((Mx - 1) / 2.) * (2.0 * Ly / (My - 1)) )
    PISM.util.init_shallow_grid(self.grid,Lx,Ly,Mx,My,PISM.NOT_PERIODIC);

  def _initPhysics(self):
    config = self.config
    config.set_flag("do_pseudo_plastic_till", False)
    basal = PISM.IceBasalResistancePlasticLaw(config)

    # irrelevant
    enthalpyconverter = PISM.EnthalpyConverter(config);

    config.set_string("ssa_flow_law", "isothermal_glen")
    config.set("ice_softness", pow(B_schoof, -config.get("Glen_exponent")))

    self.modeldata.setPhysics(basal,enthalpyconverter)


  def _initSSACoefficients(self):
    self._allocStdSSACoefficients()
    self._allocateBCs()
    vecs = self.modeldata.vecs

    vecs.bc_mask.set(0)
    vecs.thickness.set(H0_schoof)
    vecs.ice_mask.set(PISM.MASK_GROUNDED)

    # The finite difference code uses the following flag to treat 
    # the non-periodic grid correctly.
    self.config.set_flag("compute_surf_grad_inward_ssa", True);
    self.config.set("epsilon_ssa", 0.0);  # don't use this lower bound

    standard_gravity = self.config.get("standard_gravity");
    ice_rho = self.config.get("ice_density");
    theta = math.atan(0.001)
    f = ice_rho*standard_gravity*H0_schoof*math.tan(theta)
    grid = self.grid
    with PISM.util.Access(comm=[vecs.tauc]):
      for (i,j) in grid.points():
        y=grid.y[j]
        vecs.tauc[i,j] = f* (abs(y/L_schoof)**m_schoof)

    bc_mask = vecs.bc_mask
    vel_bc  = vecs.vel_bc
    surface = vecs.surface
    bed     = vecs.bed
    grid = self.grid
    with PISM.util.Access(comm=[surface,bed,vel_bc,bc_mask]):
      for (i,j) in grid.points():
        x=grid.x[i]; y=grid.y[j]
        (bed_ij,junk,u,v) = PISM.exactI(m_schoof,x,y)
        bed[i,j] = bed_ij
        surface[i,j] = bed_ij + H0_schoof

        edge = ( (j == 0) or (j == grid.My - 1) ) or ( (i==0) or (i==grid.Mx-1) );
        if (edge):
          bc_mask[i,j] = 1;
          vel_bc[i,j].u = u;
          vel_bc[i,j].v = v;
      
  def exactSolution( self, i, j, x, y ):
    (j1,j2,u,v) = PISM.exactI(m_schoof,x,y)
    return [u,v]

# The main code for a run follows:
if __name__ == '__main__':
  context = PISM.Context()

  PISM.set_abort_on_sigint(True)

  for o in PISM.OptionsGroup(context.com,"","Test J"):
    Mx = PISM.optionsInt("-Mx","Number of grid points in x-direction",default=11)
    My = PISM.optionsInt("-My","Number of grid points in y-direction",default=61)
    output_file = PISM.optionsString("-o","output file",default="testi.nc")
    verbosity = PISM.optionsInt("-verbose","verbosity level",default=3)

  PISM.setVerbosityLevel(verbosity)

  tc = testi(Mx,My)
  tc.run(output_file)



