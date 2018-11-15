#! /usr/bin/env python
#
# Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016, 2018 Ed Bueler and Constantine Khroulev and David Maxwell
#
# This file is part of PISM.
#
# PISM is free software; you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation; either version 3 of the License, or (at your option) any later
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
import math

m_schoof = 10        # (pure number)
L_schoof = 40e3      # meters
aspect_schoof = 0.05  # (pure)
H0_schoof = aspect_schoof * L_schoof
# = 2000 m THICKNESS
B_schoof = 3.7e8     # Pa s^{1/3}; hardness
# given on p. 239 of Schoof; why so big?
p_schoof = 4.0 / 3.0   # = 1 + 1/n


class testi(PISM.ssa.SSAExactTestCase):

    def _initGrid(self):
        Mx = self.Mx
        My = self.My
        Ly = 3 * L_schoof   # 300.0 km half-width (L=40.0km in Schoof's choice of variables)
        Lx = max(60.0e3, ((Mx - 1) / 2.) * (2.0 * Ly / (My - 1)))
        self.grid = PISM.IceGrid.Shallow(PISM.Context().ctx, Lx, Ly, 0, 0,
                                         Mx, My,
                                         PISM.CELL_CORNER,
                                         PISM.NOT_PERIODIC)

    def _initPhysics(self):
        config = self.config
        config.set_boolean("basal_resistance.pseudo_plastic.enabled", False)

        # irrelevant
        enthalpyconverter = PISM.EnthalpyConverter(config)

        config.set_string("stress_balance.ssa.flow_law", "isothermal_glen")
        config.set_double("flow_law.isothermal_Glen.ice_softness", pow(
            B_schoof, -config.get_double("stress_balance.ssa.Glen_exponent")))

        self.modeldata.setPhysics(enthalpyconverter)

    def _initSSACoefficients(self):
        self._allocStdSSACoefficients()
        self._allocateBCs()
        vecs = self.modeldata.vecs

        vecs.bc_mask.set(0)
        vecs.thk.set(H0_schoof)
        vecs.mask.set(PISM.MASK_GROUNDED)

        # The finite difference code uses the following flag to treat
        # the non-periodic grid correctly.
        self.config.set_boolean("stress_balance.ssa.compute_surface_gradient_inward", True)
        self.config.set_double("stress_balance.ssa.epsilon", 0.0)  # don't use this lower bound

        standard_gravity = self.config.get_double("constants.standard_gravity")
        ice_rho = self.config.get_double("constants.ice.density")
        theta = math.atan(0.001)
        f = ice_rho * standard_gravity * H0_schoof * math.tan(theta)
        grid = self.grid
        with PISM.vec.Access(comm=[vecs.tauc]):
            for (i, j) in grid.points():
                y = grid.y(j)
                vecs.tauc[i, j] = f * (abs(y / L_schoof) ** m_schoof)

        bc_mask = vecs.bc_mask
        vel_bc = vecs.vel_bc
        surface = vecs.surface_altitude
        bed = vecs.bedrock_altitude
        grid = self.grid
        with PISM.vec.Access(comm=[surface, bed, vel_bc, bc_mask]):
            for (i, j) in grid.points():
                p = PISM.exactI(m_schoof, grid.x(i), grid.y(j))
                bed[i, j] = p.bed
                surface[i, j] = p.bed + H0_schoof

                edge = ((j == 0) or (j == grid.My() - 1)) or ((i == 0) or (i == grid.Mx() - 1))
                if edge:
                    bc_mask[i, j] = 1
                    vel_bc[i, j].u = p.u
                    vel_bc[i, j].v = p.v

    def exactSolution(self, i, j, x, y):
        p = PISM.exactI(m_schoof, x, y)
        return [p.u, p.v]


# The main code for a run follows:
if __name__ == '__main__':
    context = PISM.Context()
    config = context.config

    PISM.set_abort_on_sigint(True)

    config.set_double("grid.Mx", 11, PISM.CONFIG_DEFAULT)
    config.set_double("grid.My", 61, PISM.CONFIG_DEFAULT)

    tc = testi(int(config.get_double("grid.Mx")),
               int(config.get_double("grid.My")))
    tc.run(config.get_string("output.file_name"))
