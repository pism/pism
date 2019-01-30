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
from PISM.util import convert

class testj(PISM.ssa.SSAExactTestCase):

    def _initGrid(self):
        halfWidth = 300.0e3
        Lx = halfWidth
        Ly = halfWidth
        ctx = PISM.Context().ctx
        self.grid = PISM.IceGrid.Shallow(ctx, Lx, Ly, 0, 0,
                                         self.Mx, self.My,
                                         PISM.CELL_CENTER,
                                         PISM.XY_PERIODIC)

    def _initPhysics(self):
        config = self.modeldata.config
        config.set_boolean("basal_resistance.pseudo_plastic.enabled", False)

        enthalpyconverter = PISM.EnthalpyConverter(config)

        config.set_string("stress_balance.ssa.flow_law", "isothermal_glen")

        self.modeldata.setPhysics(enthalpyconverter)

    def _initSSACoefficients(self):
        self._allocStdSSACoefficients()
        self._allocateBCs()

        vecs = self.modeldata.vecs

        vecs.tauc.set(0.0)  # irrelevant for test J
        # ensures that the ice is floating (max. thickness if 770 m)
        vecs.bedrock_altitude.set(-1000.0)
        vecs.mask.set(PISM.MASK_FLOATING)
        vecs.bc_mask.set(0)  # No dirichlet data.

        EC = PISM.EnthalpyConverter(PISM.Context().config)
        enth0 = EC.enthalpy(273.15, 0.01, 0)  # 0.01 water fraction
        vecs.enthalpy.set(enth0)

        ocean_rho = self.config.get_double("constants.sea_water.density")
        ice_rho = self.config.get_double("constants.ice.density")

        # The PISM.vec.Access object ensures that we call beginAccess for each
        # variable in 'vars', and that endAccess is called for each one on exiting
        # the 'with' block.

        with PISM.vec.Access(comm=[vecs.land_ice_thickness,
                                   vecs.surface_altitude,
                                   vecs.bc_mask,
                                   vecs.vel_bc]):
            grid = self.grid
            for (i, j) in grid.points():
                p = PISM.exactJ(grid.x(i), grid.y(j))
                vecs.land_ice_thickness[i, j] = p.H
                vecs.surface_altitude[i, j] = (1.0 - ice_rho / ocean_rho) * p.H  # // FIXME task #7297

                # special case at center point (Dirichlet BC)
                if (i == grid.Mx() // 2) and (j == grid.My() // 2):
                    vecs.bc_mask[i, j] = 1
                    vecs.vel_bc[i, j] = [p.u, p.v]

    def _initSSA(self):
        # Test J has a viscosity that is independent of velocity.  So we force a
        # constant viscosity by settting the strength_extension
        # thickness larger than the given ice thickness. (max = 770m).

        nu0 = convert(30.0, "MPa year", "Pa s")
        H0 = 500.0              # 500 m typical thickness

        ssa = self.ssa
        ssa.strength_extension.set_notional_strength(nu0 * H0)
        ssa.strength_extension.set_min_thickness(800.)

    def exactSolution(self, i, j, x, y):
        p = PISM.exactJ(x, y)
        return [p.u, p.v]


# The main code for a run follows:
if __name__ == '__main__':
    context = PISM.Context()
    config = context.config

    tc = testj(int(config.get_double("grid.Mx")), int(config.get_double("grid.My")))
    tc.run(config.get_string("output.file_name"))
