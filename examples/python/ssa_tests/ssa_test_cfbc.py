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

context = PISM.Context()

H0 = 600.          # meters
V0 = convert(300, "m/year", "m/second")
C = 2.45e-18     # "typical constant ice parameter"
T = 400          # time used to compute the calving front location

Q0 = V0 * H0
Hc1 = 4. * C / Q0
Hc2 = 1. / (H0 ** 4)


def H_exact(x):
    return (Hc1 * x + Hc2) ** (-1 / 4.)


def u_exact(x):
    return Q0 / H_exact(x)


class test_cfbc(PISM.ssa.SSAExactTestCase):

    def _initGrid(self):
        self.grid = None
        halfWidth = 250.0e3  # 500.0 km length
        Lx = halfWidth
        Ly = halfWidth
        self.grid = PISM.IceGrid.Shallow(PISM.Context().ctx, Lx, Ly, 0, 0,
                                         self.Mx, self.My,
                                         PISM.CELL_CENTER,
                                         PISM.Y_PERIODIC)

    def _initPhysics(self):
        config = self.config

        config.set_double("flow_law.isothermal_Glen.ice_softness",
                          pow(1.9e8, -config.get_double("stress_balance.ssa.Glen_exponent")))
        config.set_boolean("stress_balance.ssa.compute_surface_gradient_inward", False)
        config.set_boolean("stress_balance.calving_front_stress_bc", True)
        config.set_string("stress_balance.ssa.flow_law", "isothermal_glen")

        enthalpyconverter = PISM.EnthalpyConverter(config)

        self.modeldata.setPhysics(enthalpyconverter)

    def _initSSACoefficients(self):
        self._allocStdSSACoefficients()
        self._allocateBCs()

        vecs = self.modeldata.vecs

        vecs.tauc.set(0.0)     # irrelevant
        vecs.bedrock_altitude.set(-1000.0)  # assures shelf is floating

        EC = PISM.EnthalpyConverter(PISM.Context().config)
        enth0 = EC.enthalpy(273.15, 0.01, 0)  # 0.01 water fraction
        vecs.enthalpy.set(enth0)

        grid = self.grid
        thickness = vecs.land_ice_thickness
        surface = vecs.surface_altitude
        bc_mask = vecs.bc_mask
        vel_bc = vecs.vel_bc
        ice_mask = vecs.mask

        ocean_rho = self.config.get_double("constants.sea_water.density")
        ice_rho = self.config.get_double("constants.ice.density")

        x_min = grid.x(0)
        with PISM.vec.Access(comm=[thickness, surface, bc_mask, vel_bc, ice_mask]):
            for (i, j) in grid.points():
                x = grid.x(i)
                if i != grid.Mx() - 1:
                    thickness[i, j] = H_exact(x - x_min)
                    ice_mask[i, j] = PISM.MASK_FLOATING
                else:
                    thickness[i, j] = 0
                    ice_mask[i, j] = PISM.MASK_ICE_FREE_OCEAN

                surface[i, j] = (1.0 - ice_rho / ocean_rho) * thickness[i, j]

                if i == 0:
                    bc_mask[i, j] = 1
                    vel_bc[i, j].u = V0
                    vel_bc[i, j].v = 0.
                else:
                    bc_mask[i, j] = 0
                    vel_bc[i, j].u = 0.
                    vel_bc[i, j].v = 0.

    def exactSolution(self, i, j, x, y):
        x_min = self.grid.x(0)
        if i != self.grid.Mx() - 1:
            u = u_exact(x - x_min)
        else:
            u = 0
        return [u, 0]


if __name__ == '__main__':

    config = PISM.Context().config

    tc = test_cfbc(int(config.get_double("grid.Mx")),
                   int(config.get_double("grid.My")))

    tc.run(config.get_string("output.file_name"))
