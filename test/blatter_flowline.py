#!/usr/bin/env python3
"""This script simplifies flow line tests using the Blatter-Pattyn solver.
"""

import numpy as np
import PISM

ctx = PISM.Context()
config = ctx.config

config.set_string("stress_balance.blatter.flow_law", "isothermal_glen")

def H_Halfar(H_max, R_max, r):
    SperA = 31556926.0
    n = 3.0
    H0 = 3600.0
    R0 = 750000.0
    t = 450 * SperA

    # alpha=(2-(n+1)*lambda)/(5*n+3)=1/9
    alpha = 1.0 / 9.0
    # beta=(1+(2*n+1)*lambda)/(5*n+3)=1/18
    beta = 1.0 / 18.0
    # t0 = (beta/Gamma) * pow((2n+1)/(n+1),n)*(pow(R0,n+1)/pow(H0,2n+1))
    t0 = 422.45 * SperA

    Rmargin = R0 * pow(t / t0, beta);

    if (r < Rmargin):
        H = H0 * pow(t / t0, -alpha) * pow(1.0 - pow(pow(t / t0, -beta) * (r / R0), (n + 1) / n),
                                           n / (2*n + 1))
        return H
    else:
        return 0.0

class BlatterFlowline(object):

    def __init__(self, Mx, Mz, Lx, Lz, mg_levels, coarsening_factor):
        P = PISM.GridParameters(config)

        P.Mx = Mx
        P.Lx = Lx
        P.x0 = 0.0

        dx = (2 * P.Lx) / (P.Mx - 1)

        P.Ly = dx
        P.My = 3
        P.y0 = 0.0

        P.registration = PISM.CELL_CORNER
        P.periodicity = PISM.Y_PERIODIC
        P.z = PISM.DoubleVector([0.0, Lz])

        P.ownership_ranges_from_options(ctx.size)

        grid = PISM.IceGrid(ctx.ctx, P)

        geometry = PISM.Geometry(grid)

        enthalpy = PISM.IceModelVec3(grid, "enthalpy", PISM.WITHOUT_GHOSTS, grid.z())
        enthalpy.set_attrs("internal", "enthalpy of ice", "J kg-1", "J kg-1", "", 0)

        yield_stress = PISM.IceModelVec2S(grid, "tauc", PISM.WITHOUT_GHOSTS)
        yield_stress.set_attrs("internal", "basal yield stress", "Pa", "Pa", "", 0)

        # this value is not important (we use an isothermal flow law)
        enthalpy.set(1e5)

        self.grid = grid
        self.geometry = geometry
        self.enthalpy = enthalpy
        self.yield_stress = yield_stress
        self.Mz = Mz
        self.mg_levels = mg_levels
        self.coarsening_factor = coarsening_factor

    def initialize(self, ice_thickness, bed_elevation, yield_stress):

        geometry = self.geometry
        tauc = self.yield_stress
        with PISM.vec.Access([geometry.bed_elevation,
                              geometry.ice_thickness,
                              tauc]):
            for (i, j) in self.grid.points():
                x = self.grid.x(i)
                geometry.bed_elevation[i, j] = bed_elevation(x)
                geometry.ice_thickness[i, j] = ice_thickness(x)
                tauc[i, j]                   = yield_stress(x)

        geometry.sea_level_elevation.set(0.0)

        geometry.ensure_consistency(0.0)

    def run(self):

        self.model = PISM.Blatter(self.grid, self.Mz, self.mg_levels, self.coarsening_factor)

        inputs = PISM.StressBalanceInputs()

        inputs.geometry           = self.geometry
        inputs.basal_yield_stress = self.yield_stress
        inputs.enthalpy           = self.enthalpy

        self.model.update(inputs, True)

    def x(self):
        return np.array(self.grid.x())

    def bed(self):
        return self.geometry.bed_elevation.numpy()[1, :]

    def surface(self):
        return self.geometry.ice_surface_elevation.numpy()[1, :]

    def ice_thickness(self):
        return self.geometry.ice_thickness.numpy()[1, :]

    def velocity(self):
        return self.model.velocity_u_sigma().numpy()[1, :, :]
