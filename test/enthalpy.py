"""Simple verification/regression tests for the PISM's enthalpy
solver. Work in progress.
"""

import PISM
import numpy as np
import pylab as plt

class EnthalpyTest(object):
    def __init__(self, Mz=101, dt=0.1):
        ctx = PISM.Context()
        self.EC = PISM.EnthalpyConverter(ctx.config)
        self.grid = ctx.newgrid()

        self.Lx = 1e5
        self.Ly = 1e5
        self.Lz = 1000.0
        self.Mx = 3
        self.My = 3
        self.Mz = Mz
        self.dt = self.grid.convert(dt, "years", "seconds")

        PISM.model.initGrid(self.grid, self.Lx, self.Ly, self.Lz,
                            self.Mx, self.My, self.Mz,
                            PISM.NOT_PERIODIC)

        self.z = np.array(self.grid.zlevels_fine)

        self.enthalpy = PISM.model.createEnthalpyVec(self.grid)


        self.strain_heating = PISM.model.createStrainHeatingVec(self.grid)

        self.u, self.v, self.w = PISM.model.create3DVelocityVecs(self.grid)

        self.esys = PISM.enthSystemCtx(ctx.config, self.enthalpy,
                                       self.grid.dx, self.grid.dy, self.dt,
                                       self.grid.dz_fine, self.grid.Mz_fine,
                                       "enth", self.EC)
        # zero ice velocity:
        self.reset_flow()
        # no strain heating:
        self.reset_strain_heating()

    def reset_flow(self):
        self.u.set(0.0)
        self.v.set(0.0)
        self.w.set(0.0)

    def reset_strain_heating(self):
        self.strain_heating.set(0.0)

    def init_column(self):
        self.esys.initThisColumn(1, 1, False, # NOT marginal (but it does not matter)
                                 self.Lz, self.u, self.v, self.w, self.strain_heating)

def dirichlet_test(dt):
    """Test the enthalpy solver with Dirichlet B.C. at the base and
    surface.
    """
    T = EnthalpyTest(dt=dt)

    E_surface = T.EC.getEnth(270.0, 0.0, 0.0)
    E_base = T.EC.getEnth(230.0, 0.0, T.EC.getPressureFromDepth(T.Lz))

    T.enthalpy.set(E_base)

    def E_exact(z):
        """Exact solution is a straight line connecting basal and surface
        conditions."""
        return E_base + (z/T.Lz) * (E_surface - E_base)

    with PISM.vec.Access(nocomm=[T.enthalpy, T.u, T.v, T.w, T.strain_heating]):
        T.init_column()

        T.esys.setDirichletSurface(E_surface)
        T.esys.setDirichletBasal(E_base)

        x = T.esys.solveThisColumn()

        plt.plot(T.z, x, label="dt={}".format(dt))
        plt.plot(T.z, E_exact(T.z))

def neumann_bc_base_test(dt):
    """Test the enthalpy solver with Neumann B.C. at the base and
    Dirichlet B.C. at the surface.
    """

    T = EnthalpyTest(dt=dt)

    E_surface = T.EC.getEnth(270.0, 0.0, 0.0)

    T.enthalpy.set(0.0)

    def E_exact(z):
        """In the cold case, the exact solution is a straight line with the
        slope -G/K (G is the geothermal flux) passing through the surface B.C.
        """
        pass

    with PISM.vec.Access(nocomm=[T.enthalpy, T.u, T.v, T.w, T.strain_heating]):
        T.init_column()

        T.esys.setDirichletSurface(E_surface)
        T.esys.setBasalHeatFlux(0.0)

        x = T.esys.solveThisColumn()

        plt.plot(T.z, x, label="dt={}".format(dt))

if __name__ == "__main__":
    plt.figure(1)
    plt.hold(True)
    for dt in [10**x for x in range(1,6)]:
        dirichlet_test(dt)
    plt.grid(True)
    plt.legend()

    plt.figure(2)
    plt.hold(True)
    for dt in [10**x for x in range(1,6)]:
        neumann_bc_base_test(dt)
    plt.grid(True)
    plt.legend()

    plt.show()
