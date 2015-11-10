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
        sys = ctx.unit_system
        config = ctx.config

        config.set_string("grid_ice_vertical_spacing_type", "equal")

        self.Lz = 1000.0
        self.z = np.linspace(0, self.Lz, Mz)

        param = PISM.GridParameters()
        param.Lx = 1e5
        param.Ly = 1e5
        param.z = PISM.DoubleVector(self.z)
        param.Mx = 3
        param.My = 3
        param.Mz = Mz

        param.ownership_ranges_from_options(1)

        self.dt = PISM.convert(sys, dt, "years", "seconds")

        self.grid = PISM.IceGrid(ctx.ctx, param)
        grid = self.grid

        self.enthalpy = PISM.model.createEnthalpyVec(grid)

        self.strain_heating = PISM.model.createStrainHeatingVec(grid)

        self.u, self.v, self.w = PISM.model.create3DVelocityVecs(grid)

        self.esys = PISM.enthSystemCtx(grid.z(), "enth",
                                       grid.dx(), grid.dy(), self.dt,
                                       config,
                                       self.enthalpy,
                                       self.u, self.v, self.w,
                                       self.strain_heating,
                                       self.EC)
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
        self.esys.init(1, 1, False, self.Lz)  # NOT marginal (but it does not matter)


def dirichlet_test(dt):
    """Test the enthalpy solver with Dirichlet B.C. at the base and
    surface.
    """
    T = EnthalpyTest(dt=dt)

    E_surface = T.EC.enthalpy(270.0, 0.0, 0.0)
    E_base = T.EC.enthalpy(230.0, 0.0, T.EC.pressure(T.Lz))

    T.enthalpy.set(E_base)

    def E_exact(z):
        """Exact solution is a straight line connecting basal and surface
        conditions."""
        return E_base + (z / T.Lz) * (E_surface - E_base)

    with PISM.vec.Access(nocomm=[T.enthalpy, T.u, T.v, T.w, T.strain_heating]):
        T.init_column()

        T.esys.set_surface_dirichlet(E_surface)
        T.esys.set_basal_dirichlet(E_base)

        x = T.esys.solve()

        plt.step(T.z, x, label="dt={}".format(dt))

    return T.z, E_exact(T.z)


def neumann_bc_base_test(dt):
    """Test the enthalpy solver with Neumann B.C. at the base and
    Dirichlet B.C. at the surface.
    """

    T = EnthalpyTest(dt=dt)

    E_surface = T.EC.enthalpy(270.0, 0.0, 0.0)

    T.enthalpy.set(0.0)

    def E_exact(z):
        """In the cold case, the exact solution is a straight line with the
        slope -G/K (G is the geothermal flux) passing through the surface B.C.
        """
        return np.zeros_like(z) + E_surface

    with PISM.vec.Access(nocomm=[T.enthalpy, T.u, T.v, T.w, T.strain_heating]):
        T.init_column()

        T.esys.set_surface_dirichlet(E_surface)
        T.esys.set_basal_heat_flux(0.0)

        x = T.esys.solve()

        plt.step(T.z, x, label="dt={}".format(dt))

    return T.z, E_exact(T.z)

def neumann_bc_surface_test(dt):
    """Test the enthalpy solver with Neumann B.C. at the surface and
    Dirichlet B.C. at the base.
    """

    T = EnthalpyTest(dt=dt)

    E_base = T.EC.enthalpy(250.0, 0.0, T.EC.pressure(T.Lz))

    T.enthalpy.set(0.0)

    def E_exact(z):
        """In the cold case, the exact solution is a straight line with the
        slope -G/K (G is the geothermal flux) passing through the base B.C.
        """
        return np.zeros_like(z) + E_base

    with PISM.vec.Access(nocomm=[T.enthalpy, T.u, T.v, T.w, T.strain_heating]):
        T.init_column()

        T.esys.set_basal_dirichlet(E_base)
        T.esys.set_surface_heat_flux(0.0)

        x = T.esys.solve()

        plt.step(T.z, x, label="dt={}".format(dt))

    return T.z, E_exact(T.z)


if __name__ == "__main__":
    plt.figure(1)
    plt.title("Testing Dirichlet B.C. (base and surface)")
    plt.hold(True)
    z, exact = None, None
    for dt in [10 ** x for x in range(1, 6)]:
        z, exact = dirichlet_test(dt)
    plt.step(z, exact, label="exact (steady state)")
    plt.grid(True)
    plt.legend(loc="best")

    plt.figure(2)
    plt.title("Testing Neumann B.C. (base)")
    plt.hold(True)
    for dt in [10 ** x for x in range(1, 6)]:
        z, exact = neumann_bc_base_test(dt)
    plt.step(z, exact, label="exact (steady state)")
    plt.grid(True)
    plt.legend(loc="best")

    plt.figure(3)
    plt.title("Testing Neumann B.C. (surface)")
    plt.hold(True)
    for dt in [10 ** x for x in range(1, 6)]:
        z, exact = neumann_bc_surface_test(dt)
    plt.step(z, exact, label="exact (steady state)")
    plt.grid(True)
    plt.legend(loc="best")

    plt.show()
