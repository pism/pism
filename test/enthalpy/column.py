"""Simple verification/regression tests for vertical diffusion and
advection of enthalpy in the column. Tests PISM's enthalpy solver.

"""

import PISM
from PISM.util import convert
import numpy as np
import pylab as plt


def log_plot(x, y, style, label):
    plt.plot(log10(x), log10(y), style, label=label)
    plt.xticks(log10(x), x)


def log_fit_plot(x, p, label):
    plt.plot(log10(x), np.polyval(p, log10(x)), label=label)


log10 = np.log10

ctx = PISM.Context()
unit_system = ctx.unit_system
config = ctx.config

config.set_string("grid.ice_vertical_spacing", "equal")

k = config.get_double("constants.ice.thermal_conductivity")
c = config.get_double("constants.ice.specific_heat_capacity")
rho = config.get_double("constants.ice.density")
K = k / c
# alpha squared
alpha2 = k / (c * rho)

EC = PISM.EnthalpyConverter(ctx.config)
pressure = np.vectorize(EC.pressure)
cts = np.vectorize(EC.enthalpy_cts)


class EnthalpyColumn(object):
    "Set up the grid and arrays needed to run column solvers"

    def __init__(self, Mz, dt):
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

        self.dt = dt

        self.grid = PISM.IceGrid(ctx.ctx, param)
        grid = self.grid

        self.enthalpy = PISM.model.createEnthalpyVec(grid)

        self.strain_heating = PISM.model.createStrainHeatingVec(grid)

        self.u, self.v, self.w = PISM.model.create3DVelocityVecs(grid)

        self.sys = PISM.enthSystemCtx(grid.z(), "energy.enthalpy",
                                      grid.dx(), grid.dy(), self.dt,
                                      config,
                                      self.enthalpy,
                                      self.u, self.v, self.w,
                                      self.strain_heating,
                                      EC)

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
        ice_thickness = self.Lz
        self.sys.init(1, 1, False, ice_thickness)


def diffusion_convergence_rate_time(title, error_func):
    "Compute the convergence rate with refinement in time."
    dts = 2.0**np.arange(10)

    max_errors = np.zeros_like(dts)
    avg_errors = np.zeros_like(dts)
    for k, dt in enumerate(dts):
        max_errors[k], avg_errors[k] = error_func(False, dt_years=dt, Mz=101)

    p_max = np.polyfit(log10(dts), log10(max_errors), 1)
    p_avg = np.polyfit(log10(dts), log10(avg_errors), 1)

    if True:
        plt.figure()
        plt.title(title + "\nTesting convergence as dt -> 0")
        log_plot(dts, max_errors, 'o', "max errors")
        log_plot(dts, avg_errors, 'o', "avg errors")
        log_fit_plot(dts, p_max, "max: dt^{}".format(p_max[0]))
        log_fit_plot(dts, p_avg, "avg: dt^{}".format(p_avg[0]))
        plt.axis('tight')
        plt.grid(True)
        plt.legend(loc="best")

    return p_max[0], p_avg[0]


def diffusion_convergence_rate_space(title, error_func):
    "Compute the convergence rate with refinement in space."
    Mz = np.array(2.0**np.arange(3, 10), dtype=int)
    dzs = 1000.0 / Mz

    max_errors = np.zeros_like(dzs)
    avg_errors = np.zeros_like(dzs)
    for k, M in enumerate(Mz):
        T = 1.0
        max_errors[k], avg_errors[k] = error_func(False,
                                                  T_final_years=T,
                                                  dt_years=T,
                                                  Mz=M)

    p_max = np.polyfit(log10(dzs), log10(max_errors), 1)
    p_avg = np.polyfit(log10(dzs), log10(avg_errors), 1)

    if True:
        plt.figure()
        plt.title(title + "\nTesting convergence as dz -> 0")
        log_plot(dzs, max_errors, 'o', "max errors")
        log_plot(dzs, avg_errors, 'o', "avg errors")
        log_fit_plot(dzs, p_max, "max: dz^{}".format(p_max[0]))
        log_fit_plot(dzs, p_avg, "avg: dz^{}".format(p_avg[0]))
        plt.axis('tight')
        plt.grid(True)
        plt.legend(loc="best")

    return p_max[0], p_avg[0]


def exact_DN(L, U0, QL):
    """Exact solution (and an initial state) for the 'Dirichlet at the base,
    Neumann at the top' setup."""
    n = 1
    lambda_n = 1.0 / L * (np.pi / 2.0 + n * np.pi)
    a = L * 25.0

    def f(z, t):
        v = a * np.exp(-lambda_n**2 * alpha2 * t) * np.sin(lambda_n * z)
        return v + (U0 + QL * z)
    return f


def errors_DN(plot_results=True, T_final_years=1000.0, dt_years=100, Mz=101):
    """Test the enthalpy solver with Dirichlet B.C. at the base and
    Neumann at the top surface.
    """
    T_final = convert(T_final_years, "years", "seconds")
    dt = convert(dt_years, "years", "seconds")

    column = EnthalpyColumn(Mz, dt)

    Lz = column.Lz
    z = np.array(column.sys.z())

    E_base = EC.enthalpy(230.0, 0.0, EC.pressure(Lz))
    E_surface = EC.enthalpy(270.0, 0.0, 0.0) - 25*Lz
    dE_surface = (E_surface - E_base) / Lz

    E_steady = E_base + dE_surface * z
    Q_surface = K * dE_surface

    E_exact = exact_DN(Lz, E_base, dE_surface)

    with PISM.vec.Access(nocomm=[column.enthalpy,
                                 column.u, column.v, column.w,
                                 column.strain_heating]):
        column.sys.fine_to_coarse(E_exact(z, 0), 1, 1, column.enthalpy)
        column.reset_flow()
        column.reset_strain_heating()

        t = 0.0
        while t < T_final:
            column.init_column()

            column.sys.set_surface_heat_flux(Q_surface)
            column.sys.set_basal_dirichlet_bc(E_base)

            x = column.sys.solve()

            column.sys.fine_to_coarse(x, 1, 1, column.enthalpy)

            t += dt

    E_exact_final = E_exact(z, t)

    if plot_results:
        t_years = convert(t, "seconds", "years")

        plt.figure()
        plt.xlabel("z, meters")
        plt.ylabel("E, J/kg")
        plt.step(z, E_exact(z, 0), color="blue", label="initial condition")
        plt.step(z, E_exact_final, color="green", label="exact solution")
        plt.step(z, cts(pressure(Lz - z)), "--", color="black", label="CTS")
        plt.step(z, E_steady, "--", color="green", label="steady state profile")
        plt.grid(True)

        plt.step(z, x, label="T={} years".format(t_years), color="red")

        plt.legend(loc="best")

    errors = E_exact(z, t) - x

    max_error = np.max(np.fabs(errors))
    avg_error = np.average(np.fabs(errors))

    return max_error, avg_error


def exact_ND(L, Q0, UL):
    """Exact solution (and an initial state) for the 'Dirichlet at the base,
    Neumann at the top' setup."""
    n = 2
    lambda_n = 1.0 / L * (-np.pi / 2.0 + n * np.pi)
    a = L * 25.0

    def f(z, t):
        v = a * np.exp(-lambda_n**2 * alpha2 * t) * np.sin(lambda_n * (L - z))
        return v + (UL + Q0 * (z - L))
    return f


def errors_ND(plot_results=True, T_final_years=1000.0, dt_years=100, Mz=101):
    """Test the enthalpy solver with Neumann B.C. at the base and
    Dirichlet B.C. at the top surface.
    """
    T_final = convert(T_final_years, "years", "seconds")
    dt = convert(dt_years, "years", "seconds")

    column = EnthalpyColumn(Mz, dt)

    Lz = column.Lz
    z = np.array(column.sys.z())

    E_base = EC.enthalpy(240.0, 0.0, EC.pressure(Lz))
    E_surface = EC.enthalpy(260.0, 0.0, 0.0)
    dE_base = (E_surface - E_base) / Lz

    E_steady = E_surface + dE_base * (z - Lz)
    Q_base = - K * dE_base

    E_exact = exact_ND(Lz, dE_base, E_surface)

    with PISM.vec.Access(nocomm=[column.enthalpy,
                                 column.u, column.v, column.w,
                                 column.strain_heating]):
        column.sys.fine_to_coarse(E_exact(z, 0), 1, 1, column.enthalpy)
        column.reset_flow()
        column.reset_strain_heating()

        t = 0.0
        while t < T_final:
            column.init_column()

            column.sys.set_basal_heat_flux(Q_base)
            column.sys.set_surface_dirichlet_bc(E_surface)

            x = column.sys.solve()

            column.sys.fine_to_coarse(x, 1, 1, column.enthalpy)

            t += dt

    E_exact_final = E_exact(z, t)

    if plot_results:
        t_years = convert(t, "seconds", "years")

        plt.figure()
        plt.xlabel("z, meters")
        plt.ylabel("E, J/kg")
        plt.step(z, E_exact(z, 0), color="blue", label="initial condition")
        plt.step(z, E_exact_final, color="green", label="exact solution")
        plt.step(z, cts(pressure(Lz - z)), "--", color="black", label="CTS")
        plt.step(z, E_steady, "--", color="green", label="steady state profile")
        plt.grid(True)

        plt.step(z, x, label="T={} years".format(t_years), color="red")

        plt.legend(loc="best")

    errors = E_exact(z, t) - x

    max_error = np.max(np.fabs(errors))
    avg_error = np.average(np.fabs(errors))

    return max_error, avg_error


def exact_advection(L, w):
    "Exact solution of the 'pure advection' problem."
    C = np.pi / L

    def f(z, t):
        return np.sin(C * (z - w * t))

    def df(z, t):
        return C * np.cos(C * (z - w * t))
    return f, df


def errors_advection_up(plot_results=True, T_final=1000.0, dt=100, Mz=101):
    """Test the enthalpy solver using a 'pure advection' problem with
    Neumann (in-flow) B.C. at the base.

    We use Dirichlet B.C. at the surface but they are irrelevant due
    to upwinding.
    """
    w = 1.0

    config.set_double("constants.ice.thermal_conductivity", 0.0)
    column = EnthalpyColumn(Mz, dt)
    config.set_double("constants.ice.thermal_conductivity", k)

    Lz = column.Lz
    z = np.array(column.sys.z())

    E_exact, dE_exact = exact_advection(Lz, w)

    with PISM.vec.Access(nocomm=[column.enthalpy,
                                 column.u, column.v, column.w,
                                 column.strain_heating]):
        column.sys.fine_to_coarse(E_exact(z, 0), 1, 1, column.enthalpy)
        column.reset_flow()
        column.w.set(w)
        column.reset_strain_heating()

        t = 0.0
        while t < T_final:
            column.init_column()

            column.sys.set_basal_neumann_bc(dE_exact(0, t+dt))
            column.sys.set_surface_dirichlet_bc(E_exact(Lz, t+dt))

            x = column.sys.solve()

            column.sys.fine_to_coarse(x, 1, 1, column.enthalpy)

            t += dt

    if plot_results:
        plt.figure()
        plt.xlabel("z, meters")
        plt.ylabel("E, J/kg")
        plt.step(z, E_exact(z, 0), color="blue", label="initial condition")
        plt.step(z, E_exact(z, t), color="green", label="exact solution")
        plt.grid(True)

        plt.step(z, x, label="T={} seconds".format(t), color="red")

        plt.legend(loc="best")

    errors = E_exact(z, t) - x

    max_error = np.max(np.fabs(errors))
    avg_error = np.average(np.fabs(errors))

    return max_error, avg_error


def errors_advection_down(plot_results=True, T_final=1000.0, dt=100, Mz=101):
    """Test the enthalpy solver using a 'pure advection' problem with
    Neumann (in-flow) B.C. at the surface.

    We use Dirichlet B.C. at the base but they are irrelevant due
    to upwinding.
    """
    w = -1.0

    config.set_double("constants.ice.thermal_conductivity", 0.0)
    column = EnthalpyColumn(Mz, dt)
    config.set_double("constants.ice.thermal_conductivity", k)

    Lz = column.Lz
    z = np.array(column.sys.z())

    E_exact, dE_exact = exact_advection(Lz, w)

    with PISM.vec.Access(nocomm=[column.enthalpy,
                                 column.u, column.v, column.w,
                                 column.strain_heating]):
        column.sys.fine_to_coarse(E_exact(z, 0), 1, 1, column.enthalpy)
        column.reset_flow()
        column.w.set(w)
        column.reset_strain_heating()

        t = 0.0
        while t < T_final:
            column.init_column()

            column.sys.set_basal_dirichlet_bc(E_exact(0, t+dt))
            column.sys.set_surface_neumann_bc(dE_exact(Lz, t+dt))

            x = column.sys.solve()

            column.sys.fine_to_coarse(x, 1, 1, column.enthalpy)

            t += dt

    if plot_results:
        plt.figure()
        plt.xlabel("z, meters")
        plt.ylabel("E, J/kg")
        plt.step(z, E_exact(z, 0), color="blue", label="initial condition")
        plt.step(z, E_exact(z, t), color="green", label="exact solution")
        plt.grid(True)

        plt.step(z, x, label="T={} seconds".format(t), color="red")

        plt.legend(loc="best")

    errors = E_exact(z, t) - x

    max_error = np.max(np.fabs(errors))
    avg_error = np.average(np.fabs(errors))

    return max_error, avg_error


def advection_convergence_rate_time(title, error_func):
    "Compute the convergence rate with refinement in time."
    dts = np.linspace(1, 101, 11)
    Mz = 1001
    T_final = 200

    max_errors = np.zeros_like(dts)
    avg_errors = np.zeros_like(dts)
    for k, dt in enumerate(dts):
        max_errors[k], avg_errors[k] = error_func(False,
                                                  T_final=T_final,
                                                  dt=dt,
                                                  Mz=Mz)

    p_max = np.polyfit(log10(dts), log10(max_errors), 1)
    p_avg = np.polyfit(log10(dts), log10(avg_errors), 1)

    if True:
        plt.figure()
        plt.title(title + "\nTesting convergence as dt -> 0")
        log_plot(dts, max_errors, 'o', "max errors")
        log_plot(dts, avg_errors, 'o', "avg errors")
        log_fit_plot(dts, p_max, "max: dt^{}".format(p_max[0]))
        log_fit_plot(dts, p_avg, "avg: dt^{}".format(p_avg[0]))
        plt.axis('tight')
        plt.grid(True)
        plt.legend(loc="best")

    return p_max[0], p_avg[0]


def advection_convergence_rate_space(title, error_func):
    "Compute the convergence rate with refinement in time."
    dt = 0.4
    Mzs = np.linspace(500, 5000, 11, dtype="i")
    T_final = 10

    dzs = 1000.0 / (Mzs - 1)

    max_errors = np.zeros_like(dzs)
    avg_errors = np.zeros_like(dzs)
    for k, M in enumerate(Mzs):
        max_errors[k], avg_errors[k] = error_func(False,
                                                  T_final=T_final,
                                                  dt=dt,
                                                  Mz=M)

    p_max = np.polyfit(log10(dzs), log10(max_errors), 1)
    p_avg = np.polyfit(log10(dzs), log10(avg_errors), 1)

    if True:
        plt.figure()
        plt.title(title + "\nTesting convergence as dz -> 0")
        log_plot(dzs, max_errors, 'o', "max errors")
        log_plot(dzs, avg_errors, 'o', "avg errors")
        log_fit_plot(dzs, p_max, "max: dz^{}".format(p_max[0]))
        log_fit_plot(dzs, p_avg, "avg: dz^{}".format(p_avg[0]))
        plt.axis('tight')
        plt.grid(True)
        plt.legend(loc="best")

    return p_max[0], p_avg[0]


def diffusion_DN_test():
    assert diffusion_convergence_rate_time("Diffusion: Dirichlet at the base, Neumann at the surface",
                                           errors_DN)[1] > 0.93
    assert diffusion_convergence_rate_space("Diffusion: Dirichlet at the base, Neumann at the surface",
                                            errors_DN)[1] > 2.0


def diffusion_ND_test():
    assert diffusion_convergence_rate_time("Diffusion: Neumann at the base, Dirichlet at the surface",
                                           errors_ND)[1] > 0.93
    assert diffusion_convergence_rate_space("Diffusion: Neumann at the base, Dirichlet at the surface",
                                            errors_ND)[1] > 2.0


def advection_up_test():
    assert advection_convergence_rate_time("Advection: Upward flow",
                                           errors_advection_up)[1] > 0.87
    assert advection_convergence_rate_space("Advection: Upward flow",
                                            errors_advection_up)[1] > 0.96


def advection_down_test():
    assert advection_convergence_rate_time("Advection: Downward flow",
                                           errors_advection_down)[1] > 0.87
    assert advection_convergence_rate_space("Advection: Downward flow",
                                            errors_advection_down)[1] > 0.96


if __name__ == "__main__":
    diffusion_ND_test()
    diffusion_DN_test()
    advection_up_test()
    advection_down_test()
    plt.show()
