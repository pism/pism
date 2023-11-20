#!/usr/bin/env python3

"""Simple verification/regression tests for the heat equation in bedrock columns. Tests
PISM's "bedrock thermal unit".

"""

import numpy as np
import PISM
from PISM.util import convert

log10 = np.log10

ctx = PISM.Context()

k = ctx.config.get_number("energy.bedrock_thermal.conductivity")
c = ctx.config.get_number("energy.bedrock_thermal.specific_heat_capacity")
rho = ctx.config.get_number("energy.bedrock_thermal.density")
K = k / c
# alpha squared
alpha2 = k / (c * rho)

# set to True to plot solutions (for debugging)
plot_solutions = False


def convergence_rate_time(error_func, plot):
    "Compute the convergence rate with refinement in time."
    dts = 2.0 ** np.arange(4, 10)

    max_errors = np.zeros_like(dts)
    avg_errors = np.zeros_like(dts)
    for k, dt in enumerate(dts):
        max_errors[k], avg_errors[k] = error_func(plot_solutions, dt_years=dt, Mz=101)

    p_max = np.polyfit(log10(dts), log10(max_errors), 1)
    p_avg = np.polyfit(log10(dts), log10(avg_errors), 1)

    if plot:
        plt.figure()
        plt.title("Heat equation in the bedrock column\nconvergence as dt -> 0")
        log_plot(dts, max_errors, "o", "max errors")
        log_plot(dts, avg_errors, "o", "avg errors")
        log_fit_plot(dts, p_max, "max: dt^{:.3}".format(p_max[0]))
        log_fit_plot(dts, p_avg, "avg: dt^{:.3}".format(p_avg[0]))
        plt.axis("tight")
        plt.grid(True)
        plt.legend(loc="best")

    return p_max[0], p_avg[0]


def convergence_rate_space(error_func, plot):
    "Compute the convergence rate with refinement in space."
    Mz = np.array(2.0 ** np.arange(2, 7), dtype=int)
    dzs = 1000.0 / Mz

    max_errors = np.zeros_like(dzs)
    avg_errors = np.zeros_like(dzs)
    for k, M in enumerate(Mz):
        T = 1000.0
        # time step has to be short enough so that errors due to the time discretization
        # are smaller than errors due to the spatial discretization
        dt = 0.001 * T
        max_errors[k], avg_errors[k] = error_func(
            plot_solutions, T_final_years=T, dt_years=dt, Mz=M
        )

    p_max = np.polyfit(log10(dzs), log10(max_errors), 1)
    p_avg = np.polyfit(log10(dzs), log10(avg_errors), 1)

    if plot:
        plt.figure()
        plt.title("Heat equation in the bedrock column\nconvergence as dz -> 0")
        log_plot(dzs, max_errors, "o", "max errors")
        log_plot(dzs, avg_errors, "o", "avg errors")
        log_fit_plot(dzs, p_max, "max: dz^{:.3}".format(p_max[0]))
        log_fit_plot(dzs, p_avg, "avg: dz^{:.3}".format(p_avg[0]))
        plt.axis("tight")
        plt.grid(True)
        plt.legend(loc="best")

    return p_max[0], p_avg[0]


def exact(L, Q_bottom, U_top):
    """Exact solution (and an initial state) for the 'Neumann at the base,
    Dirichlet at the top' setup."""
    n = 2
    lambda_n = 1.0 / L * (-np.pi / 2.0 + n * np.pi)
    a = L * 25.0

    def f(z, t):
        v = a * np.exp(-(lambda_n**2) * alpha2 * t) * np.sin(lambda_n * z)
        return v + (U_top + Q_bottom * z)

    return f


def errors(plot_results=True, T_final_years=1000.0, dt_years=100, Mz=101):
    """Test the bedrock temperature solver with Neumann B.C. at the base and
    Dirichlet B.C. at the top surface.
    """
    T_final = convert(T_final_years, "years", "seconds")
    dt = convert(dt_years, "years", "seconds")

    Lz = 1000.0
    dz = Lz / (Mz - 1.0)

    column = PISM.BedrockColumn("btu", ctx.config, dz, int(Mz))

    z = np.linspace(-Lz, 0, Mz)

    T_base = 240.0  # Kelvin
    T_surface = 260.0  # Kelvin
    dT_base = (T_surface - T_base) / Lz

    T_steady = T_base + dT_base * (z - (-Lz))
    Q_base = -K * dT_base

    T_exact = exact(Lz, dT_base, T_surface)

    t = 0.0
    # initial condition
    x = T_exact(z, t)
    while t < T_final:
        x = column.solve(dt, Q_base, T_surface, x)
        t += dt

    T_exact_final = T_exact(z, t)

    if plot_results:
        t_years = convert(t, "seconds", "years")

        plt.figure()
        plt.xlabel("z, meters")
        plt.ylabel("T, K")
        plt.step(z, T_exact(z, 0), color="blue", label="initial condition")
        plt.step(z, T_exact_final, color="green", label="exact solution")
        plt.step(z, T_steady, "--", color="black", label="steady state profile")
        plt.grid(True)

        plt.step(z, x, label="T={} years".format(t_years), color="red")

        plt.legend(loc="best")

    errors = T_exact(z, t) - x

    max_error = np.max(np.fabs(errors))
    avg_error = np.average(np.fabs(errors))

    return max_error, avg_error


def test(plot=False):
    assert convergence_rate_time(errors, plot)[1] > 0.94
    assert convergence_rate_space(errors, plot)[1] > 1.89


if __name__ == "__main__":
    import pylab as plt

    def log_plot(x, y, style, label):
        plt.plot(log10(x), log10(y), style, label=label)
        plt.xticks(log10(x), x)

    def log_fit_plot(x, p, label):
        plt.plot(log10(x), np.polyval(p, log10(x)), label=label)

    test(plot=True)
    plt.show()
