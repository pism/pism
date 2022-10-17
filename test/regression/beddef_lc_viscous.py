#!/usr/bin/env python3

"""Solves equation (16) of BuelerLingleBrown to obtain the steady
state viscous plate deflection corresponding to a given disc load.

Used as a verification (and regression) test for LingleClarkSerial::bootstrap().
"""

import PISM
import numpy as np
from PISM.util import convert

config = PISM.Context().config
log = PISM.Context().log

config.set_number("bed_deformation.lc.grid_size_factor", 2)
config.set_flag("bed_deformation.lc.elastic_model", False)

# constants
standard_gravity = config.get_number("constants.standard_gravity")
ice_density = config.get_number("constants.ice.density")
mantle_density = config.get_number("bed_deformation.mantle_density")
mantle_viscosity = config.get_number("bed_deformation.mantle_viscosity")
lithosphere_flexural_rigidity = config.get_number("bed_deformation.lithosphere_flexural_rigidity")

# disc load parameters
disc_radius = convert(1000, "km", "m")
disc_thickness = 1000.0         # meters
# domain size
Lx = 2 * disc_radius
# time to use for the comparison
T = convert(1e6,   "years", "second")
t_final = convert(20000, "years", "second")
dt = convert(500,   "years", "second")


def deflection(time, radius, disc_thickness, disc_radius):
    """Compute the viscous plate deflection. See formula (17) in 'Fast
    computation of a viscoelastic deformable Earth model for ice-sheet
    simulations' by Bueler, Lingle, and Brown, 2007.
    """
    return PISM.viscDisc(time, disc_thickness, disc_radius, radius,
                         mantle_density, ice_density, standard_gravity,
                         lithosphere_flexural_rigidity, mantle_viscosity)


def exact(dics_radius, disc_thickness, t, L, N):
    "Evaluate the exact solution at N points."
    r = np.linspace(0, L, N)
    z = [deflection(t, rr, disc_thickness, disc_radius) for rr in r]
    return (r, z)


def modeled_time_dependent(dics_radius, disc_thickness, t_end, L, Nx, dt):
    "Use the LingleClark class to compute plate deflection."

    Ny = Nx
    Mx = int(2 * Nx - 1)
    My = int(2 * Ny - 1)

    ctx = PISM.Context().ctx
    grid = PISM.IceGrid.Shallow(ctx, L, L, 0, 0, Mx, My, PISM.CELL_CORNER, PISM.NOT_PERIODIC)

    bed_model = PISM.LingleClark(grid)

    ice_thickness = PISM.Scalar(grid, "thk")

    bed = PISM.Scalar(grid, "topg")

    bed_uplift = PISM.Scalar(grid, "uplift")

    sea_level = PISM.Scalar(grid, "sea_level")

    # start with a flat bed, no ice, and no uplift
    bed.set(0.0)
    bed_uplift.set(0.0)
    ice_thickness.set(0.0)
    sea_level.set(-1000.0)

    bed_model.bootstrap(bed, bed_uplift, ice_thickness, sea_level)

    # add the disc load
    with PISM.vec.Access(nocomm=ice_thickness):
        for (i, j) in grid.points():
            r = PISM.radius(grid, i, j)
            if r <= disc_radius:
                ice_thickness[i, j] = disc_thickness

    t = 0.0
    while t < t_end:
        # make sure we don't step past t_end
        if t + dt > t_end:
            dt = t_end - t

        bed_model.step(ice_thickness, sea_level, dt)

        t += dt
        log.message(2, ".")
    log.message(2, "\n")

    # extract half of the x grid
    r = grid.x()[Nx-1:]

    # extract values along the x direction (along the radius of the disc)
    z = bed_model.bed_elevation().numpy()[Ny-1, Nx-1:]

    return r, z


def modeled_steady_state(dics_radius, disc_thickness, time, L, Nx):
    "Use the LingleClark class to compute plate deflection."

    Ny = Nx
    Mx = int(2 * Nx - 1)
    My = int(2 * Ny - 1)

    ctx = PISM.Context().ctx
    grid = PISM.IceGrid.Shallow(ctx, L, L, 0, 0, Mx, My, PISM.CELL_CORNER, PISM.NOT_PERIODIC)

    bed_model = PISM.LingleClark(grid)

    ice_thickness = PISM.Scalar(grid, "thk")

    bed = PISM.Scalar(grid, "topg")
    bed.set(0.0)

    bed_uplift = PISM.Scalar(grid, "uplift")
    bed_uplift.set(0.0)

    sea_level = PISM.Scalar(grid, "sea_level")
    sea_level.set(-1000.0)

    with PISM.vec.Access(nocomm=ice_thickness):
        for (i, j) in grid.points():
            r = PISM.radius(grid, i, j)
            if r <= disc_radius:
                ice_thickness[i, j] = disc_thickness
            else:
                ice_thickness[i, j] = 0.0

    bed_model.bootstrap(bed, bed_uplift, ice_thickness, sea_level)

    # extract half of the x grid
    r = grid.x()[Nx-1:]

    # extract values along the x direction (along the radius of the disc)
    z = bed_model.total_displacement().numpy()[Ny-1, Nx-1:]

    return r, z


def compare_steady_state(N):
    "Compare eact and modeled results."
    r_exact, z_exact = exact(disc_radius, disc_thickness, T, Lx, N)
    r, z = modeled_steady_state(disc_radius, disc_thickness, T, Lx, N)

    diff_origin = np.fabs(z_exact[0] - z[0])
    diff_average = np.average(np.fabs(z_exact - z))
    diff_max = np.max(np.fabs(z_exact - z))

    return diff_origin, diff_max, diff_average


def compare_time_dependent(N):
    "Compare exact and modeled results."
    r_exact, z_exact = exact(disc_radius, disc_thickness, t_final, Lx, N)

    dx = r_exact[1] - r_exact[0]
    log.message(2, "N = {}, dx = {} km\n".format(N, dx/1000.0))

    r, z = modeled_time_dependent(disc_radius, disc_thickness, t_final, Lx, N, dt)

    diff_origin = np.fabs(z_exact[0] - z[0])
    diff_average = np.average(np.fabs(z_exact - z))
    diff_max = np.max(np.fabs(z_exact - z))

    return diff_origin, diff_max, diff_average, dx


def time_dependent_test():
    "Time dependent bed deformation (disc load)"
    diff = np.array([compare_time_dependent(n)[:3] for n in [34, 67]])

    stored = [[0.04099917, 5.05854,    0.93909436],
              [0.05710513, 4.14329508, 0.71246272]]

    return np.testing.assert_almost_equal(diff, stored)


def steady_state_test():
    "Steady state bed deformation (disc load)"
    Ns = 10 * np.arange(1, 5) + 1
    diff = np.array([compare_steady_state(n) for n in Ns])

    stored = [[ 0.0399697,  15.71882867,  3.80458833],
              [ 0.04592036, 11.43876195,  1.94967725],
              [ 0.04357962,  9.7207298,   1.76262896],
              [ 0.04019595,  7.71929661,  1.38746767]]

    return np.testing.assert_almost_equal(diff, stored)


def verify_steady_state():
    "Set up a grid refinement study and produce convergence plots."

    Ns = 101 + 10 * np.arange(0, 10)

    diff = np.array([compare_steady_state(n) for n in Ns])

    plt.figure()
    plt.title("Steady state")

    d = np.log10(diff)
    log_n = np.log10(1.0 / Ns)
    for j, label in enumerate(["origin", "max", "average"]):
        p = np.polyfit(log_n, d[:, j], 1)
        plt.plot(log_n, d[:, j], marker='o', label="{}, O(dx^{:.2})".format(label, p[0]))
        plt.plot(log_n, np.polyval(p, log_n), ls="--")

    plt.legend()
    plt.grid(True)
    plt.xlabel("log10(1/N)")
    plt.ylabel("log10(error)")
    plt.title("Convergence rates for the steady-state problem")
    plt.show()


def verify_time_dependent():
    "Set up a spatial grid refinement study and produce convergence plots."

    dxs = [15, 30, 60, 125, 250, 500]
    # in km, same as in BuelerLingleBrown, figure 4
    #
    # Note that we compute max and average errors differently here, so
    # convergence rates produced by this scrips are not directly
    # comparable to the paper.

    Ns = [int(Lx / (1000 * dx)) + 1 for dx in dxs]

    diff = np.array([compare_time_dependent(n) for n in Ns])

    plt.figure()
    plt.title("Time-dependent")

    d = np.log10(diff)
    dx = diff[:, 3] / 1000.0    # convert to km
    log_dx = np.log10(dx)
    for j, label in enumerate(["origin", "max", "average"]):
        p = np.polyfit(log_dx, d[:, j], 1)
        plt.plot(log_dx, d[:, j], marker='o', label="{}, O(dx^{:.2})".format(label, p[0]))
        plt.plot(log_dx, np.polyval(p, log_dx), ls="--")
        plt.xticks(log_dx, ["{:.0f}".format(x) for x in dx])

    plt.legend()
    plt.grid(True)
    plt.xlabel("dx, km")
    plt.ylabel("log10(error)")
    plt.title("Convergence rates for the time-dependent problem")
    plt.show()

if __name__ == "__main__":
    import pylab as plt
    log.message(2, "  Creating convergence plots (spatial refinement)...\n")
    log.message(2, "  1. Steady state problem...\n")
    verify_steady_state()
    log.message(2, "  2. Time-dependent problem...\n")
    verify_time_dependent()
