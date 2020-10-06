#!/usr/bin/env python3
import numpy as np

import PISM

# silence initialization messages
PISM.Context().log.set_threshold(1)

def triangle_ridge_grid(dx=5e4, dy=5e4):
    "Allocate the grid for the synthetic geometry test."
    x_min, x_max   = -100e3, 100e3
    y_min, y_max   = -100e3, 100e3

    x0 = (x_max + x_min) / 2.0
    y0 = (y_max + y_min) / 2.0

    Lx = (x_max - x_min) / 2.0
    Ly = (y_max - y_min) / 2.0
    Mx = int((x_max - x_min) / dx)
    My = int((y_max - y_min) / dy)

    return PISM.IceGrid_Shallow(PISM.Context().ctx,
                                Lx, Ly,
                                x0, y0,
                                Mx, My,
                                PISM.CELL_CORNER, PISM.NOT_PERIODIC)

def triangle_ridge(x, A=500.0, d=50e3):
    "Create the 'triangle ridge' topography"
    return np.maximum(A * (1 - np.fabs(x) / d), 0)

def triangle_ridge_exact(x, u, Cw, tau, A=500.0, d=50e3):
    """Exact precipitation for the triangle ridge setup.

    See equations 44, 45, 46 in Smith and Barstad (2004).
    """
    assert d > 0

    C = Cw * u * A / d
    Ut = u * tau

    xc = Ut * np.log(2 - np.exp(-d / Ut))

    def P(x):
        if x < 0 and x >= -d:
            return C * (1.0 - np.exp(-(x + d) / Ut))
        elif x >= 0 and x <= xc:
            return C * (np.exp(-x / Ut) * (2 - np.exp(-d / Ut)) - 1)
        else:
            return 0

    try:
        return np.array([P(t) for t in x])
    except TypeError:
        return P(x)

def run_model(grid, orography):
    "Run the PISM implementation of the model to compare to the Python version."

    model    = PISM.AtmosphereOrographicPrecipitation(grid, PISM.AtmosphereUniform(grid))
    geometry = PISM.Geometry(grid)

    with PISM.vec.Access(nocomm=geometry.ice_thickness):
        for i,j in grid.points():
            geometry.ice_thickness[i, j] = orography[j, i]

    geometry.bed_elevation.set(0.0)
    geometry.sea_level_elevation.set(0.0)
    geometry.ice_area_specific_volume.set(0.0)

    # compute surface elevation from ice thickness and bed elevation
    geometry.ensure_consistency(0)

    model.init(geometry)
    model.update(geometry, 0, 1)

    config = PISM.Context().config
    water_density = config.get_number("constants.fresh_water.density")

    # convert from kg / (m^2 s) to mm/s
    return model.mean_precipitation().numpy() / (1e-3 * water_density)

def max_error(spacing, wind_direction):
    # Set conversion time to zero (we could set fallout time to zero instead: it does not
    # matter which one is zero)

    wind_speed = 15

    config = PISM.Context().config

    # set wind speed and direction
    config.set_number("atmosphere.orographic_precipitation.wind_speed", wind_speed)
    config.set_number("atmosphere.orographic_precipitation.wind_direction", wind_direction)

    # set conversion time to zero
    config.set_number("atmosphere.orographic_precipitation.conversion_time", 0.0)
    # eliminate the effect of airflow dynamics
    config.set_number("atmosphere.orographic_precipitation.water_vapor_scale_height", 0.0)
    # eliminate the effect of the Coriolis force
    config.set_number("atmosphere.orographic_precipitation.coriolis_latitude", 0.0)

    # get constants needed to compute the exact solution
    tau      = config.get_number("atmosphere.orographic_precipitation.fallout_time")
    Theta_m  = config.get_number("atmosphere.orographic_precipitation.moist_adiabatic_lapse_rate")
    rho_Sref = config.get_number("atmosphere.orographic_precipitation.reference_density")
    gamma    = config.get_number("atmosphere.orographic_precipitation.lapse_rate")
    Cw       = rho_Sref * Theta_m / gamma

    if wind_direction == 90 or wind_direction == 270:
        # east or west
        grid = triangle_ridge_grid(dx=spacing)
        t = np.array(grid.x())

        h = triangle_ridge(t)
        orography = np.tile(h, (grid.My(), 1))

        P = run_model(grid, orography)
        P = P[grid.My() // 2, :]
    else:
        # north or south
        grid = triangle_ridge_grid(dy=spacing)
        t = np.array(grid.y())

        h = triangle_ridge(t)
        orography = np.tile(h, (grid.Mx(), 1)).T

        P = run_model(grid, orography)
        P = P[:, grid.Mx() // 2]

    if wind_direction == 0 or wind_direction == 90:
        P_exact = triangle_ridge_exact(-t, wind_speed, Cw, tau)
    else:
        P_exact = triangle_ridge_exact(t, wind_speed, Cw, tau)

    return np.max(np.fabs(P - P_exact))

def convergence_rate(dxs, error, wind_direction, plot):
    errors = [error(dx, wind_direction) for dx in dxs]

    p = np.polyfit(np.log10(dxs), np.log10(errors), 1)

    if plot:
        import pylab as plt

        direction = {0 : "north", 90 : "east", 180 : "south", 270 : "west"}

        def log_plot(x, y, style, label):
            plt.plot(np.log10(x), np.log10(y), style, label=label)
            plt.xticks(np.log10(x), x)

        def log_fit_plot(x, p, label):
            plt.plot(np.log10(x), np.polyval(p, np.log10(x)), label=label)

        plt.figure()
        plt.title("Precipitation errors (wind direction: {})".format(direction[wind_direction]))
        log_fit_plot(dxs, p, "polynomial fit (dx^{:1.4})".format(p[0]))
        log_plot(dxs, errors, 'o', "errors")
        plt.legend()
        plt.grid()
        plt.xlabel("grid spacing (meters)")
        plt.ylabel("log10(error)")
        plt.show()

    return p[0]

def ltop_test(dxs=[2000, 1000, 500], plot=False):
    "Orographic precipitation (triangle ridge test case)"

    assert convergence_rate(dxs, max_error,   0, plot) > 1.99
    assert convergence_rate(dxs, max_error,  90, plot) > 1.99
    assert convergence_rate(dxs, max_error, 180, plot) > 1.99
    assert convergence_rate(dxs, max_error, 270, plot) > 1.99

if __name__ == "__main__":
    ltop_test(dxs=[2000, 1000, 500, 250, 125], plot=True)
