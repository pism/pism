#!/usr/bin/env python

import numpy as np

import PISM

class Constants(object):
    "Spatially-constant inputs of the LT model"

    tau_c = 1000.0
    "conversion time [s]"

    tau_f = 1000.0
    "fallout time [s]"

    P0 = 0.0
    'Background precipitation rate [mm hr-1]'

    P_scale = 1.0
    'Precipitation scale factor'

    Nm = 0.005
    'moist stability frequency [s-1]'

    Hw = 2500
    'Water vapor scale height [m]'

    latitude = 45.0
    "Latitude used to compute the Coriolis force"

    direction = 270.0
    "Wind direction, 0 is north, 270 is west"

    speed = 15.0
    "Wind speed"

    f = None
    "Coriolis force"

    u = None
    "u component of the wind velocity"

    v = None
    "v component of the wind velocity"

    Cw = None
    "Uplift sensitivity factor [kg m-3]"

    Theta_m = -6.5
    "moist adiabatic lapse rate [K / km]"

    rho_Sref = 7.4e-3
    "reference density [kg m-3]"

    gamma = -5.8
    "adiabatic lapse rate [K / km]"

    def update(self):
        "Update derived constants"

        self.f = 2 * 7.2921e-5 * np.sin(self.latitude * np.pi / 180.0)

        self.u = -np.sin(self.direction * 2 * np.pi / 360) * self.speed
        self.v = np.cos(self.direction * 2 * np.pi / 360) * self.speed

        self.Cw = self.rho_Sref * self.Theta_m / self.gamma

    def __init__(self):
        self.update()

def orographic_precipitation(orography, dx, dy, truncate=True):

    constants = Constants()

    eps = 1e-18

    ny, nx = orography.shape

    pad = 250

    h = np.pad(orography, pad, 'constant')
    ny, nx = h.shape

    h_hat = np.fft.fft2(h)

    x_freq = np.fft.fftfreq(nx, (ny * dx) / (2 * np.pi * nx))
    y_freq = np.fft.fftfreq(ny, (nx * dy) / (2 * np.pi * ny))

    kx, ky = np.meshgrid(x_freq, y_freq)

    # Intrinsic frequency sigma = U*k + V*l
    u0 = constants.u
    v0 = constants.v

    # $\sigma = U k + V l$, see paragraph after eq 5.
    sigma = u0 * kx + v0 * ky

    denominator = sigma**2 - constants.f**2
    denominator[np.logical_and(np.fabs(denominator) < eps, denominator >= 0)] = eps
    denominator[np.logical_and(np.fabs(denominator) < eps, denominator  < 0)] = -eps

    m_squared = (constants.Nm**2 - sigma**2) * (kx**2 + ky**2) / denominator

    m = np.sqrt(np.array(m_squared, dtype=np.complex))

    # Regularization
    nonzero = np.logical_and(m_squared >= 0, sigma != 0)
    m[nonzero] *= np.sign(sigma[nonzero])

    P_hat = h_hat * (constants.Cw * 1j * sigma /
                     ((1 - 1j * m * constants.Hw) *
                      (1 + 1j * sigma * constants.tau_c) *
                      (1 + 1j * sigma * constants.tau_f)))

    # Convert from wave domain back to space domain
    P = np.real(np.fft.ifft2(P_hat))

    # Remove padding
    if pad > 0:
        P = P[pad:-pad, pad:-pad]

    # convert to mm hr-1
    P *= 3600

    # Add background precipitation
    P += constants.P0

    # Truncate
    if truncate:
        P[P < 0] = 0.0

    P *= constants.P_scale

    return P

def create_grid():
    "Allocate the grid for the synthetic geometry test."
    x_min   = -100e3
    x_max   = 200e3
    y_min   = -150e3
    y_max   = 150e3
    dx      = 750
    dy      = 750

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

def gaussian_bump(x, y, h_max=500.0,
                  x0=-25e3, y0=0.0, sigma_x=15e3, sigma_y=15e3):
    "Create the setup needed to reproduce Fig 4c in SB2004"
    X, Y = np.meshgrid(x, y)
    surface = h_max * np.exp(-(((X - x0)**2 / (2 * sigma_x**2)) +
                              ((Y - y0)**2 / (2 * sigma_y**2))))
    return surface

def orographic_precipitation_pism(grid, surface_elevation):
    "Run the PISM implementation of the model to compare to the Python version."

    config = PISM.Context().config

    c = Constants()

    # set parameters to match Constants above
    config.set_boolean("atmosphere.orographic_precipitation.truncate", True)
    config.set_double("atmosphere.orographic_precipitation.conversion_time", c.tau_c)
    config.set_double("atmosphere.orographic_precipitation.coriolis_latitude", c.latitude)
    config.set_double("atmosphere.orographic_precipitation.fallout_time", c.tau_f)
    config.set_double("atmosphere.orographic_precipitation.grid_size_factor", 2)
    config.set_double("atmosphere.orographic_precipitation.lapse_rate", c.gamma)
    config.set_double("atmosphere.orographic_precipitation.moist_adiabatic_lapse_rate", c.Theta_m)
    config.set_double("atmosphere.orographic_precipitation.moist_stability_frequency", c.Nm)
    config.set_double("atmosphere.orographic_precipitation.reference_density", c.rho_Sref)
    config.set_double("atmosphere.orographic_precipitation.water_vapor_scale_height", c.Hw)
    config.set_double("atmosphere.orographic_precipitation.wind_direction", c.direction)
    config.set_double("atmosphere.orographic_precipitation.wind_speed", c.speed)

    config.set_double("atmosphere.orographic_precipitation.background_precip_post", 0.0)
    config.set_double("atmosphere.orographic_precipitation.background_precip_pre", c.P0)
    config.set_double("atmosphere.orographic_precipitation.scale_factor", c.P_scale)

    uniform  = PISM.AtmosphereUniform(grid)
    model    = PISM.AtmosphereOrographicPrecipitation(grid, uniform)
    geometry = PISM.Geometry(grid)

    H          = geometry.ice_thickness.allocate_proc0_copy()
    H.get()[:] = surface_elevation
    geometry.ice_thickness.get_from_proc0(H.get())

    geometry.bed_elevation.set(0.0)
    geometry.sea_level_elevation.set(0.0)
    geometry.ice_area_specific_volume.set(0.0)
    geometry.cell_area.set(grid.dx() * grid.dy())

    # compute surface elevation from ice thickness and bed elevation
    geometry.ensure_consistency(0)

    model.init(geometry)
    model.update(geometry, 0, 1)

    return model.mean_precipitation().numpy() * 3600

def orographic_precipitation_test():
    "Comparing PISM's orographic precipitation model to a reference implementation"

    grid = create_grid()
    orography = gaussian_bump(grid.x(), grid.y())

    P_python = orographic_precipitation(orography, grid.dx(), grid.dy())

    P_pism = orographic_precipitation_pism(grid, orography)

    np.testing.assert_allclose(P_python, P_pism, rtol=0, atol=1e-4)

if __name__ == "__main__":

    grid = create_grid()

    x = grid.x()
    y = grid.y()

    orography = gaussian_bump(x, y)

    P_python = orographic_precipitation(orography, grid.dx(), grid.dy())

    P_pism = orographic_precipitation_pism(grid, orography)

    import pylab as plt

    def plot(v, title):
        plt.figure()
        plt.contour(x, y, v, levels=np.linspace(0.025, 2.025, 6))
        plt.colorbar()
        plt.contour(x, y, orography, colors="black", linewidths=0.5)
        plt.title(title)

    plot(P_pism, "PISM")
    plot(P_python, "Python")

    plt.figure()
    plt.title("difference")
    plt.imshow(P_pism - P_python)
    plt.colorbar()

    plt.show()
