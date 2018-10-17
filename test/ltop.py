#!/usr/bin/env python

import PISM

import pylab as plt
import numpy as np

def gaussian_bump(x, y, h_max=500.0,
                  x0=-25e3, y0=0.0, sigma_x=15e3, sigma_y=15e3):
    "Create the setup needed to reproduce Fig 4c in SB2004"
    X, Y = np.meshgrid(x, y)
    surface = h_max * np.exp(-(((X - x0)**2 / (2 * sigma_x**2)) +
                              ((Y - y0)**2 / (2 * sigma_y**2))))
    return surface

def grid():
    "Allocate the grid for the synthetic geometry test."
    x_min   = -100e3
    x_max   = 200e3
    y_min   = -150e3
    y_max   = 150e3
    dx      = 750
    dy      = 750

    Lx = (x_max - x_min) / 2.0
    Ly = (y_max - y_min) / 2.0
    Mx = int((x_max - x_min) / dx)
    My = int((y_max - y_min) / dy)

    ctx = PISM.Context()

    return PISM.IceGrid_Shallow(ctx.ctx, Lx, Ly, 0, 0, Mx, My, PISM.CELL_CENTER, PISM.NOT_PERIODIC)

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

    def update(self):
        "Update derived constants"

        Theta_m  = -6.5     # K / km
        rho_Sref = 7.4e-3   # kg m-3
        gamma    = -5.8     # K / km

        self.f = 2 * 7.2921e-5 * np.sin(self.latitude * np.pi / 180.0)

        self.u = -np.sin(self.direction * 2 * np.pi / 360) * self.speed
        self.v = np.cos(self.direction * 2 * np.pi / 360) * self.speed

        self.Cw = rho_Sref * Theta_m / gamma

    def __init__(self):
        self.update()

def orographic_precipitation(dx, dy, orography, constants, truncate):
    """Calculates orographic precipitation following Smith & Barstad
    (2004) (Python implementation).

    """
    eps = 1e-18
    pad = 250
    spy = 31556925.9747

    ny, nx = orography.shape

    h = np.pad(orography, pad, 'constant')
    nx, ny = h.shape

    h_hat = np.fft.fft2(h)

    x_n_value = np.fft.fftfreq(ny, (1.0 / ny))
    y_n_value = np.fft.fftfreq(nx, (1.0 / nx))

    x_len = nx * dx
    y_len = ny * dy
    kx_line = 2 * np.pi * x_n_value / float(x_len)
    ky_line = 2 * np.pi * y_n_value / float(y_len)

    kx, ky = np.meshgrid(kx_line, ky_line)

    # Intrinsic frequency sigma = U*k + V*l
    u0 = constants.u
    v0 = constants.v

    # $\sigma = U k + V l$, see paragraph after eq 5.
    sigma = u0 * kx + v0 * ky

    sigma_c = sigma**2 - constants.f**2

    # The vertical wave number
    # Eqn. 12
    # Regularization
    sigma_c[np.logical_and(np.fabs(sigma_c) < eps, sigma_c >= 0)] = eps
    sigma_c[np.logical_and(np.fabs(sigma_c) < eps, sigma_c  < 0)] = -eps

    sigma_c[sigma_c < 0] = eps
    
    # $(k_x^2 + k_y^2) (N_m^2 - \sigma^2) / \sigma^2$
    m_squared = (constants.Nm**2 - sigma**2) * (kx**2 + ky**2) / sigma_c
    m_squared[m_squared < 0] = 0
    # where does this -1 come from?
    m = np.sqrt(-m_squared)

    # Regularization
    m[np.logical_and(m_squared >= 0, sigma == 0)] = np.sqrt(m_squared[np.logical_and(m_squared >= 0, sigma == 0)])
    
    m[np.logical_and(m_squared >= 0, sigma != 0)] = np.sqrt(m_squared[np.logical_and(
       m_squared >= 0, sigma != 0)]) * np.sign(sigma[np.logical_and(m_squared >= 0, sigma != 0)])

    P_hat = (constants.Cw * 1j * sigma * h_hat /
             ((1 - 1j * m * constants.Hw) *
              (1 + 1j * sigma * constants.tau_c) *
              (1 + 1j * sigma * constants.tau_f)))

    # Converting from wave domain back to space domain
    P = np.fft.ifft2(P_hat)

    # Remove padding
    P = P[pad:-pad, pad:-pad]
    P = np.multiply(np.real(P), 3600)   # convert to mm hr-1

    # Add background precipitation
    P += constants.P0

    # Truncate
    if truncate:
        P[P < 0] = 0.0

    P *= constants.P_scale

    return P

def ltop_python(grid, truncate=True):
    "Run the Python implementation of the model on the synthetic geometry."
    return orographic_precipitation(grid.dx(),
                                    grid.dy(),
                                    gaussian_bump(grid.x(), grid.y()),
                                    Constants(),
                                    truncate)

def ltop_pism(grid):
    "Run the PISM implementation of the model to compare to the Python version."
    uniform  = PISM.AtmosphereUniform(grid)
    model    = PISM.AtmosphereOrographicPrecipitation(grid, uniform)
    geometry = PISM.Geometry(grid)

    H = geometry.ice_thickness.allocate_proc0_copy()
    H.get()[:] = gaussian_bump(grid.x(), grid.y())
    geometry.ice_thickness.get_from_proc0(H.get())

    geometry.bed_elevation.set(0.0)
    geometry.sea_level_elevation.set(0.0)
    geometry.ice_area_specific_volume.set(0.0)
    geometry.cell_area.set(grid.dx() * grid.dy())

    # compute surface elevation from ice thickness and bed elevation
    geometry.ensure_consistency(0.0)

    model.init(geometry)
    model.update(geometry, 0, 1)

    return model.mean_precipitation().numpy()

def comparison_test():

    g = grid()

    P_python = ltop_python(g)

    P_pism = ltop_pism(g)

    return P_python, P_pism
