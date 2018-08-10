#!/usr/bin/env python

import PISM
from math import cos, pi
import numpy as np

ctx = PISM.Context()
config = ctx.config

class Constants(object):
    "Spatially-constant inputs of the model"

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
    """
    Calculates orographic precipitation following Smith & Barstad (2004).

    """
    logger.info('Computing orographic precipitation')

    eps = 1e-18
    pad = 250

    ny, nx = orography.shape
    logger.debug('Raster shape before padding ({},{})'.format(nx, ny))

    padded_orography = np.pad(orography, pad, 'constant')
    nx, ny = padded_orography.shape
    logger.debug('Raster shape after padding ({},{})'.format(ny, nx))

    logger.info('Fourier transform orography')
    padded_orography_fft = np.fft.fft2(padded_orography)

    x_n_value = np.fft.fftfreq(ny, (1.0 / ny))
    y_n_value = np.fft.fftfreq(nx, (1.0 / nx))

    x_len = nx * dx
    y_len = ny * dy
    kx_line = np.divide(np.multiply(2.0 * np.pi, x_n_value), x_len)
    ky_line = np.divide(
        np.multiply(
            2.0 * np.pi,
            y_n_value),
        y_len)[
        np.newaxis].T

    kx = np.tile(kx_line, (nx, 1))
    ky = np.tile(ky_line, (1, ny))

    # Intrinsic frequency sigma = U*k + V*l
    u0 = constants.u
    v0 = constants.v

    logger.info('Calculate sigma')
    sigma = np.add(np.multiply(kx, u0), np.multiply(ky, v0))
    sigma_sqr_reg = sigma ** 2
    m_denom = np.power(sigma, 2.) - constants.f**2

    sigma_sqr_reg[
        np.logical_and(
            np.fabs(sigma_sqr_reg) < eps,
            np.fabs(
                sigma_sqr_reg >= 0))] = eps
    sigma_sqr_reg[
        np.logical_and(
            np.fabs(sigma_sqr_reg) < eps,
            np.fabs(
                sigma_sqr_reg < 0))] = -eps

    # The vertical wave number
    # Eqn. 12
    # Regularization
    m_denom[
        np.logical_and(
            np.fabs(m_denom) < eps,
            np.fabs(m_denom) >= 0)] = eps
    m_denom[
        np.logical_and(
            np.fabs(m_denom) < eps,
            np.fabs(m_denom) < 0)] = -eps

    m1 = np.divide(
        np.subtract(
            constants.Nm**2,
            np.power(
                sigma,
                2.)),
        m_denom)
    m2 = np.add(np.power(kx, 2.), np.power(ky, 2.))
    m_sqr = np.multiply(m1, m2)
    logger.info('Calculating m')
    m = np.sqrt(-1 * m_sqr)
    # Regularization
    m[np.logical_and(m_sqr >= 0, sigma == 0)] = np.sqrt(
        m_sqr[np.logical_and(m_sqr >= 0, sigma == 0)])
    m[np.logical_and(m_sqr >= 0, sigma != 0)] = np.sqrt(m_sqr[np.logical_and(
        m_sqr >= 0, sigma != 0)]) * np.sign(sigma[np.logical_and(m_sqr >= 0, sigma != 0)])
    # Numerator in Eqn. 49
    P_karot_num = np.multiply(np.multiply(np.multiply(
        constants.Cw, 1j), sigma), padded_orography_fft)
    P_karot_denom_Hw = np.subtract(1, np.multiply(
        np.multiply(constants.Hw, m), 1j))
    P_karot_denom_tauc = np.add(1, np.multiply(np.multiply(
        sigma, constants.tau_c), 1j))
    P_karot_denom_tauf = np.add(1, np.multiply(np.multiply(
        sigma, constants.tau_f), 1j))
    # Denominator in Eqn. 49
    P_karot_denom = np.multiply(
        P_karot_denom_Hw, np.multiply(
            P_karot_denom_tauc, P_karot_denom_tauf))
    P_karot = np.divide(P_karot_num, P_karot_denom)

    # Converting from wave domain back to space domain
    logger.info('Performing inverse Fourier transform')
    P = np.fft.ifft2(P_karot)
    spy = 31556925.9747
    logger.info('De-pad array')
    P = P[pad:-pad, pad:-pad]
    P = np.multiply(np.real(P), 3600)   # mm hr-1
    # Add background precip
    P0 = constants.P0
    logger.info('Adding background precpipitation {} mm hr-1'.format(P0))
    P += P0
    # Truncation

    if truncate:
        logger.info('Truncate precipitation')
        P[P < 0] = 0
    P_scale = constants.P_scale
    logger.info('Scale precipitation P = P * {}'.format(P_scale))
    P *= P_scale

    return sigma, m, P

def gaussian_bump(xmin, xmax, ymin, ymax, dx, dy, h_max=500.0,
                  x0=-25e3, y0=0.0, sigma_x=15e3, sigma_y=15e3):
    "Create the setup needed to reproduce Fig 4c in SB2004"
    # Reproduce Fig 4c in SB2004
    x = np.arange(xmin, xmax, dx)
    y = np.arange(ymin, ymax, dy)
    X, Y = np.meshgrid(x, y)
    Orography = h_max * np.exp(-(((X - x0)**2 / (2 * sigma_x**2)) +
                                 ((Y - y0)**2 / (2 * sigma_y**2))))
    return X, Y, Orography


def initialize_thickness(thickness, H):
    grid = thickness.grid()
    with PISM.vec.Access(nocomm=[thickness]):
        for (i, j) in grid.points():
            thickness[i, j] = orography[i, j]

def allocate(grid):
    H = PISM.model.createIceThicknessVec(grid)
    bed = PISM.model.createBedrockElevationVec(grid)
    sea_level = PISM.IceModelVec2S(grid, "sea_level", PISM.WITHOUT_GHOSTS)

    return H, bed, sea_level

def create_grid():
    P = PISM.GridParameters(config)
    P.horizontal_size_from_options()
    P.horizontal_extent_from_options()
    P.vertical_grid_from_options(config)
    P.ownership_ranges_from_options(ctx.size)

    return PISM.IceGrid(ctx.ctx, P)

def run(plot, pause, save):

    # set grid defaults
    config.set_double("grid.Mx", 400)
    config.set_double("grid.My", 400)

    config.set_double("grid.Lx", 300e3)
    config.set_double("grid.Ly", 200e3)

    config.set_double("grid.Mz", 2)
    config.set_double("grid.Lz", 1000)

    config.set_string("atmosphere.orographic_precipitation.file", "~/pism-olympics/data_sets/climate_forcing/ltop_climate_olympics_1000m_dir_220_kg_m-2_yr-1.nc")

    config.set_double("atmosphere.orographic_precipitation.background_precip_post", 1)
    
    grid = create_grid()

    thickness, bed, sea_level = allocate(grid)

    # set initial geometry
    bed.set(0.0)
    thickness.set(0.0)
    sea_level.set(0.0)
    initialize_thickness(thickness, orography)
    g = PISM.Geometry(grid)
    op = PISM.AtmosphereOrographicPrecipitation(grid)
    op.init(g)
    op.update(g, 0.1, 0.1)
    p = op.mean_precipitation()
    p.dump('foo.nc')

if __name__ == "__main__":

    x_min   = -100e3
    x_max   = 200e3
    y_min   = -150e3
    y_max   = 150e3
    dx      = 750
    dy      = 750
    x0      = -25e3
    y0      = 0.0
    sigma_x = 15e3
    sigma_y = 15e3

    _, _, orography = gaussian_bump(x_min, x_max, y_min, y_max, dx, dy,
                                   x0=x0, y0=y0, sigma_x=sigma_x, sigma_y=sigma_y)

    plot = PISM.OptionBool("-plot", "Plot bed elevation and uplift.")
    save = PISM.OptionBool("-save", "Save final states of the bed elevation and uplift.")
    pause = PISM.OptionBool("-pause", "Pause for 5 seconds to look at runtime 2D plots.")

    run(plot, pause, save)

