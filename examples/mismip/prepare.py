#!/usr/bin/env python

try:
    from netCDF4 import Dataset as NC
except:
    from netCDF3 import Dataset as NC

import MISMIP

import numpy as np

def topg(x, experiment):
    """Computes bed elevation as a function of x.
    Experiment can be '1a', '1b', '2a', '2b', '3a', '3b'.
    """

    return np.tile(-MISMIP.b(experiment, np.abs(x)), (3, 1))

def smb(x):
    """Computes the surface mass balance."""
    return np.tile(np.zeros_like(x) + MISMIP.a(), (3, 1))

def ice_surface_temp(x):
    """Computes the ice surface temperature (irrelevant)."""
    return np.tile(np.zeros_like(x) + 273.15, (3, 1))

def x(mismip_mode, N=None):
    if mismip_mode in (1, 2):
        return np.linspace(-MISMIP.L(), MISMIP.L(), 2 * MISMIP.N(mismip_mode) + 1)

    return return np.linspace(-MISMIP.L(), MISMIP.L(), N)

def y(x):
    """Computes y coordinates giving the 1:1 aspect ratio.
    Takes cross-flow grid periodicity into account."""
    dx = x[1] - x[0]
    My = 3.0
    Ly = dx * My / 2.0
    dy = 2 * Ly / My
    return np.array([-dy, 0, dy])

