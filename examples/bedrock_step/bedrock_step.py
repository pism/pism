#!/usr/bin/env python3
"""Steady state geometry from

A. H. Jarosch, C. G. Schoof, and F. S. Anslow, “Restoring mass
conservation to shallow ice flow models over complex terrain,” The
Cryosphere, vol. 7, Art. no. 1, Feb. 2013.

Most of the code below is borrowed from the supplement

http://www.the-cryosphere.net/7/229/2013/tc-7-229-2013-supplement.zip

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>

"""

import numpy as np

A      = 1e-16                  # Pa^(-3) / year
b_0    = 500.0                  # m
g      = 9.81                   # m / s^2
mdot_0 = 2.0                    # m / year
n      = 3.0                    # no units
rho    = 910.0                  # kg / m^3
x_m    = 20000.0                # m
x_s    = 7000.0                 # m

def bed_elevation(x):
    """Bed elevation in meters"""
    B = np.zeros_like(x)
    B[x < x_s] = b_0

    return B

def surface(x):
    """Steady state solution from section 6, in meters"""
    s_x_s_x_m = s_eval_x_s_x_m(x, x_s, x_m, n, A, mdot_0, rho, g)

    # combine solutions
    s = s_eval_x_s(x, x_s, x_m, n, A, mdot_0, rho, g, b_0) + b_0
    s[x >= x_s] = s_x_s_x_m[x >= x_s]

    # correct s
    s[x > x_m] = 0.

    return s

def accumulation(x):
    "Accumulation/ablation rate in m/year."
    # Eq. 54
    mdot = ((n*mdot_0)/(x_m**(2.*n-1.)))*x**(n-1.)*(abs(x_m-x)**(n-1.))*(x_m-2.*x)

    mdot[x>x_m] = 0.0

    return mdot

def s_eval_x_s_x_m(x,x_s,x_m,n,A,mdot_0,rho,g):
    "Surface elevation on [x_s, x_m], in meters"
    # Eq. 56
    s_x_s_x_m = (((2.*n+2.)*(n+2.)**(1./n)*mdot_0**(1./n))/(2.**(1./n)*6*n*A**(1./n)*rho*g*x_m**((2.*n-1)/n))*(x_m+2.*x)*(x_m-x)**2.)**(n/(2.*n+2.))

    return s_x_s_x_m

def s_eval_x_s(x,x_s,x_m,n,A,mdot_0,rho,g,b_0):
    "Surface elevation on [0, x_s), in meters"
    # Eq. 58
    h_splus = (((2.*n+2.)*(n+2.)**(1./n)*mdot_0**(1./n))/(2.**(1./n)*6*n*A**(1./n)*rho*g*x_m**((2.*n-1.)/n))*(x_m+2.*x_s)*(x_m-x_s)**2.)**(n/(2.*n+2.))

    # Eq. 59
    h_sminus = np.maximum(h_splus - b_0, 0.0)

    # Eq. 57
    return (h_sminus**((2.*n+2.)/n)-h_splus**((2.*n+2.)/n)+((2.*n+2.)*(n+2.)**(1./n)*mdot_0**(1./n))/(2.**(1./n)*6*n*A**(1./n)*rho*g*x_m**((2.*n-1.)/n))*(x_m+2.*x)*(x_m-x)**2.)**(n/(2.*n+2.))

def create_pism_input(filename):
    "Create a NetCDF file that can be used with PISM"
    import netCDF4 as NC

    def tile(x):
        return np.tile(x, [3, 1])

    dx = 200.0
    x = np.arange(-2 * x_m, 2 * x_m + dx, dx)

    with NC.Dataset(filename, "w") as f:
        f.createDimension("x", len(x))
        f.createDimension("y", 3)

        xv = f.createVariable("x", np.float64, ('x', ))
        xv.units = "meter"
        xv[:] = x

        yv = f.createVariable("y", np.float64, ('y', ))
        yv.units = "meter"
        yv[:] = [-100*dx, 0, 100*dx]

        b = f.createVariable("bed_elevation", np.float64, ('y', 'x'))
        b.units = "meter"
        b.standard_name = "bedrock_altitude"
        b[:] = tile(bed_elevation(np.abs(x)))

        s = f.createVariable("surface_elevation", np.float64, ('y', 'x'))
        s.units = "meter"
        s.standard_name = "surface_altitude"
        s[:] = tile(surface(np.abs(x)))

        ice_density = 910.0
        M = f.createVariable("climatic_mass_balance", np.float64, ('y', 'x'))
        M.units = "kg m-2 year-1"
        M.standard_name = "land_ice_surface_specific_mass_balance_flux"
        M[:] = tile(accumulation(np.abs(x))) * ice_density

        T = f.createVariable("ice_surface_temp", np.float64, ('y', 'x'))
        T.units = "Celsius"
        T.long_name = "ice surface temperature"
        T[:] = 0.0

if __name__ == "__main__":
    import sys
    create_pism_input(sys.argv[1])
