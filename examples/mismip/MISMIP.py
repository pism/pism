#!/usr/bin/env python
import numpy as np
from scipy.integrate import odeint

"""This module contains MISMIP constants and parameters, as well as functions
computing theoretical steady state profiles corresponding to various MISMIP
experiments.

It should not be cluttered with plotting or NetCDF output code.
"""

def L():
    "The length of the MISMIP domain."
    return 1800e3

def N(mode):
    "Number of grid points corresponding to a MISMIP 'mode.'"
    if mode == 1:
        return 150

    if mode == 2:
        return 1500

    raise ValueError("invalid mode (%s)" % mode)

def A(experiment, step):
    """Ice softness parameter for given experiment and step."""
    A1 = np.array([4.6416e-24,  2.1544e-24,  1.0e-24,
                   4.6416e-25,  2.1544e-25,  1.0e-25,
                   4.6416e-26,  2.1544e-26,  1.0e-26])
    # Values of A to be used in experiment 1.

    A3a = np.array([3.0e-25, 2.5e-25, 2.0e-25,
                    1.5e-25, 1.0e-25, 5.0e-26,
                    2.5e-26, 5.0e-26, 1.0e-25,
                    1.5e-25, 2.0e-25, 2.5e-25,
                    3.0e-25])
    # Values of A to be used in experiment 3a.

    A3b = np.array([1.6e-24, 1.4e-24, 1.2e-24,
                    1.0e-24, 8.0e-25, 6.0e-25,
                    4.0e-25, 2.0e-25, 4.0e-25,
                    6.0e-25, 8.0e-25, 1.0e-24,
                    1.2e-24, 1.4e-24, 1.6e-24])
    # Values of A to be used in experiment 3b.

    try:
        if experiment in ("1a", "1b", "2a", "2b"):
            return A1[step - 1]

        if experiment == "3a":
            return A3a[step - 1]

        if experiment == "3b":
            return A3b[step - 1]
    except:
        raise ValueError("invalid step (%s) for experiment %s" % (step, experiment))

    raise ValueError("invalid experiment (%s)" % experiment)

def time_interval(experiment, step):
    """Returns the time interval for an experiment 3 step."""
    T3a = np.array([3.0e4, 1.5e4, 1.5e4,
                    1.5e4, 1.5e4, 3.0e4,
                    3.0e4, 1.5e4, 1.5e4,
                    3.0e4, 3.0e4, 3.0e4,
                    1.5e4])
    # Time intervals to be used in experiment 3a.

    T3b = np.array([3.0e4, 1.5e4, 1.5e4,
                    1.5e4, 1.5e4, 1.5e4,
                    1.5e4, 3.0e4, 1.5e4,
                    1.5e4, 1.5e4, 1.5e4,
                    1.5e4, 3.0e4, 1.5e4])
    # Time intervals to be used in experiment 3b.

    try:
        if experiment == "3a":
            return T3a[step - 1]

        if experiment == "3b":
            return T3b[step - 1]
    except:
        raise ValueError("invalid step (%s) for experiment %s" % (step, experiment))

    raise ValueError("invalid experiment (%s)" % experiment)

def rho_i():
    "Ice density"
    return 900.0

def rho_w():
    "Water density"
    return 1000.0

def g():
    "Acceleration due to gravity"
    return 9.81

def n():
    "Glen exponent"
    return 3.0

def a():
    "Accumulation rate (m/s)"
    secpera = 31556926.  # (= 365.2422 days/a)
    return 0.3 / secpera

def m(experiment):
    "Sliding law exponent"
    if experiment in ("1a", "2a", "3a"):
        return 1/3.0

    if experiment in ("1b", "2b", "3b"):
        return 1.0

    raise ValueError("invalid experiment (%s)" % experiment)

def C(experiment):
    "Sliding law coefficient"
    if experiment in ("1a", "2a", "3a"):
        return 7.624e6

    if experiment in ("1b", "2b", "3b"):
        return 7.2082e10

    raise ValueError("invalid experiment (%s)" % experiment)

def b(experiment, x):
    "Bed depth below sea level. (-b(x) = topg(x))"

    if experiment in ("1a", "1b", "2a", "2b"):
        return -720. + 778.5 * (x / 7.5e5)

    if experiment in ("3a", "3b"):
        xx = x / 7.5e5
        return -(729. - 2184.8 * xx**2. + 1031.72 * xx**4. - 151.72 * xx**6.)

    raise ValueError("invalid experiment (%s)" % experiment)

def b_slope(experiment, x):
    """The x-derivative of b(experiment, x)."""

    if experiment in ("1a", "1b", "2a", "2b"):
        return 778.5 / 7.5e5

    if experiment in ("3a", "3b"):
        xx = x / 7.5e5
        return -(- 2184.8 * (2./7.5e5) * xx
                 + 1031.72 * (4./7.5e5) * xx**3.
                 - 151.72 * (6./7.5e5) * xx**5.)

    raise ValueError("invalid experiment (%s)" % experiment)

def cold_function(experiment, step, x, theta=0.0):
    """Evaluates function whose zeros define x_g in 'cold' steady marine sheet problem."""
    r = rho_i() / rho_w()
    h_f = r**(-1.) * b(experiment, x)
    b_x = b_slope(experiment, x);
    s = a() * x;
    rho_g = rho_i() * g()
    return (theta * a()
            + C(experiment) * s**(m(experiment) + 1.0) / (rho_g * h_f**(m(experiment) + 2.))
            - theta * s * b_x / h_f
            - A(experiment, step) * (rho_g * (1.0 - r) / 4.0)**n() * h_f**(n() + 1.0))

def x_g(experiment, step, theta=0.0):
    """Computes the theoretical grounding line location using Newton's method."""

    # set the initial guess
    if experiment in ("3a", "3b"):
        x = 800.0e3
    else:
        x = 1270.0e3

    delta_x = 10. # Finite difference step size (metres) for gradient calculation
    tolf = 1.e-4  # Tolerance for finding zeros
    eps = np.finfo(float).eps
    normf = tolf + eps
    toldelta = 1.e1                     # Newton step size tolerance
    dx = toldelta + 1.0

    # this is just a shortcut
    def F(x):
        return cold_function(experiment, step, x, theta)

    while (normf > tolf) or (abs(dx) > toldelta):
        f = F(x)
        normf = abs(f);
        grad = (F(x + delta_x) - f) / delta_x
        dx = -f / grad
        x = x + dx
    print "x_g = %.3f km" % (x / 1.e3)

    return x

# def thickness(experiment, step, theta=0.0, x_grid):
#     xg = x_g(experiment, step, theta)

#     def surface(h, x):
#         b_x = b_slope(experiment, x)
#         rho_g = rho_i() * g()
#         s = a() * np.abs(x)
#         return b_x - (C(experiment) / rho_g) * s**m(experiment) / h**(m(experiment) + 1)
