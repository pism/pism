#!/usr/bin/env python
import numpy as np

"""This module contains MISMIP constants and parameters, as well as functions
computing theoretical steady state profiles corresponding to various MISMIP
experiments.

It should not be cluttered with plotting or NetCDF output code.
"""


def secpera():
    "Number of seconds per year."
    return 3.15569259747e7


def L():
    "The length of the MISMIP domain."
    return 1800e3


def N(mode):
    "Number of grid spaces corresponding to a MISMIP 'mode.'"
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
    # Values of A to be used in experiments 1 and 2.

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


def run_length(experiment, step):
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

    return 3e4


def rho_i():
    "Ice density"
    return 900.0


def rho_w():
    "Water density"
    return 1000.0


def g():
    """Acceleration due to gravity. (Table 2 on page 19 of mismip_4.pdf
    uses this value, i.e. g = 9.8 m s-2.)"""
    return 9.8


def n():
    "Glen exponent"
    return 3.0


def a():
    "Accumulation rate (m/s)"
    return 0.3 / secpera()


def m(experiment):
    "Sliding law exponent"
    if experiment in ("1a", "2a", "3a"):
        return 1 / 3.0

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
        return -(729. - 2184.8 * xx ** 2. + 1031.72 * xx ** 4. - 151.72 * xx ** 6.)

    raise ValueError("invalid experiment (%s)" % experiment)


def b_slope(experiment, x):
    """The x-derivative of b(experiment, x)."""

    if experiment in ("1a", "1b", "2a", "2b"):
        return 778.5 / 7.5e5

    if experiment in ("3a", "3b"):
        xx = x / 7.5e5
        return -(- 2184.8 * (2. / 7.5e5) * xx
                 + 1031.72 * (4. / 7.5e5) * xx ** 3.
                 - 151.72 * (6. / 7.5e5) * xx ** 5.)

    raise ValueError("invalid experiment (%s)" % experiment)


def cold_function(experiment, step, x, theta=0.0):
    """Evaluates function whose zeros define x_g in 'cold' steady marine sheet problem."""
    r = rho_i() / rho_w()
    h_f = r ** (-1.) * b(experiment, x)
    b_x = b_slope(experiment, x)
    s = a() * x
    rho_g = rho_i() * g()
    return (theta * a()
            + C(experiment) * s ** (m(experiment) + 1.0) / (rho_g * h_f ** (m(experiment) + 2.))
            - theta * s * b_x / h_f
            - A(experiment, step) * (rho_g * (1.0 - r) / 4.0) ** n() * h_f ** (n() + 1.0))


def x_g(experiment, step, theta=0.0):
    """Computes the theoretical grounding line location using Newton's method."""

    # set the initial guess
    if experiment in ("3a", "3b"):
        x = 800.0e3
    else:
        x = 1270.0e3

    delta_x = 10.  # Finite difference step size (metres) for gradient calculation
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
        normf = abs(f)
        grad = (F(x + delta_x) - f) / delta_x
        dx = -f / grad
        x = x + dx

    return x


def thickness(experiment, step, x, theta=0.0):
    """Compute ice thickness for x > 0.
    """
    # compute the grounding line position
    xg = x_g(experiment, step, theta)

    def surface(h, x):
        b_x = b_slope(experiment, x)
        rho_g = rho_i() * g()
        s = a() * np.abs(x)
        return b_x - (C(experiment) / rho_g) * s ** m(experiment) / h ** (m(experiment) + 1)

    # extract the grounded part of the grid
    x_grounded = x[x < xg]

    # We will integrate from the grounding line inland. odeint requires that
    # the first point in x_grid be the one corresponding to the initial
    # condition; append it and reverse the order.
    x_grid = np.append(xg, x_grounded[::-1])

    # use thickness at the grounding line as the initial condition
    h_f = b(experiment, xg) * rho_w() / rho_i()

    import scipy.integrate
    thk_grounded = scipy.integrate.odeint(surface, [h_f], x_grid, atol=1.e-9, rtol=1.e-9)

    # now 'result' contains thickness in reverse order, including the grounding
    # line point (which is not on the grid); discard it and reverse the order.
    thk_grounded = np.squeeze(thk_grounded)[:0:-1]

    # extract the floating part of the grid
    x_floating = x[x >= xg]

    # compute the flux through the grounding line
    q_0 = a() * xg

    # Calculate ice thickness for shelf from van der Veen (1986)
    r = rho_i() / rho_w()
    rho_g = rho_i() * g()
    numer = h_f * (q_0 + a() * (x_floating - xg))
    base = q_0 ** (n() + 1) + h_f ** (n() + 1) * ((1 - r) * rho_g / 4) ** n() * A(experiment, step) \
        * ((q_0 + a() * (x_floating - xg)) ** (n() + 1) - q_0 ** (n() + 1)) / a()
    thk_floating = numer / (base ** (1.0 / (n() + 1)))

    return np.r_[thk_grounded, thk_floating]


def plot_profile(experiment, step, out_file):
    from pylab import figure, subplot, hold, plot, xlabel, ylabel, text, title, axis, vlines, savefig

    if out_file is None:
        out_file = "MISMIP_%s_A%d.pdf" % (experiment, step)

    xg = x_g(experiment, step)

    x = np.linspace(0, L(), N(2) + 1)
    thk = thickness(experiment, step, x)
    x_grounded, thk_grounded = x[x < xg],  thk[x < xg]
    x_floating, thk_floating = x[x >= xg], thk[x >= xg]

    figure(1)
    ax = subplot(111)
    hold(True)
    plot(x / 1e3, np.zeros_like(x), ls='dotted', color='red')
    plot(x / 1e3, -b(experiment, x), color='black')
    plot(x / 1e3, np.r_[thk_grounded - b(experiment, x_grounded),
                        thk_floating * (1 - rho_i() / rho_w())],
         color='blue')
    plot(x_floating / 1e3, -thk_floating * (rho_i() / rho_w()), color='blue')
    _, _, ymin, ymax = axis(xmin=0, xmax=x.max() / 1e3)
    vlines(xg / 1e3, ymin, ymax, linestyles='dashed', color='black')

    xlabel('distance from the summit, km')
    ylabel('elevation, m')
    text(0.6, 0.9, "$x_g$ (theory) = %4.0f km" % (xg / 1e3),
         color='black', transform=ax.transAxes)
    title("MISMIP experiment %s, step %d" % (experiment, step))
    savefig(out_file)


if __name__ == "__main__":
    from optparse import OptionParser

    parser = OptionParser()

    parser.usage = "%prog [options]"
    parser.description = "Plots the theoretical geometry profile corresponding to MISMIP experiment and step."
    parser.add_option("-e", "--experiment", dest="experiment", type="string",
                      default='1a',
                      help="MISMIP experiments (one of '1a', '1b', '2a', '2b', '3a', '3b')")
    parser.add_option("-s", "--step", dest="step", type="int", default=1,
                      help="MISMIP step number")
    parser.add_option("-o", dest="out_file", help="output file name")

    (opts, args) = parser.parse_args()

    plot_profile(opts.experiment, opts.step, opts.out_file)
