#!/usr/bin/env python

# Computes and plots exact solution for "test V", in preparation for
# implementing C version for PISM verification.

from pylab import *

secpera = 3.15569259747e7               # seconds per year
rho_sw = 1028                           # sea water density
rho_ice = 910                           # ice density
standard_gravity = 9.81                 # g
B0 = 1.9e8                              # ice hardness

# "typical constant ice parameter" as defined in the paper and in Van der
# Veen's "Fundamentals of Glacier Dynamics", 1999
C0 = (rho_ice * standard_gravity * (1.0 - rho_ice/rho_sw) / (4 * B0))**3

C = 2.45e-18                            # from the paper; does not match

# upstream ice thickness
H0 = 600.0                              # meters
# upstream ice velocity
v0 = 300.0/secpera                      # 300 meters/year
# upstream ice flux
Q0 = H0 * v0;

Mx = 201
x = linspace(0, 400e3, Mx)

def H(x):
    """Ice thickness."""    
    return (4 * C / Q0 * x + 1 / H0**4)**(-0.25)

def v(x):
    """Ice velocity."""    
    return Q0 / H(x)

def x_c(t):
    """Location of the calving front."""
    return Q0 / (4*C) * ((3*C*t + 1/H0**3)**(4.0/3.0) - 1/H0**4)

def plot_xc(t_years):
    """Plot the location of the calving front."""
    x = x_c(t_years * secpera)/1000.0   # convert to km
    a = axis()
    y_min = a[2]
    y_max = a[3]

    old_hold = hold(True)
    plot([x, x], [y_min, y_max], '--g')
    hold(old_hold)

figure(1)
subplot(211)
plot(x/1000, H(x))
plot_xc(300)
ylabel("m")
title("Ice thickness")
grid(True)

subplot(212)
plot(x/1000, v(x) * secpera)
plot_xc(300)
xlabel("km")
ylabel("m/year")
title("Horizontal ice velocity")
grid(True)

show()
