#!/usr/bin/env python

# Computes and plots exact solution "test M", in preparation for implementing C
# version for PISM verification.  Preliminary tests are possible with gridded
# version of this output, even before a C implementation.

from pylab import *
from scipy.optimize import fsolve
from scipy.integrate import odeint

SperA = 31556926.0
g = 9.81
rho = 910.0    # density of ice; kg/m^3
rhow = 1028.0   # density of ocean water; kg/m^3
barB = 3.7e8    # strength of shelf; Pa s^(1/3); as in Schoof 2006; compare 1.9e8
#   from MacAyeal et al 1996
Rc = 600.0e3  # calving front at 600 km
Rg = 300.0e3  # grounding line at 300 km
H0 = 500.0    # uniform 500 m thickness
ug = 100.0 / SperA  # velocity across grounding line is 100 m/year

# compute physical constant in ODE
Q = (1.0 - rho / rhow) * rho * g * Rc * H0 / (2.0 * barB)
print('physical constant: Q = %f' % Q)

# solve FF(x)=0 to get strain rate alpha'; x is RHS of ODE


def FF(x, alpha, r):
    DD = x * x + x * (alpha / r) + (alpha / r) ** 2
    return Q * DD ** (1. / 3.) - 2.0 * r * x - alpha


# ODE we are solving is   alpha' = GG(alpha,r)
qual = []  # quality information goes here


def GG(alpha, r):
    # heuristic: guess is about 1/7 th of solution to a nearby problem
    guess = 0.15 * ((Q / r) ** 3 - alpha[0] / r)
    result = fsolve(FF, guess, args=(alpha[0], r))
    qual.append([r, guess, result, abs(guess - result) / abs(result)])
    return [result]


# build vector of locations to solve (and plot)
dr = 1000.0  # plot on 1 km grid
r = linspace(Rg, Rc, ((Rc - Rg) / dr) + 1)

# solve:
print('solving with odeint from scipy.integrate; it reports: ', end=' ')
alpha = odeint(GG, [ug], r, printmessg=1,
               rtol=0.0,            # ask for no change in digits
               atol=0.00001 / SperA)  # ask for abs tol of 0.01 mm/year

# info to evaluate solution
qual = array(qual)
rused = qual[:, 0]
strainrate = qual[:, 2]
guesserr = qual[:, 3]
print('maximum relative error in initial guess for fsolve() is %f'
      % max(guesserr))

print("at calving front:  alpha(Rc) = %f  (m/year)"
      % (alpha[-1] * SperA))
# print "   (last r used = %f km; strain rate at last r = %f  (1/a))" \
#   % (rused[-1] / 1000.0, strainrate[-1] * SperA)

# plot velocity solution alpha(r), and also strain rate alpha'(r)
figure(1)
subplot(211)
plot(r / 1000.0, alpha * SperA, 'k', linewidth=3)
ylabel(r'velocity (m/year)', size=14)
subplot(212)
plot(rused / 1000.0, strainrate * SperA, 'ko-')
axis([Rg / 1000.0, Rc / 1000.0, 0, 1.1 * max(strainrate * SperA)])
xlabel(r'r (km)', size=14)
ylabel(r'strain rate (1/a)', size=14)
# print "saving figure 'combinedM.png'"  # optional: save as PNG
#savefig('combinedM.png', dpi=300, facecolor='w', edgecolor='w')
print('close figure to end')
show()
