#!/usr/bin/env python3

import sympy as sp
from sympy import exp, sin, cos, pi
from sympy.core import S

from blatter import *

from manufactured_solutions import *

# ice softness
sp.var("A, s_0, rho, g, beta, mu")

nx, ny, nz = sp.var("n_(x:z)")

N = sp.Matrix([nx, ny, nz])

def surface():
    return s_0 - alpha * x**2

sp.var("phi_(1:6)")
# source term; see equation (37) in Tezaur et al
f_u = S(16) / 3 * A * mu**4 * (-2 * phi_4**2 * phi_5
                               + 24 * phi_3 * phi_4 * (phi_1 + 2 * alpha * x**2)
                               - 6 * x**3 * phi_1**3 * phi_2 * phi_3
                               - 18 * x**2 * phi_1**2 * phi_2 * phi_4**2
                               - 6 * x * phi_1 * phi_3 * phi_5)

# See equation (38) in Tezaur et al
psi_1 = z - s
psi_2 = 4 * A * alpha**3 * rho**3 * g**3 * x
psi_3 = 4 * x**3 * phi_1**5 * phi_2**2
psi_4 = 8 * alpha * x**3 * phi_1**3 * phi_2 - (2 * H * alpha * rho * g) / beta + 3 * x * phi_2 * (phi_1**4 - H**4)
psi_5 = (56 * alpha * x**2 * phi_1**3 * phi_2 + 48 * alpha**2 * x**4 * phi_1**2 * phi_2
         + 6 * phi_2 * (phi_1**4 - H**4))

mu_xz = S(1) / 2 * (A * phi_4**2 + A * x * phi_1 * phi_3)**(-S(1) / 3)

def exact_xz():

    s_x = s.diff(x)

    np = n + 1
    nm = n - 1

    rhog = rho * g

    u = (2 * A * rhog**n / np) * ((s - z)**np - H**np) * abs(s_x)**nm * s_x - rhog * H * s_x
    v = S(0)

    return u, v

def source_xz_albany():
    "The source term for the XZ test case as defined in the Albany/FELIX paper."

    return f_u, S(0)

f_lat = -4 * phi_4 * mu

f_top = -4 * phi_4 * mu * n_x - 4 * phi_2 * x**2 * x**2 * phi_1**3 * mu * n_y

def print_xz(header=False):

    u0, v0 = exact_xz()
    f_u, f_v = source_xz_albany()

    args = ["x", "z", "A", "rho", "g", "s_0", "alpha", "H"]
    print_exact(u0, v0, name="exact_xz", args=args, header=header)
    print_source(f_u, f_v, name="source_xz", args=args, header=header)

def main(header=False):
    if header:
        print("#include <cmath>")
        print("")
        print('#include "pism/util/Vector2.hh"')
    else:
        print('#include "manufactured_solutions.hh"')

    print("")
    print("namespace pism {")

    print_xz(header)

    print("")
    print("} // end of namespace pism")

if False and __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument("--header", dest="header", action="store_true",
                        help="print function declarations for the header file")

    options = parser.parse_args()

    main(options.header)
