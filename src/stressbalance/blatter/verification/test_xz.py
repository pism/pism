#!/usr/bin/env python3

import sympy as sp
from sympy.core import S

from blatter import x, y, z, s
from manufactured_solutions import define, declare

return_template = """
  return {{
    {},
    {}
  }};"""

# misc. variables
sp.var("A, s_0, rho, g, beta, mu, n, alpha, H", positive=True)
sp.var("phi_(1:6)")

surface = s_0 - alpha * x**2

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

def exact():
    s_x = s.diff(x)

    np = n + 1
    nm = n - 1

    rhog = rho * g

    u = (2 * A * rhog**n / np) * ((s - z)**np - H**np) * abs(s_x)**nm * s_x - rhog * H * s_x / beta
    v = S(0)

    return u, v

nx, ny, nz = sp.var("n_(x:z)")
N = sp.Matrix([nx, ny, nz])

# outward-pointing normals to the top and bottom surfaces
n_mag = sp.sqrt((2 * alpha * x)**2 + 1)
nx_top = 2 * alpha * x / n_mag
nz_top = 1 / n_mag
nx_bed = - nx_top
nz_bed = - nz_top

# extra term for vertical cliffs
f_lat = -4 * phi_4 * mu

# extra term for the top surface
f_top = (- 4 * phi_4 * mu * nx
         - 4 * phi_2 * x**2 * phi_1**3 * mu * nz)

# extra term for the bottom surface
f_bed = (- 4 * phi_4 * mu * nx
         - 4 * phi_2 * x**2 * phi_1**3 * mu * nz
         + 2 * H * alpha * rho * g * x - beta * x**2 * phi_2 * (phi_1**4 - H**4))

def print_var(var, name):
    print("  double " + sp.ccode(var, assign_to=name))

def print_exact(args):
    "Print the code computing the exact solution"
    arguments = ", ".join(["double " + x for x in args])

    U, _ = exact()
    U = U.subs({s: surface, n: 3}).doit()

    print("")
    print("Vector2 blatter_xz_exact({args}) {{".format(args=arguments))
    print(return_template.format(sp.ccode(U), 0.0))
    print("}")

def print_source(args, header=False):
    "Print the code computing the source term"
    arguments = ", ".join(["double " + x for x in args])

    print("")
    print("Vector2 blatter_xz_source({args}) {{".format(args=arguments))
    print_var(psi_1.subs(s, surface), "phi_1")
    print_var(psi_2, "phi_2")
    print_var(psi_3, "phi_3")
    print_var(psi_4, "phi_4")
    print_var(psi_5, "phi_5")
    print_var(mu_xz, "mu")
    print(return_template.format(sp.ccode(f_u), 0.0))
    print("}")

def print_source_bed(args, header=False):
    "Print the code computing the extra term at the base"
    arguments = ", ".join(["double " + x for x in args])

    print("")
    print("Vector2 blatter_xz_source_bed({args}) {{".format(args=arguments))
    print_var(psi_1.subs(s, surface), "phi_1")
    print_var(psi_2, "phi_2")
    print_var(psi_3, "phi_3")
    print_var(psi_4, "phi_4")
    print_var(mu_xz, "mu")
    print_var(nx_bed, nx)
    print_var(nz_bed, nz)
    print(return_template.format(sp.ccode(f_bed), 0.0))
    print("}")

def print_source_surface(args, header=False):
    "Print the code computing the extra term at the top surface"
    arguments = ", ".join(["double " + x for x in args])

    print("")
    print("Vector2 blatter_xz_source_surface({args}) {{".format(args=arguments))
    print_var(psi_1.subs(s, surface), "phi_1")
    print_var(psi_2, "phi_2")
    print_var(psi_3, "phi_3")
    print_var(psi_4, "phi_4")
    print_var(mu_xz, "mu")
    print_var(nx_top, nx)
    print_var(nz_top, nz)
    print(return_template.format(sp.ccode(f_top), 0.0))
    print("}")

def print_code(header=False):
    "Print all the code needed by the XZ verification test"
    args = ["x", "z", "A", "rho", "g", "s_0", "alpha", "H", "beta"]

    if header:
        declare(name="blatter_xz_exact", args=args)
        declare(name="blatter_xz_source", args=args)
        declare(name="blatter_xz_source_bed", args=args)
        declare(name="blatter_xz_source_surface", args=args)
        return

    print_exact(args)
    print_source(args)
    print_source_bed(args)
    print_source_surface(args)
