import sympy
from sympy import sin, cos, pi, S

from blatter import x, y, z, B, M, source_term, eta
from manufactured_solutions import define, declare

sympy.var("rho_i, rho_w, g, L")

nx, ny, nz = sympy.var("n_(x:z)")

N = sympy.Matrix([nx, ny, nz])

# Glen exponent of 1 (constant viscosity)
n = 1

def exact():
    """X-Z verification test for lateral boundary conditions

    Constant bed elevation and ice thickness, constant hardness, constant viscosity.

    """

    u0 = (rho_i - rho_w) * g * L / (2 * B * pi) * sin(pi * x / L) * z
    v0 = S(0)

    return u0, v0

def surface_bc():
    u0, v0 = exact()

    return (eta(u0, v0, n) * M(u0, v0).row(0) * N)[0].subs({nx: 0, ny: 0, nz : 1})

def lateral_bc():

    u0, v0 = exact()

    return (eta(u0, v0, n) * M(u0, v0).row(0) * N)[0].subs({nx: 1, ny: 0, nz : 0, x : L})

def print_code(header=False):
    args = ["x", "z", "B", "L", "rho_i", "rho_w"]
    if header:
        declare(name="blatter_xz_cfbc_exact", args=args)
        declare(name="blatter_xz_cfbc_source", args=args)
        declare(name="blatter_xz_cfbc_surface", args=args)
        declare(name="blatter_xz_cfbc_base", args=args)
        return

    u0, v0 = exact()
    define(u0, v0, name="blatter_xz_cfbc_exact", args=args)

    f_u, f_v = source_term(eta(u0, v0, n), u0, v0)
    define(f_u, f_v, name="blatter_xz_cfbc_source", args=args)

    f_s = surface_bc()
    define(f_s, S(0), name="blatter_xz_cfbc_surface", args=args)
    define(-f_s, S(0), name="blatter_xz_cfbc_base", args=args)
