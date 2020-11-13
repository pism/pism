import sympy
from sympy import S

from blatter import x, y, z, B, source_term, eta
from manufactured_solutions import define, declare

sympy.var("R_0 H_0 t_0 t s rho_i g C_0 C_1 C_2", positive=True)
h = sympy.Function("h", positive=True)(x)

# Glen exponents n
n = 3

# s = 1 corresponds to t = t_0
s = 1

A = 3.1689e-24
constants = {H_0 : 1000, # ice thickness at the center at t_0
             R_0: 750e3, # radius at t_0
             s: 1,       # (t / t_0)^(-1/11), s = 1 means use t_0
             g: 9.81,    # acceleration due to gravity
             rho_i: 910, # ice density
             B: A**(-1.0/3.0)} # ice hardness

s0 = (t / t_0)**S("-1/11")

c_0 = H_0 * s
c_1 = s / R_0
c_2 = 2 * B**(-3) * (rho_i * g)**3 / 4

# Combine definitions to simplify substitutions
definitions = {C_0: c_0, C_1: c_1, C_2: c_2}

H = C_0 * (1 - (C_1 * x)**S("4/3"))**S("3/7")

def exact():
    """X-Z verification test for lateral boundary conditions

    Constant bed elevation and ice thickness, constant hardness, constant viscosity.

    """

    u0 = -C_2 * (h**4 - (h - z)**4) * h.diff(x)**3
    v0 = S(0)

    return u0, v0

def print_code(header=False):
    args = ["x", "z", "H_0", "R_0", "rho_i", "g", "B"]
    if header:
        declare(name="blatter_xz_halfar_exact", args=args)
        declare(name="blatter_xz_halfar_source", args=args)
        return

    u0, v0 = exact()
    u0 = u0.subs(h, H).subs(definitions).doit()
    define(u0, v0, name="blatter_xz_halfar_exact", args=args)

    f_u, f_v = source_term(eta(u0, v0, n), u0, v0)
    f_u = f_u.subs(definitions).doit()
    define(f_u, f_v, name="blatter_xz_halfar_source", args=args)
