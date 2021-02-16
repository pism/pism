import sympy
from sympy import S

from blatter import x, y, z, B, source_term, eta, M
from blatter_codegen import define, declare

sympy.var("R_0 H_0 rho_i g C_0 C_1 C_2", positive=True)
h = sympy.Function("h", positive=True)(x)

nx, ny, nz = sympy.var("n_(x:z)")

N = sympy.Matrix([nx, ny, nz])

# Glen exponents n
n = 3

def parameters(H_0, R_0, rho_i, g, B):
    # s = 1 corresponds to t = t_0
    s = 1

    c0 = H_0 * s
    c1 = s / R_0
    c2 = 2 * B**(-3) * (rho_i * g)**3 / 4

    return {C_0 : c0, C_1 : c1, C_2 : c2}

def H(x):
    return C_0 * (1 - (C_1 * x)**S("4/3"))**S("3/7")

def u_exact():
    """X-Z verification test using the Halfar dome geometry.
    """

    u0 = -C_2 * (h**4 - (h - z)**4) * h.diff(x)**3
    v0 = S(0)

    return u0, v0

def lateral_bc(u0, v0):
    N_right = {nx: 1, ny: 0, nz : 0}
    return (2 * eta(u0, v0, n) * M(u0, v0).row(0) * N)[0].subs(N_right), 0.0

def print_code(header=False):
    constants = ["H_0", "R_0", "rho_i", "g", "B"]
    coords = ["x", "z"]
    if header:
        declare(name="blatter_xz_halfar_exact", args=coords + constants)
        declare(name="blatter_xz_halfar_source", args=coords + constants)
        declare(name="blatter_xz_halfar_lateral", args=coords + constants)
        return

    definitions = parameters(H_0, R_0, rho_i, g, B)

    u0, v0 = u_exact()
    u0 = u0.subs(h, H(x)).subs(definitions).doit()
    define(u0, v0, name="blatter_xz_halfar_exact", args=coords + constants)

    f_u, f_v = source_term(eta(u0, v0, n), u0, v0)
    f_u = f_u.subs(definitions).doit()
    define(f_u, f_v, name="blatter_xz_halfar_source", args=coords + constants)

    f_u, f_v = lateral_bc(u0, v0)
    f_u = f_u.subs(definitions).doit()
    define(f_u, f_v, name="blatter_xz_halfar_lateral", args=coords + constants)
