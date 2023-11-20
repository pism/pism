import sympy
from blatter import B, M, eta, source_term, x, y, z
from blatter_codegen import declare, define
from sympy import S, cos, pi, sin

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

    return (2 * eta(u0, v0, n) * M(u0, v0).row(0) * N)[0].subs({nx: 0, ny: 0, nz: 1})


def lateral_bc():
    u0, v0 = exact()

    return (2 * eta(u0, v0, n) * M(u0, v0).row(0) * N)[0].subs(
        {nx: 1, ny: 0, nz: 0, x: L}
    )


def print_code(header=False):
    args = ["x", "z", "B", "L", "rho_i", "rho_w", "g"]
    source_args = ["x", "z", "L", "rho_i", "rho_w", "g"]
    surface_args = ["x", "L", "rho_i", "rho_w", "g"]
    if header:
        declare(name="blatter_xz_cfbc_exact", args=args)
        declare(name="blatter_xz_cfbc_source", args=source_args)
        declare(name="blatter_xz_cfbc_surface", args=surface_args)
        declare(name="blatter_xz_cfbc_base", args=surface_args)
        return

    u0, v0 = exact()
    define(u0, v0, name="blatter_xz_cfbc_exact", args=args)

    f_u, f_v = source_term(eta(u0, v0, n), u0, v0)
    define(f_u, f_v, name="blatter_xz_cfbc_source", args=source_args)

    f_s = surface_bc()
    define(f_s, S(0), name="blatter_xz_cfbc_surface", args=surface_args)
    define(-f_s, S(0), name="blatter_xz_cfbc_base", args=surface_args)
