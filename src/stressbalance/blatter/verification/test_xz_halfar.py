import sympy as sp
from sympy import S

from blatter import x, y, z, B, source_term, eta, M
from blatter_codegen import define, declare, print_header, print_var, print_footer

sp.var("R_0 H_0 rho_i g C_0 C_1 C_2", positive=True)
h = sp.Function("h", positive=True)(x)

u = sp.Function("u")(x, z)
v = S(0)
u_y = S(0)
sp.var("u_x u_xx u_z u_xz u_zz h0 h_x h_xx")
subs = {h.diff(x, 2): h_xx,
        h.diff(x): h_x,
        h: h0,
        u.diff(x, 2): u_xx,
        u.diff(x): u_x,
        u.diff(x).diff(z): u_xz,
        u.diff(y): u_y,
        u.diff(z): u_z,
        u.diff(z, 2): u_zz}

nx, ny, nz = sp.var("n_(x:z)")

N = sp.Matrix([nx, ny, nz])

# Glen exponents n
n = 3

def constants():
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

def surface_bc(u0, v0, surface):
    ds = surface.diff(x)
    # normalized x component of the downward-pointing normal vector
    n_s_norm = (ds**2 + 1**2)**S("1/2")
    nx_s = - ds / n_s_norm
    ny_s = 0
    nz_s = 1 / n_s_norm

    N_surface = {nx: nx_s, ny: ny_s, nz : nz_s}

    return (2 * eta(u0, v0, n) * M(u0, v0).row(0) * N)[0].subs(N_surface)

def lateral_bc(u0, v0):
    N_right = {nx: 1, ny: 0, nz : 0}
    return (2 * eta(u0, v0, n) * M(u0, v0).row(0) * N)[0].subs(N_right), 0.0

def print_code(header=False):
    constants = ["H_0", "R_0", "rho_i", "g", "B"]
    coords = ["x", "z"]

    if header:
        declare(name="blatter_xz_halfar_exact", args=coords + constants)
        declare(name="blatter_xz_halfar_source", args=coords + constants)
        declare(name="blatter_xz_halfar_source_lateral", args=coords + constants)
        declare(name="blatter_xz_halfar_source_surface", args=["x"] + constants)
        return

    print_exact(coords + constants)
    print_source(coords + constants)
    print_source_lateral(coords + constants)
    print_source_surface(["x"] + constants)

def print_source_surface(args):
    "Print the code computing the extra term at the top surface"

    f_top = surface_bc(u, v, h)

    # take advantage of the fact that u_z = 0 at z = h
    f_top = f_top.subs(subs).subs(u_z, 0)

    u0, _ = u_exact()

    U_x = u0.diff(x).subs(subs)
    U_z = u0.diff(z).subs(subs)

    print_header("blatter_xz_halfar_source_surface", args)

    for key, value in constants().items():
        print_var(value, key)
    print_var(H(x), h0)
    print_var(h0, z)
    print_var(H(x).diff(x), h_x)
    print_var(H(x).diff(x, 2), h_xx)
    print_var(U_x, u_x)

    print_footer(f_top, 0.0)

def print_source_lateral(args):
    "Print the code computing the extra term at the right boundary"

    f_lat, _ = lateral_bc(u, v)
    f_lat = f_lat.subs(subs).factor()

    u0, _ = u_exact()

    U_x = u0.diff(x).subs(subs)
    U_z = u0.diff(z).subs(subs)

    print_header("blatter_xz_halfar_source_lateral", args)

    for key, value in constants().items():
        print_var(value, key)
    print_var(H(x), h0)
    print_var(H(x).diff(x), h_x)
    print_var(H(x).diff(x, 2), h_xx)
    print_var(U_x, u_x)
    print_var(U_z, u_z)

    print_footer(f_lat, 0.0)

def print_exact(args):
    u0, v0 = u_exact()

    u0 = u0.subs(subs)

    print_header("blatter_xz_halfar_exact", args)

    for key, value in constants().items():
        print_var(value, key)
    print_var(H(x), h0)
    print_var(H(x).diff(x), h_x)

    print_footer(u0, v0)

def print_source(args):
    f, _ = source_term(eta(u, v, 3), u, v)
    f = f.subs(subs)

    u0, _ = u_exact()

    U_x = u0.diff(x).subs(subs)
    U_z = u0.diff(z).subs(subs)

    print_header("blatter_xz_halfar_source", args)

    for key, value in constants().items():
        print_var(value, key)
    print_var(H(x), h0)
    print_var(H(x).diff(x), h_x)
    print_var(H(x).diff(x, 2), h_xx)
    print_var(U_x, u_x)
    print_var(U_z, u_z)
    print_var(U_x.diff(x), u_xx)
    print_var(U_x.diff(z), u_xz)
    print_var(U_z.diff(z), u_zz)

    print_footer(f, 0.0)
