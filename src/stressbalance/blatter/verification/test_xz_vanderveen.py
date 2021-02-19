import sympy as sp
from sympy import S, Eq, solve

from blatter import x, y, z, B, source_term, eta, M, grad
from blatter_codegen import define, declare, print_header, print_var, print_footer

sp.var("alpha H_0 Q_0 rho_i g C", positive=True)

x_p = sp.symbols("x_p", positive=True)

nx, ny, nz = sp.var("n_(x:z)")

N = sp.Matrix([nx, ny, nz])

H = sp.Function("H", positive=True)(x)
H_exact = (4 * C * x / Q_0 + H_0**(-4))**S("-1/4")

u = sp.Function("u", positive=True)(x)
u_exact = Q_0 / H

v = S(0)
v_exact = v

def C0():
    f_u, _ = source_term(eta(u_exact, v_exact, 3), u_exact, v_exact)

    # surface elevation
    s = alpha * H

    # "driving stress" term
    f_taud = rho_i * g * s.diff(x)

    equation = sp.Eq(f_u, f_taud).subs(H, H_exact).doit().subs(x, x_p)

    return sp.solve(equation, C)[0]

def source_lateral():
    "Lateral BC"
    N_right = {nx: 1, ny: 0, nz: 0}

    I = 2 * eta(u_exact, v_exact, 3) * M(u_exact, v_exact).row(0) * N

    I = I[0].subs(N_right).subs(H, H_exact).doit()

    # tell SymPy that x is positive to make it simplify the expression
    I = I.subs(x, x_p)

    I = I.subs(x_p, x).subs(H_exact, H)

    return I, 0

def basal_beta():
    # surface elevation
    s = alpha * H

    b = sp.Function("b", real=True)(x)
    nx_b, ny_b, nz_b = list(grad(b - z))

    db = (s - H).subs(H, H_exact).factor(fraction=False).diff(x).subs(H_exact, H)

    # normalized x component of the downward-pointing normal vector
    norm = (db**2 + 1**2)**S("1/2")
    nx_n = db / norm
    nz_n = -1 / norm

    N_base = {nx: nx_n, ny: 0, nz: nz_n}

    I_base = 2 * eta(u_exact, v_exact, 3) * M(u_exact, v_exact).row(0) * N
    I_base = I_base[0].subs(N_base).doit()

    BC_basal = I_base.subs(H, H_exact).doit().subs(H_exact, H)

    sp.var("beta", positive=True)
    eq_beta = sp.Eq(BC_basal + beta * u, 0)

    return solve(eq_beta, beta)[0].subs(u, Q_0 / H)

def surface_bc():
    # surface elevation
    s = alpha * H

    ds = s.diff(x)
    # normalized x component of the upward-pointing normal vector
    n_s_norm = (ds**2 + 1**2)**S("1/2")
    nx_s = - ds / n_s_norm
    ny_s = 0
    nz_s = 1 / n_s_norm

    N_surface = {nx: nx_s, ny: ny_s, nz : nz_s}

    f_s = (2 * eta(u_exact, v_exact, 3) * M(u_exact, v_exact).row(0) * N)

    f_s = f_s[0].subs(N_surface).subs(H, H_exact).doit()

    f_s = f_s.subs(H_exact, H).factor()

    return f_s, 0

def print_code(header=False):
    args = ["x", "alpha", "H_0", "Q_0", "rho_i", "g", "B"]

    if header:
        declare(name="blatter_xz_vanderveen_thickness", args=args, return_type="double")
        declare(name="blatter_xz_vanderveen_exact", args=args)
        declare(name="blatter_xz_vanderveen_source_lateral", args=args)
        declare(name="blatter_xz_vanderveen_source_surface", args=args)
        declare(name="blatter_xz_vanderveen_beta", args=args, return_type="double")
        return

    print_thickness(args)
    print_exact(args)
    print_source_lateral(args)
    print_source_surface(args)
    print_basal_beta(args)

def print_exact(args):

    print_header("blatter_xz_vanderveen_exact", args)

    h0 = sp.var("thickness")

    print_var(C0(), C)
    print_var(H_exact, h0)

    print_footer(u_exact.subs(H, h0), v_exact)

def print_source_lateral(args):
    "Print the code computing the extra term at the right boundary"

    f_lat, _ = source_lateral()

    print_header("blatter_xz_vanderveen_source_lateral", args)

    h0 = sp.var("thickness")

    print_var(C0(), C)
    print_var(H_exact, h0)

    print_footer(f_lat.subs(H, h0), 0.0)

def print_source_surface(args):
    "Print the code computing the extra term at the top surface"

    f_top, _ = surface_bc()

    print_header("blatter_xz_vanderveen_source_surface", args)

    h0 = sp.var("thickness")

    print_var(C0(), C)
    print_var(H_exact, h0)

    print_footer(f_top.subs(H, h0), 0.0)

def print_basal_beta(args):
    "Print the code computing basal sliding coefficient"

    beta = basal_beta()

    print_header("blatter_xz_vanderveen_beta", args, return_type="double")

    h0 = sp.var("thickness")

    print_var(C0(), C)
    print_var(H_exact, h0)

    print_footer(beta.subs(H, h0))


def print_thickness(args):
    "Print the code computing the ice thickness"

    print_header("blatter_xz_vanderveen_thickness", args, return_type="double")

    print_var(C0(), C)

    print_footer(H_exact)
