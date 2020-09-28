import sympy as sp
from sympy.core import S

# coordinate variables
sp.var("x, y, z", real=True)

# velocity components
u = sp.Function("u")(x, y, z)
v = sp.Function("v")(x, y, z)

# ice viscosity, defined this way to use compact expressions for partial derivatives
eta_u, eta_v = [f(u, v) for f in sp.symbols("eta_(u:v)", cls=sp.Function)]
class Eta(sp.Function):
    def fdiff(self, argindex):
        # partial derivatives of the ice viscosity eta
        return [eta_u, eta_v][argindex - 1]

def div(f):
    "Divergence"
    return f[0].diff(x) + f[1].diff(y) + f[2].diff(z)

def grad(f):
    "Gradient"
    return sp.Matrix([f.diff(x), f.diff(y), f.diff(z)])

# ice hardness
B = sp.var("B", positive=True)

s = sp.Function("s")(x, y)

def second_invariant(u, v):
    """Second invariant of the strain rate tensor"""
    ux = u.diff(x)
    uy = u.diff(y)
    uz = u.diff(z)

    vx = v.diff(x)
    vy = v.diff(y)
    vz = v.diff(z)

    # In the BP approximation we assume that dw/dx << du/dz and dw/dy << dv/dz, so these
    # terms are omitted.
    wx = S(0)
    wy = S(0)

    return ux**2 + vy**2 + ux * vy + S(1) / 4 * ((uy + vx)**2 + (uz + wx)**2 + (vz + wy)**2)

def eta(u, v, n):
    "Ice viscosity"

    n = S(n)

    gamma = second_invariant(u, v)

    return S(1) / 2 * B * gamma**((1 - n) / (2 * n))

def M(u, v):
    "'Effective' strain rate tensor corresponding to the Blatter-Pattyn system"
    return sp.Matrix([[4 * u.diff(x) + 2 * v.diff(y), u.diff(y) + v.diff(x), u.diff(z)],
                      [u.diff(y) + v.diff(x), 2 * u.diff(x) + 4 * v.diff(y), v.diff(z)]])

def source_term(E, u, v):
    "Compute the source term required by a given 'exact' solution."
    # we assume that the system is written in this form:
    #
    # -div(eta*M) + f = 0
    #
    # where eta is ice viscosity and M is the effective strain rate tensor
    #
    # in the "regular" context
    #
    # f = rho * g * grad(s)

    f_u = div(E * M(u, v).row(0))
    f_v = div(E * M(u, v).row(1))

    return f_u, f_v
