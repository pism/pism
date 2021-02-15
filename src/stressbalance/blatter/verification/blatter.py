import sympy as sp
from sympy.core import S

# References:
#
# Lipscomb2019: W. H. Lipscomb et al., “Description and evaluation of the Community Ice
# Sheet Model (CISM) v2.1,” Geoscientific Model Development, vol. 12, no. 1, Art. no. 1,
# Jan. 2019, doi: 10.5194/gmd-12-387-2019.
#
# Greve2009: R. Greve and H. Blatter, Dynamics of Ice Sheets and Glaciers. Springer Berlin
# Heidelberg, 2009.
#

# coordinate variables
sp.var("x, y, z", real=True)

# velocity components
u = sp.Function("u")(x, y, z)
v = sp.Function("v")(x, y, z)

def div(f):
    "Divergence"
    return f[0].diff(x) + f[1].diff(y) + f[2].diff(z)

def grad(f):
    "Gradient"
    return sp.Matrix([f.diff(x), f.diff(y), f.diff(z)])

# ice hardness
B = sp.var("B", positive=True)

# ice surface elevation
s = sp.Function("s")(x, y)

def edot(x1, x2):
    """Elements of the strain rate tensor"""
    # express w_z using incompressibility:
    if x1 == z and x2 == z:
        return -(u.diff(x) + v.diff(y))

    # In the BP approximation partial derivatives w_x and w_y are omitted.
    #
    # See Greve2009, section 5.3
    V = {x : u, y : v, z : S(0)}

    return (S(1) / 2 * (V[x1].diff(x2) + V[x2].diff(x1))).factor()

def second_invariant(U, V):
    "Second invariant of the strain rate tensor"

    # strain rate tensor
    D = sp.Matrix([[edot(x,x), edot(x,y), edot(x,z)],
                   [edot(y,x), edot(y,y), edot(y,z)],
                   [edot(z,x), edot(z,y), edot(z,z)]])

    # second invariant (Greve2009, equation 2.42)
    #
    # Note: D.trace() is zero.
    II = S(1) / 2 * ((D**2).trace() - D.trace()**2)

    return II.subs({u: U, v: V}).doit()

def eta(u, v, n):
    "Ice viscosity"

    n = S(n)

    gamma = second_invariant(u, v)

    # Greve2009, equation 4.22
    return S(1) / 2 * B * gamma**((1 - n) / (2 * n))

def M(U, V):
    "'Effective' strain rate tensor corresponding to the Blatter-Pattyn approximation"

    # See equation 8 in Lipscomb2019, compare to equation 5.70 in Greve2009
    M = sp.Matrix([[2 * edot(x, x) + edot(y, y), edot(x, y), edot(x, z)],
                   [edot(x, y), edot(x, x) + 2 * edot(y, y), edot(y, z)]])

    return M.subs({u: U, v: V}).doit()

def source_term(E, u, v):
    "Compute the source term required by a given 'exact' solution."
    # we assume that the system is written in this form:
    #
    # -div(2 * eta * M) + f = 0
    #
    # where eta is ice viscosity and M is the effective strain rate tensor
    #
    # in the "regular" context
    #
    # f = rho * g * grad(s)

    f_u = div(2 * E * M(u, v).row(0))
    f_v = div(2 * E * M(u, v).row(1))

    return f_u, f_v
