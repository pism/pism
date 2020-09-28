from sympy import exp, sin, cos, pi

from blatter import x, y, z, source_term, eta
from manufactured_solutions import define, declare

def exact():
    """X-Y verification test

    Constant bed elevation and ice thickness, constant hardness.

    """

    u0 = exp(x) * sin(2 * pi * y)
    v0 = exp(x) * cos(2 * pi * y)

    return u0, v0

def print_code(header=False):
    if header:
        declare(name="blatter_xy_exact", args=["x", "y"])
        declare(name="blatter_xy_source", args=["x", "y", "B"])
        return

    u0, v0 = exact()
    f_u, f_v = source_term(eta(u0, v0, 3), u0, v0)

    def cleanup(expr):
        return expr.factor().collect([sin(2*pi*y), cos(2*pi*y)])

    define(u0, v0, name="blatter_xy_exact", args=["x", "y"])
    define(cleanup(f_u), cleanup(f_v), name="blatter_xy_source", args=["x", "y", "B"])
