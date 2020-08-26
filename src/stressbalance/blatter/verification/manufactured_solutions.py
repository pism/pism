#!/usr/bin/env python3

import sympy as sp
from sympy import exp, sin, cos, pi
from sympy.core import S

from blatter import *

nx, ny, nz = sp.var("n_(x:z)")

N = sp.Matrix([nx, ny, nz])

def print_source(f_u, f_v, name="func", args=[], header=False):
    arguments = ", ".join(["double " + x for x in args])

    print("")
    print("Vector2 {name}({arguments}){semicolon}".format(arguments=arguments,
                                                          name=name,
                                                          semicolon=";" if header else ""))
    if header:
        return

    print("{")

    tmps, (u, v) = sp.cse([f_u, f_v])

    for variable, value in tmps:
        print("  double " + sp.ccode(value, assign_to=variable))

    print("  return {")
    print("    {},".format(sp.ccode(u, standard="c99")))
    print("    {}".format(sp.ccode(v, standard="c99")))
    print("  };")

    print("}")

def print_exact(u, v, name="exact", args=[], header=False):
    arguments = ", ".join(["double " + x for x in args])
    print("")
    print("Vector2 {name}({arguments}){semicolon}".format(arguments=arguments, name=name,
                                                          semicolon=";" if header else ""))
    if header:
        return
    print("{")
    print("  return {")
    print("    {},".format(sp.ccode(u, standard="c99")))
    print("    {}".format(sp.ccode(v, standard="c99")))
    print("  };")
    print("}")

def exact_xy():
    """X-Y verification test

    Constant bed elevation and ice thickness, constant hardness, periodic boundary
    conditions in X and Y directions.
    """

    u0 = exp(x) * sin(2 * pi * y)
    v0 = exp(x) * cos(2 * pi * y)

    return u0, v0

def source_xy_albany():

    half   = S(1) / 2
    ex = exp(x)
    sin2piy = sin(2 * pi * y)
    cos2piy = cos(2 * pi * y)
    pi2 = 2 * pi

    mut = ex * ((1 + pi2**2 - pi2) * sin2piy**2 + S(1)/4 * (pi2 + 1)**2 * cos2piy**2)**half
    mut_x = mut
    mut_y = 3 * half * ((pi * (1 + pi2**2 - 4*pi) * cos2piy * sin2piy * ex) /
                        ((1 + pi2**2 - pi2) * sin2piy**2 + S(1)/4 * (pi2 + 1)**2 * cos2piy**2))

    e_xx = ex * sin2piy
    e_yy = -pi2 * ex * sin2piy
    e_xy = half * (pi2 + 1) * ex * cos2piy

    mu = half * A**(-1/n) * mut**(1/n - 1)

    f_u = (2 * mu * ex * sin2piy * (2 - 3*pi - 2*pi**2) +
           A**(-1/n) * (1/n - 1) * mut**(1/n - 2) * (mut_x * (2 * e_xx + e_yy) + mut_y * e_xy))
    f_v = (2 * mu * ex * cos2piy * (3*pi + half - 8 * pi**2) +
           A**(-1/n) * (1/n - 1) + mut**(1/n - 2) * (mut_x * e_xy + mut_y * (e_xx + 2 * e_xy)))

    return f_u, f_v

def print_xy(header=False):
    u0, v0 = exact_xy()
    f_u, f_v = source_term(eta(u0, v0), u0, v0)

    def cleanup(expr):
        return expr.factor().collect([sin(2*pi*y), cos(2*pi*y)])

    print_exact(u0, v0, name="exact_xy", args=["x", "y"], header=header)
    print_source(cleanup(f_u), cleanup(f_v), name="source_xy", args=["x", "y", "B"], header=header)

def print_xz(header=False):

    u0, v0 = exact_xz()
    f_u, f_v = source_term(eta(u0, v0), u0, v0)

    args = ["x", "z", "B", "rhog", "s0", "alpha", "H"]
    print_exact(u0, v0, name="exact_xz", args=args, header=header)
    print_source(f_u, f_v, name="source_xz", args=args, header=header)

def main(header=False):
    if header:
        print("#include <cmath>")
        print("")
        print('#include "pism/util/Vector2.hh"')
    else:
        print('#include "manufactured_solutions.hh"')

    print("")
    print("namespace pism {")

    print_xy(header)
    print_xz(header)

    print("")
    print("} // end of namespace pism")

if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument("--header", dest="header", action="store_true",
                        help="print function declarations for the header file")

    options = parser.parse_args()

    main(options.header)
