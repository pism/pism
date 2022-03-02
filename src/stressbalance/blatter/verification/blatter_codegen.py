#!/usr/bin/env python3
import sympy as sp

func_template = "{return_type} {name}({arguments})"

return_template = """  return {{
    {},
    {}
  }};"""

def code(a, **kwargs):
    return sp.ccode(a, standard="c99", **kwargs)

def join(args):
    return ", ".join(["double " + x for x in args])

def print_var(var, name):
    print("  double " + code(var, assign_to=name))

def print_header(name, args, return_type="Vector2d"):
    print("")
    print((func_template + " {{").format(return_type=return_type,
                                         name=name,
                                         arguments=join(args)))

def print_footer(a, b=None):
    if b is not None:
        print(return_template.format(code(a), code(b)))
    else:
        print("  return {};".format(code(a)))
    print("}")

def declare(name, args, return_type="Vector2d"):
    print("")
    print((func_template + ";").format(return_type=return_type,
                                       name=name,
                                       arguments=join(args)))

def define(f_u, f_v, name, args):
    print("")
    print(func_template.format(return_type="Vector2d",
                               name=name,
                               arguments=join(args)))
    print("{")

    tmps, (u, v) = sp.cse([f_u, f_v])

    for variable, value in tmps:
        print_var(value, variable)

    print(return_template.format(code(u),
                                 code(v)))

    print("}")
