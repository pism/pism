#!/usr/bin/env python3
import sympy as sp

func_template = "Vector2 {name}({arguments})"

return_template = """  return {{
    {},
    {}
  }};"""

def join(args):
    return ", ".join(["double " + x for x in args])

def print_var(var, name):
    print("  double " + sp.ccode(var, assign_to=name, standard="c99"))

def print_header(name, args):
    print("")
    print((func_template + " {{").format(name=name, arguments=join(args)))

def declare(name, args):
    print("")
    print((func_template + ";").format(name=name, arguments=join(args)))

def define(f_u, f_v, name, args):
    print("")
    print(func_template.format(name=name, arguments=join(args)))
    print("{")

    tmps, (u, v) = sp.cse([f_u, f_v])

    for variable, value in tmps:
        print_var(value, variable)

    print(return_template.format(sp.ccode(u, standard="c99"),
                                 sp.ccode(v, standard="c99")))

    print("}")
