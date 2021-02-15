#!/usr/bin/env python3
import sympy as sp

func_template = "Vector2 {name}({arguments})"

def declare(name, args):
    arguments = ", ".join(["double " + x for x in args])
    print("")
    print((func_template + ";").format(arguments=arguments, name=name))

def define(f_u, f_v, name, args):
    arguments = ", ".join(["double " + x for x in args])

    print("")
    print(func_template.format(arguments=arguments, name=name))

    print("{")

    tmps, (u, v) = sp.cse([f_u, f_v])

    for variable, value in tmps:
        print("  double " + sp.ccode(value, assign_to=variable))

    print("  return {")
    print("    {},".format(sp.ccode(u, standard="c99")))
    print("    {}".format(sp.ccode(v, standard="c99")))
    print("  };")

    print("}")
