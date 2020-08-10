
"""This script documents an alternative way of writing the Blatter-Pattyn system which
makes it easier to see that it can be thought of as a "convection-diffusion" equation with
the diffusivity and "velocity" depending on viscosity.

This will help me implement the streamline-upwind Petrov-Galerkin method that should
eliminate spurious oscillations. (I hope so, anyway.)
"""

from blatter import *

u = sp.Function("u")(x, y, z)
v = sp.Function("v")(x, y, z)

# Re-define ice viscosity as an "undefined function" to keep sympy from expanding it.
eta   = sp.Function("eta")(x, y, z)
eta_x = eta.diff(x)
eta_y = eta.diff(y)
eta_z = eta.diff(z)

s = sp.Function("s")(x, y)
sp.var("rho, g", real=True)

# "velocity" vectors for the x-component
A = [sp.Matrix([4 * eta_x, eta_y, eta_z]),
     sp.Matrix([2 * eta_y, eta_x, 0])]

# "velocity" vectors for the y-component
B = [sp.Matrix([eta_y, 2 * eta_x, 0]),
     sp.Matrix([eta_x, 4 * eta_y, eta_z])]

def diffusive_part(n):
    return -eta * div(M(u, v).row(n))

def convective_part(n):
    return -(A[n].T * grad(u) + B[n].T * grad(v))[0]

def body_force(n):
    return rho * g * grad(s)[n]

def eq_convection(n):
    "The BP system in the convection-diffusion form"
    return convective_part(n) + diffusive_part(n) + body_force(n)

def eq_standard(n):
    "The BP system in the 'standard' form"
    return - div(eta * M(u, v).row(n)) + body_force(n)

def are_equivalent(n):
    "Return true if the two forms of the BP system are equivalent."
    return (eq_convection(n) - eq_standard(n)).expand() == 0
