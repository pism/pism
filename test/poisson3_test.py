#!/usr/bin/env python

"""This script is used to test Poisson3, the 3D Q1 FEM solver for the Poisson equation.

"""

import PISM
import PISM.testing as pt
import numpy as np

config = PISM.Context().config

def run(M):

    grid = pt.shallow_grid(Mx=M, My=M, Lx=1, Ly=1)

    pp = PISM.Poisson3(grid, M)

    inp = PISM.StressBalanceInputs()

    pp.update(inp, False)

    return pp

def write(model):

    output = PISM.util.prepare_output(config.get_string("output.file_name"))

    model.exact().write(output)
    model.solution().write(output)

    output.close()

def verify():

    N = 4

    dxs = []
    errors = []
    for k in range(N):
        M = 5 * 2**k

        P = run(M)

        dxs += [2.0 / M]
        errors += [P.error()]

    f = np.polyfit(np.log(dxs), np.log(errors), 1)
    F = np.exp(np.polyval(f, np.log(dxs)))

    write(P)

    if PISM.Context().rank == 0:

        import pylab as plt
        plt.loglog(dxs, errors, ".-")
        plt.loglog(dxs, F, "--", label=f"dx^{f[0]:2.3f}")
        plt.grid()
        plt.legend()
        plt.xlabel("dx")
        plt.ylabel("error")
        plt.show()

if __name__ == "__main__":
    verify()
