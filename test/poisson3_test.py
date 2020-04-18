#!/usr/bin/env python

"""This script is used to test Poisson3, the 3D Q1 FEM solver for the Poisson equation.

"""

import PISM
import PISM.testing as pt
import numpy as np
import pylab as plt

config = PISM.Context().config

n_levels = 2

N = 6

L = 1.0

def run(Mx, Mz, L):

    grid = pt.shallow_grid(Mx=Mx, My=Mx, Lx=L, Ly=L)

    pp = PISM.Poisson3(grid, Mz, n_levels)

    inp = PISM.StressBalanceInputs()

    pp.update(inp, False)

    return pp

def write(model, filename):

    output = PISM.util.prepare_output(filename)

    model.exact().write(output)
    model.solution().write(output)

    output.close()

def verify_1():
    "Verification test checking the code excluding ice-free elements"

    dxs = []
    errors = []
    for k in range(1, N):
        M = 2 * 2**k + 1

        Mx = M + 2
        Mz = M
        dx = 2.0 * L / (M - 1)

        P = run(Mx, Mz, L + dx)

        dxs += [dx]
        errors += [P.error()]

    f = np.polyfit(np.log(dxs), np.log(errors), 1)
    F = np.exp(np.polyval(f, np.log(dxs)))

    write(P, "output-ice-free.nc")

    if PISM.Context().rank == 0:

        plt.figure()
        plt.title("ice-free")
        plt.loglog(dxs, errors, ".-", label="errors")
        plt.loglog(dxs, F, "--", label=f"polyfit: dx^{f[0]:2.3f}")
        plt.legend()
        plt.xlabel("dx")
        plt.ylabel("error")

def verify_2():
    "Verification test with the whole domain filled with ice"

    dxs = []
    errors = []
    for k in range(1, N):
        M = 2 * 2**k + 1

        dx = 2.0 * L / (M - 1)

        P = run(M, M, L)

        dxs += [dx]
        errors += [P.error()]

    f = np.polyfit(np.log(dxs), np.log(errors), 1)
    F = np.exp(np.polyval(f, np.log(dxs)))

    write(P, "output-icy.nc")

    if PISM.Context().rank == 0:

        plt.figure()
        plt.title("icy")
        plt.loglog(dxs, errors, ".-", label="errors")
        plt.loglog(dxs, F, "--", label=f"polyfit: dx^{f[0]:2.3f}")
        plt.legend()
        plt.xlabel("dx")
        plt.ylabel("error")

if __name__ == "__main__":
    verify_1()
    verify_2()
    plt.show()
