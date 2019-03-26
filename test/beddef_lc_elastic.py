#!/usr/bin/env python

import numpy as np
from scipy.integrate import dblquad

import PISM
from PISM.util import convert


ge = PISM.greens_elastic()

def ge_integrand(eta, xi, dx, dy, p, q):
    xi_shift  = p * dx - xi
    eta_shift = q * dy - eta
    r         = np.sqrt(xi_shift * xi_shift + eta_shift * eta_shift)

    return ge(r)

def f(dx, dy, p, q):
    return dblquad(ge_integrand,
                   -dx/2.0, dx/2.0,
                   lambda x: -dy / 2.0, lambda x: dy / 2.0,
                   (dx, dy, p, q), epsrel=1e-8)[0]

def lrm(Mx, My, dx, dy):
    "Compute the load response matrix, taking advantage of its symmetry."
    Mx2 = int(Mx) // 2
    My2 = int(My) // 2

    a = np.zeros((My, Mx))

    # top half
    for j in range(My2 + 1):
        # top left quarter
        for i in range(Mx2 + 1):
            p = Mx2 - i
            q = My2 - j

            a[j, i] = f(dx, dy, p, q)

        # top right quarter
        for i in range(Mx2 + 1, Mx):
            a[j, i] = a[j, 2 * Mx2 - i]

    # bottom half
    for j in range(My2 + 1, My):
        for i in range(Mx):
            a[j, i] = a[2 * My2 - j, i]

    return a

def lrm_naive(Mx, My, dx, dy):
    "Naive LRM computation used to test the optimized one"
    Mx2 = int(Mx) // 2
    My2 = int(My) // 2

    a = np.zeros((My, Mx))

    # top half
    for j in range(My):
        for i in range(Mx):
            p = abs(Mx2 - i)
            q = abs(My2 - j)

            a[j, i] = f(dx, dy, p, q)

    return a

Lx = convert(2000, "km", "m")

Mx = 101
My = Mx

def test_elastic(grid):
    bed_model = PISM.LingleClark(grid)

    ice_thickness = PISM.IceModelVec2S(grid, "thk", PISM.WITHOUT_GHOSTS)

    bed = PISM.IceModelVec2S(grid, "topg", PISM.WITHOUT_GHOSTS)

    bed_uplift = PISM.IceModelVec2S(grid, "uplift", PISM.WITHOUT_GHOSTS)

    sea_level = PISM.IceModelVec2S(grid, "sea_level", PISM.WITHOUT_GHOSTS)

    # start with a flat bed, no ice, and no uplift
    bed.set(0.0)
    bed_uplift.set(0.0)
    ice_thickness.set(0.0)
    sea_level.set(-1000.0)

    bed_model.bootstrap(bed, bed_uplift, ice_thickness, sea_level)

    Mx2 = int(grid.Mx()) // 2
    My2 = int(grid.My()) // 2

    # add the disc load
    with PISM.vec.Access(nocomm=ice_thickness):
        for (i, j) in grid.points():
            if i == Mx2 and j == My2:
            # if abs(i - Mx2) < 4 and abs(j - My2) < 4:
                ice_thickness[i, j] = 1000.0

    # dt of zero disables the viscous part of the model, so all we get is the elastic
    # response
    bed_model.step(ice_thickness, sea_level, 0)

    return ice_thickness.numpy(), bed_model.bed_elevation().numpy(), bed_model.total_displacement().numpy()

def elastic_model(thk, LRM, rho):
    return np.fft.ifft2(np.fft.fft2(rho * thk) * np.fft.fft2(LRM))


if __name__ == "__main__":
    ctx = PISM.Context()
    ctx.config.set_boolean("bed_deformation.lc.elastic_model", True)
    ctx.config.set_double("bed_deformation.lc.grid_size_factor", 4)

    grid = PISM.IceGrid.Shallow(ctx.ctx, Lx, Lx, 0, 0, Mx, My, PISM.CELL_CORNER, PISM.NOT_PERIODIC)

    dx = grid.dx()
    dy = grid.dy()

    # load_response_matrix = lrm(Mx, My, dx, dy)

    H, b, db = test_elastic(grid)

    import pylab as plt
    size=(5,5)
    plt.figure(figsize=size)
    plt.imshow(H)
    plt.title("thickness")
    plt.colorbar()

    plt.figure(figsize=size)
    plt.imshow(b)
    plt.title("bed")
    plt.colorbar()

    plt.figure(figsize=size)
    plt.imshow(db)
    plt.title("displacement")
    plt.colorbar()

    plt.show()
