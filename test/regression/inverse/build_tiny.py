#! /usr/bin/env python
#

import sys
import petsc4py
petsc4py.init(sys.argv)
from petsc4py import PETSc

import PISM
import math

# Default constants that  may get overridden later.

Ly = 25e3  # 25 km
Lx = 50e3  # 50 km
Lz = 4000
My = 13
Mx = 23
Mz = 21

sea_level = 0      # m sea level elevation
H0 = 60.           # ice thickness at cliff
alpha = 0.008       # constant surface slope
Lext = 15e3          # width of strip beyond cliff
Lstream_x = 50e3
Lstream_y = 30e3

Hext = 0.    # m ice thickeness beyond the cliff

tauc_hi = 2e6       # Pa
tauc_lo = 1e4       # Pa
tauc_free_bedrock = 0  # Will get set later

EC = PISM.EnthalpyConverter(PISM.Context().config)
enth0 = EC.enthalpy(273.15, 0.01, 0)  # 0.01 water fraction
bed0 = 0


def geometry(x, y):
    x0 = -Lx + Lext
    if x < x0:
        return (0, Hext)
    return (0, H0 + alpha * (x - x0))


def stream_tauc(x, y):
    x0 = -Lx + Lext
    if x < x0:
        return tauc_free_bedrock
    if x < x0 + Lstream_x:
        if abs(y) < Lstream_y / 2:
            return tauc_lo
    return tauc_hi


# The main code for a run follows:
if __name__ == '__main__':
    PISM.set_abort_on_sigint(True)
    context = PISM.Context()

    Mx = PISM.optionsInt("-Mx", "Number of grid points in x-direction", default=Mx)
    My = PISM.optionsInt("-My", "Number of grid points in y-direction", default=My)
    output_filename = PISM.optionsString("-o", "output file", default="tiny.nc")
    verbosity = PISM.optionsInt("-verbose", "verbosity level", default=3)

    # Build the grid.
    config = PISM.Context().config

    p = PISM.GridParameters(config)
    p.Mx = Mx
    p.My = My
    p.Lx = Lx
    p.Ly = Ly
    z = PISM.IceGrid.compute_vertical_levels(Lz, Mz, PISM.EQUAL, 4.0)
    p.z = PISM.DoubleVector(z)
    p.ownership_ranges_from_options(context.size)
    p.periodicity = PISM.NOT_PERIODIC
    grid = PISM.IceGrid(context.ctx, p)

    vecs = PISM.model.ModelVecs(grid.variables())
    vecs.add(PISM.model.createIceSurfaceVec(grid))
    vecs.add(PISM.model.createIceThicknessVec(grid))
    vecs.add(PISM.model.createBedrockElevationVec(grid))
    vecs.add(PISM.model.createYieldStressVec(grid), 'tauc')
    vecs.add(PISM.model.createEnthalpyVec(grid), 'enthalpy')
    vecs.add(PISM.model.createIceMaskVec(grid))
    vecs.add(PISM.model.createNoModelMaskVec(grid), 'no_model_mask')
    vecs.add(PISM.model.create2dVelocityVec(grid,  name='_ssa_bc', desc='SSA Dirichlet BC'))

    # Set constant coefficients.
    vecs.enthalpy.set(enth0)

    # Build the continent
    bed = vecs.bedrock_altitude
    thickness = vecs.land_ice_thickness

    with PISM.vec.Access(comm=[bed, thickness]):
        for (i, j) in grid.points():
            x = grid.x(i)
            y = grid.y(j)
            (b, t) = geometry(x, y)
            bed[i, j] = b
            thickness[i, j] = t

    # Compute mask and surface elevation from geometry variables.
    gc = PISM.GeometryCalculator(grid.ctx().config())
    gc.compute(sea_level, bed, thickness, vecs.mask, vecs.surface_altitude)

    tauc = vecs.tauc
    mask = vecs.mask
    tauc_free_bedrock = config.get_double('high_tauc')
    with PISM.vec.Access(comm=tauc, nocomm=mask):
        for (i, j) in grid.points():
            tauc[i, j] = stream_tauc(grid.x(i), grid.y(j))

    vecs.vel_ssa_bc.set(0.0)
    no_model_mask = vecs.no_model_mask
    no_model_mask.set(0)
    with PISM.vec.Access(comm=[no_model_mask]):
        for (i, j) in grid.points():
            if (i == 0) or (i == grid.Mx() - 1) or (j == 0) or (j == grid.My() - 1):
                no_model_mask[i, j] = 1

    pio = PISM.PIO(grid.com, "netcdf3")
    pio.open(output_filename, PISM.PISM_READWRITE_MOVE)
    PISM.define_time(pio,
                     grid.ctx().config().get_string("time_dimension_name"),
                     "365_day",
                     "seconds since 1-1-1",
                     grid.ctx().unit_system())
    PISM.append_time(pio,
                     grid.ctx().config().get_string("time_dimension_name"),
                     0.0)
    pio.close()
    vecs.writeall(output_filename)
    PISM.util.writeProvenance(output_filename)
