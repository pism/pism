#!/usr/bin/env python

"Test the pointwise isostasy model."

import numpy
import PISM

ctx  = PISM.Context().ctx

def exact(ice_thickness_change):
    "Exact deflection corresponding to a given thickness change."

    config = PISM.Context().config

    ice_density    = config.get_double("constants.ice.density")
    mantle_density = config.get_double("bed_deformation.mantle_density")

    return -(ice_density / mantle_density) * ice_thickness_change

def bed_def_iso(ice_thickness_change):
    "Use the pointwise isostasy model to compute plate deflection."

    # grid size and domain size are irrelevant
    M = 3
    L = 1000e3

    grid = PISM.IceGrid.Shallow(ctx, L, L, 0, 0, M, M, PISM.CELL_CORNER, PISM.NOT_PERIODIC)

    geometry = PISM.Geometry(grid)
    geometry.bed_elevation.set(0.0)
    geometry.sea_level_elevation.set(0.0)
    geometry.ice_thickness.set(0.0)
    geometry.ice_area_specific_volume.set(0.0)
    geometry.ensure_consistency(0.0)

    # uplift is required (but not used)
    bed_uplift = PISM.IceModelVec2S(grid, "uplift", PISM.WITHOUT_GHOSTS)
    bed_uplift.set(0.0)

    bed_model = PISM.PointwiseIsostasy(grid)
    bed_model.bootstrap(geometry.bed_elevation,
                        bed_uplift,
                        geometry.ice_thickness,
                        geometry.sea_level_elevation)

    geometry.ice_thickness.set(ice_thickness_change)

    # time step duration is irrelevant
    bed_model.update(geometry.ice_thickness, geometry.sea_level_elevation, 0, 1)

    return bed_model.bed_elevation().numpy()[0,0]

def beddef_iso_test():
    "Test the pointwise isostasy model"

    dH = 1000.0

    numpy.testing.assert_almost_equal(bed_def_iso(dH), exact(dH))
