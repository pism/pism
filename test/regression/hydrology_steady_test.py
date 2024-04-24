#!/usr/bin/env python3
from unittest import TestCase
import numpy as np
import os

import PISM
ctx = PISM.Context()
ctx.log.set_threshold(1)

class SteadyHydrology(TestCase):
    def setUp(self):
        # domain size
        C = 1000.0
        x0 = 0.0
        y0 = 0.0
        Lx = C*0.5
        Ly = C*0.5
        # grid size
        Mx = 101
        My = 101

        grid = PISM.Grid.Shallow(ctx.ctx, Lx, Ly, x0, y0, Mx, My,
                                    PISM.CELL_CENTER, PISM.NOT_PERIODIC)
        self.grid = grid

        geometry = PISM.Geometry(grid)
        self.geometry = geometry

        with PISM.vec.Access(nocomm=[geometry.bed_elevation, geometry.ice_thickness]):
            for (i, j) in grid.points():
                x = grid.x(i)
                y = grid.y(j)
                geometry.bed_elevation[i, j] = (x / C)**2 + y / C + 0.25

                if geometry.bed_elevation[i, j] >= 0:
                    geometry.ice_thickness[i, j] = 1000.0
        geometry.sea_level_elevation.set(0.0)
        geometry.ensure_consistency(0.0)

        surface_input_rate = PISM.Scalar(grid, "water_input_rate")
        surface_input_rate.metadata(0).long_name("water input rate").units("kg m^-2 s^-1")
        self.surface_input_rate = surface_input_rate

        self.surface_input_file = "hydrology_steady_surface_input.nc"
        surface_input_rate.dump(self.surface_input_file)

        # center of the patch of non-zero input
        cx = 0.0
        cy = 0.25 * C
        # size of the patch
        R = 0.125 * C
        with PISM.vec.Access(nocomm=surface_input_rate):
            for (i, j) in grid.points():
                x = grid.x(i)
                y = grid.y(j)
                if abs(x - cx) < R and abs(y - cy) < R:
                    surface_input_rate[i, j] = 1000.0
                else:
                    surface_input_rate[i, j] = 0

        ctx.config.set_flag("hydrology.add_water_input_to_till_storage", False)
        ctx.config.set_string("hydrology.surface_input.file", self.surface_input_file)

        self.model = PISM.SteadyState(grid)

        zero = PISM.Scalar(grid, "zero")
        zero.set(0.0)
        self.zero = zero

        self.model.init(zero, zero, zero)

    def tearDown(self):
        os.remove(self.surface_input_file)

    def divergence_theorem_test(self):
        "Test that the total input equals the total flux through the boundary."
        grid = self.grid

        inputs = PISM.HydrologyInputs()

        inputs.no_model_mask = None
        inputs.geometry = self.geometry
        inputs.basal_melt_rate = self.zero
        inputs.ice_sliding_speed = self.zero
        inputs.surface_input_rate = self.surface_input_rate

        dt = self.model.max_timestep(0).value()
        self.model.update(0, dt, inputs)

        flux_magnitude = PISM.Scalar(grid, "flux_magnitude")

        PISM.compute_magnitude(self.model.flux(), flux_magnitude)

        # Compute the total input. This approximates a double integral, hence
        # the "dx dy" factor.
        total_input = PISM.sum(self.model.surface_input_rate()) * (grid.dx() * grid.dy())

        # Compute the total flux through the grounding line. This is not
        # exactly what we want, but it's close. It would be better to use the
        # flux on the staggered grid as computed internally, but the flux
        # magnitude will do, especially in the dx=dy case. This approximates a
        # line integral over the grounding line, hence the dx factor below.
        total_flux = 0.0
        cell_type = self.geometry.cell_type
        with PISM.vec.Access(nocomm=[cell_type, flux_magnitude]):
            for (i, j) in grid.points():
                if cell_type.ice_free_ocean(i, j) and cell_type.next_to_grounded_ice(i, j):
                   total_flux += flux_magnitude[i, j] * grid.dx()

        total_flux = PISM.GlobalSum(ctx.com, total_flux)

        # This is the relative error. Note that it is not sensitive to the
        # value of hydrology.steady.volume_ratio.
        relative_error = np.fabs(total_input - total_flux) / total_input

        assert relative_error < 1e-5
        ctx.log.message(1, "relative error: {}\n".format(relative_error))

    def write_results(self):
        geometry = self.geometry
        model = self.model

        output_file = ctx.config.get_string("output.file")

        f = PISM.util.prepare_output(output_file)

        geometry.bed_elevation.write(f)
        geometry.cell_type.write(f)
        geometry.ice_thickness.write(f)

        model.surface_input_rate().write(f)
        model.flux().write(f)

        Q = PISM.Scalar(self.grid, "water_flux_magnitude", PISM.WITHOUT_GHOSTS)
        Q.set_attrs("", "magnitude of the water flux", "m2 s-1", "m2 s-1", "", 0)
        Q.set_to_magnitude(model.flux())
        Q.write(f)

        f.close()

if __name__ == "__main__":

    t = SteadyHydrology()
    t.setUp()
    t.divergence_theorem_test()
    t.write_results()
