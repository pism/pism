#!/usr/bin/env python3

"""This script re-implements key features of the `pismr` executable and shows how to write
a custom ocean model component in Python.

Note the

try:
    # some code
except Exception:
    traceback.print_exc()
    raise

blocks below. They are needed to get decent error messages when something in
`PythonOceanModel` fails.

"""

import PISM
import traceback

class PythonOceanModel(PISM.PyOceanModel):
    def __init__(self, grid, standalone=True):
        super().__init__()

        # 2D arrays provided by the base class and used to provide ocean fields to the
        # rest of PISM
        self.shelf_base_temperature = None
        self.shelf_base_mass_flux = None
        self.water_column_pressure = None

        # When used as a part of a PISM run the storage fields (above) are allocated by
        # PISM. The "standalone" mode makes it possible to test this class on its own
        # (e.g. in a Jupyter notebook).
        if standalone:
            self.allocate(grid)

        self.time = grid.ctx().time()
        self.log = grid.ctx().log()
        self.config = grid.ctx().config()

    def init(self, geometry):
        "Initialize the ocean model. This is where you'd read data from an input file, etc."
        try:
            melt_rate = self.config.get_number("ocean.constant.melt_rate", "m year-1")

            self.log.message(2, "* Initializing the Python constant ocean model...\n");
            self.log.message(2, f"  Sub-shelf melt rate set to {melt_rate} m/year.\n");
        except Exception:
            traceback.print_exc()
            raise

    def melting_point_temperature(self, depth, result):
        "Compute melting point temperature [K] of ice at the depth `depth` [meters]."

        T0          = self.config.get_number("constants.fresh_water.melting_point_temperature")
        beta_CC     = self.config.get_number("constants.ice.beta_Clausius_Clapeyron")
        g           = self.config.get_number("constants.standard_gravity")
        ice_density = self.config.get_number("constants.ice.density")

        with PISM.vec.Access([depth, result]):
            for (i, j) in result.grid().points():
                pressure = ice_density * g * depth[i, j]
                result[i, j] = T0 - beta_CC * pressure;

    def update(self, geometry, t, dt):
        """Perform a time step from `t` to `t + dt` [seconds]. Needs to update its state and
        data provided to PISM.

        """
        try:
            start_date = self.time.date(t)
            end_date = self.time.date(t + dt)
            self.log.message(2, f"Python ocean model: time step ({start_date}, {end_date})\n")

            melt_rate     = self.config.get_number("ocean.constant.melt_rate", "m second-1")
            ice_density   = self.config.get_number("constants.ice.density")
            water_density = self.config.get_number("constants.sea_water.density")
            g             = self.config.get_number("constants.standard_gravity")
            mass_flux     = melt_rate * ice_density

            self.shelf_base_mass_flux.set(mass_flux)

            PISM.compute_average_water_column_pressure(geometry, ice_density, water_density, g,
                                                       self.water_column_pressure);

            self.melting_point_temperature(geometry.ice_thickness, self.shelf_base_temperature)
        except Exception:
            traceback.print_exc()
            raise

    def max_timestep(self, t):
        "Return the maximum time step allowed at time `t` [seconds]"
        try:
            # no time step restriction
            return PISM.MaxTimestep("ocean constant")
        except Exception:
            traceback.print_exc()
            raise

    def define_model_state(self, output):
        "Define model state variables in the file `output`."
        try:
            # This model does not have a state but this code shows how to define the state
            # in models that do.
            self.shelf_base_temperature.define(output, PISM.PISM_DOUBLE)
            self.shelf_base_mass_flux.define(output, PISM.PISM_DOUBLE)
        except Exception:
            traceback.print_exc()
            raise

    def write_model_state(self, output):
        "Write model state variables to the file `output`."
        try:
            # This model does not have a state but this code shows how to save the state
            # in models that do.
            self.shelf_base_temperature.write(output)
            self.shelf_base_mass_flux.write(output)
        except Exception:
            traceback.print_exc()
            raise

def main():
    PISM.set_abort_on_sigint(True)
    context = PISM.Context().ctx

    usage = \
      """  pism.py -i IN.nc [-bootstrap] [OTHER PISM & PETSc OPTIONS]\n"
      where:
        -i                   IN.nc is input file in NetCDF format: contains PISM-written model state
        -bootstrap           enable heuristics to produce an initial state from an incomplete input
      notes:
        * option -i is required
      """

    if PISM.show_usage_check_req_opts(context.log(), "PISM (basic evolution run mode)" ,
                                      ["-i"], usage):
        return

    grid = PISM.Grid.FromOptions(context)

    model = PISM.IceModel(grid, context)

    # Allocate the ocean model. The `ocean` below is stored in a variable to ensure that
    # it does not get de-allocated before `model` is done.
    #
    # This
    #
    # model.set_python_ocean_model(PythonOceanModel(grid, standalone=False))
    #
    # would not work.

    ocean = PythonOceanModel(grid, standalone=False)

    # Tell `model` to use the ocean model implemented here. This overrides the command
    # line option `-ocean` (if present).
    model.set_python_ocean_model(ocean)

    model.init()

    model.run()

    model.save_results()

def standalone_test():
    context = PISM.Context().ctx

    Lx = 10e4
    Mx = 100
    grid = PISM.Grid.Shallow(context, Lx, Lx, 0, 0, Mx, Mx, PISM.CELL_CENTER, PISM.NOT_PERIODIC)
    grid.report_parameters()

    geometry = PISM.Geometry(grid)

    ocean = PythonOceanModel(grid, standalone=True)

    with PISM.vec.Access(geometry.ice_thickness):
        for (i, j) in grid.points():
            r = PISM.radius(grid, i, j)
            H = max(1000 - (r / 2000)**2, 0)
            geometry.ice_thickness[i, j] = H

    geometry.bed_elevation.set(0.0)
    geometry.ensure_consistency(0.0)

    ocean.init(geometry)

    ocean.update(geometry, 0, 86400)

    ocean.shelf_base_temperature.dump("shelf_base_temperature.nc")

if __name__ == "__main__":
    main()
