#! /usr/bin/env python
#
# Copyright (C) 2011, 2012, 2014, 2015, 2016, 2017, 2018 David Maxwell and Constantine Khroulev
#
# This file is part of PISM.
#
# PISM is free software; you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation; either version 3 of the License, or (at your option) any later
# version.
#
# PISM is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License
# along with PISM; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

import PISM
from petsc4py import PETSc
import os

context = PISM.Context()
ctx = context.ctx
config = context.config

PISM.set_abort_on_sigint(True)

usage = """\
sia_forward.py -i IN.nc [-o file.nc]
  where:
    -i      IN.nc is input file in NetCDF format: contains PISM-written model state
  notes:
    * -i is required
"""

PISM.show_usage_check_req_opts(ctx.log(), "sia_forward.py", ["-i"], usage)

input_filename = config.get_string("input.file")
if len(input_filename) == 0:
    import sys
    sys.exit(1)

config.set_string("output.file_name", "sia_" + os.path.basename(input_filename), PISM.CONFIG_DEFAULT)

output_file = config.get_string("output.file_name")
is_regional = PISM.OptionBool("-regional", "Compute SIA using regional model semantics")

registration = PISM.CELL_CENTER
if is_regional:
    registration = PISM.CELL_CORNER

input_file = PISM.PIO(ctx.com(), "netcdf3", input_filename, PISM.PISM_READONLY)
grid = PISM.IceGrid.FromFile(ctx, input_file, "enthalpy", registration)

config.set_boolean("basal_resistance.pseudo_plastic.enabled", False)

enthalpyconverter = PISM.EnthalpyConverter(config)

modeldata = PISM.model.ModelData(grid)
modeldata.setPhysics(enthalpyconverter)

vecs = modeldata.vecs

vecs.add(PISM.model.createIceSurfaceVec(grid))
vecs.add(PISM.model.createIceThicknessVec(grid))
vecs.add(PISM.model.createBedrockElevationVec(grid))
vecs.add(PISM.model.createEnthalpyVec(grid))
vecs.add(PISM.model.createIceMaskVec(grid))

# Read in the PISM state variables that are used directly in the SSA solver
for v in [vecs.thk, vecs.topg, vecs.enthalpy]:
    v.regrid(input_file, critical=True)

# variables mask and surface are computed from the geometry previously read
sea_level = PISM.model.createSeaLevelVec(grid)
sea_level.set(0.0)
gc = PISM.GeometryCalculator(config)
gc.compute(sea_level, vecs.topg, vecs.thk, vecs.mask, vecs.surface_altitude)

# If running in regional mode, load in regional variables
if is_regional:
    vecs.add(PISM.model.createNoModelMask(grid))
    vecs.no_model_mask.regrid(input_file, critical=True)

    if PISM.util.fileHasVariable(input_file, 'usurfstore'):
        vecs.add(PISM.model.createIceSurfaceStoreVec(grid))
        vecs.usurfstore.regrid(input_file, critical=True)
    else:
        vecs.add(vecs.surface, 'usurfstore')

    solver = PISM.SIAFD_Regional
else:
    solver = PISM.SIAFD

PISM.verbPrintf(2, context.com, "* Computing SIA velocities...\n")
vel_sia = PISM.sia.computeSIASurfaceVelocities(modeldata, siasolver=solver)

PISM.verbPrintf(2, context.com, "* Saving results to %s...\n" % output_file)
pio = PISM.util.prepare_output(output_file)
pio.close()

# Save time & command line & results
PISM.util.writeProvenance(output_file)
vel_sia.write(output_file)
