#! /usr/bin/env python
#
# Copyright (C) 2011, 2012 David Maxwell
# 
# This file is part of PISM.
# 
# PISM is free software; you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation; either version 2 of the License, or (at your option) any later
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
config = context.config

PISM.set_abort_on_sigint(True)

usage = \
"""  sia_forward.py -i IN.nc [-o file.nc]
  where:
    -i      IN.nc is input file in NetCDF format: contains PISM-written model state
  notes:
    * -i is required
  """

PISM.verbosityLevelFromOptions()
PISM.show_usage_check_req_opts(context.com,"sia.py",["-i"],usage)

for o in PISM.OptionsGroup(context.com,"","sia"):
  input_file = PISM.optionsString("-i","input file")
  output_file = PISM.optionsString("-o","output file",default="sia_"+os.path.basename(input_file))
  is_regional = PISM.optionsFlag("-regional","Compute SIA using regional model semantics",default=False)
  verbosity = PISM.optionsInt("-verbose","verbosity level",default=2)
  PISM.set_config_from_options(context.com,config)


periodicity = PISM.XY_PERIODIC
if is_regional:
  periodicity=PISM.NOT_PERIODIC
grid = PISM.Context().newgrid()
PISM.util.init_grid_from_file(grid,input_file,periodicity);

basal = PISM.IceBasalResistancePlasticLaw(config)

enthalpyconverter = PISM.EnthalpyConverter(config)
if PISM.getVerbosityLevel() >3:
  enthalpyconverter.viewConstants(PETSc.Viewer.STDOUT())


modeldata = PISM.model.ModelData(grid)
modeldata.setPhysics(basal,enthalpyconverter)

vecs = modeldata.vecs;

vecs.add( PISM.util.standardIceSurfaceVec( grid ), 'surface' )
vecs.add( PISM.util.standardIceThicknessVec( grid ), 'thickness' )
vecs.add( PISM.util.standardBedrockElevationVec( grid ), 'bed' )
vecs.add( PISM.util.standardEnthalpyVec( grid ), 'enthalpy' )
vecs.add( PISM.util.standardIceMask( grid ), 'ice_mask' )

# Read in the PISM state variables that are used directly in the SSA solver
for v in [vecs.thickness, vecs.bed, vecs.enthalpy]:
  v.regrid(input_file,critical=True)

# variables mask and surface are computed from the geometry previously read
sea_level = 0 # FIXME setFromOption?
gc = PISM.GeometryCalculator(sea_level, config)
gc.compute(vecs.bed,vecs.thickness,vecs.ice_mask,vecs.surface)

# If running in regional mode, load in regional variables
if is_regional:
  vecs.add( PISM.util.standardNoModelMask( grid ), 'no_model_mask' )
  vecs.no_model_mask.regrid(input_file,critical=True)
  
  if PISM.util.fileHasVariable(input_file,'usurfstore'):
    vecs.add( PISM.util.standardIceSurfaceStoreVec(grid) )
    vecs.usurfstore.regrid(input_file,critical=True)
  else:        
    vecs.add( vecs.surface, 'usurfstore')
    vecs.setPISMVarsName('usurfstore','usurfstore')
    
  solver = PISM.SIAFD_Regional
else:
  solver = PISM.SIAFD


vel_sia =   PISM.sia.computeSIASurfaceVelocities(modeldata,siasolver=solver)


pio = PISM.PIO(grid.com,grid.rank,"netcdf3")
pio.open(output_file,PISM.NC_WRITE,False)
pio.def_time(grid.config.get_string("time_dimension_name"),
             grid.config.get_string("calendar"), grid.time.units())
pio.append_time(grid.config.get_string("time_dimension_name"),grid.time.current())
pio.close()

# Save time & command line & results
PISM.util.writeProvenance(output_file)
vel_sia.write(output_file)
