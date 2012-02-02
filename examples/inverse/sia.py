#! /usr/bin/env python
#
# Copyright (C) 2011 David Maxwell
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

import sys, os, petsc4py
petsc4py.init(sys.argv)
from petsc4py import PETSc
import PISM

context = PISM.Context()
config = context.config

PISM.set_abort_on_sigint(True)

usage = \
"""  sia.py -i IN.nc [-o file.nc]
  where:
    -i      IN.nc is input file in NetCDF format: contains PISM-written model state
  notes:
    * -i is required
  """

PISM.show_usage_check_req_opts(context.com,"sia.py",["-i"],usage)

for o in PISM.OptionsGroup(context.com,"","sia"):
  bootfile = PISM.optionsString("-i","input file")
  output_file = PISM.optionsString("-o","output file",default="sia_"+os.path.basename(bootfile))

  verbosity = PISM.optionsInt("-verbose","verbosity level",default=2)
  PISM.set_config_from_options(context.com,config)


grid = PISM.Context().newgrid()
PISM.util.init_grid_from_file(grid,bootfile,
                                periodicity=PISM.XY_PERIODIC);

enthalpyconverter = PISM.EnthalpyConverter(config)
if PISM.getVerbosityLevel() >3:
  enthalpyconverter.viewConstants(PETSc.Viewer.STDOUT())

if PISM.optionsIsSet("-ssa_glen"):
  ice = PISM.CustomGlenIce(com,"",config,enthalpyconverter)
  B_schoof = 3.7e8;     # Pa s^{1/3}; hardness 
  ice.setHardness(B_schoof)
else:
  ice =  PISM.GPBLDIce(grid.com, "", config,enthalpyconverter)
ice.setFromOptions()

surface    = PISM.util.standardIceSurfaceVec( grid )
thickness  = PISM.util.standardIceThicknessVec( grid )
bed        = PISM.util.standardBedrockElevationVec( grid )
tauc       = PISM.util.standardYieldStressVec( grid )
enthalpy   = PISM.util.standardEnthalpyVec( grid )
ice_mask   = PISM.util.standardIceMask( grid )
v = [surface,thickness,bed,tauc,enthalpy,ice_mask]
pv =PISM.PISMVars()
for var in v:
  var.regrid(bootfile,True)
  pv.add(var)

sia = PISM.SIAFD(grid,ice,enthalpyconverter,config)
sia.init(pv)

zero_sliding = PISM.IceModelVec2V()
zero_sliding.create(grid, 'basal_velocity', False )
zero_sliding.set(0.)

zero_D2 = PISM.IceModelVec2S();
zero_D2.create(grid, 'D2', True, PISM.util.WIDE_STENCIL)
zero_D2.set(0)

sia.update(zero_sliding,zero_D2,False)
(u,v) = sia.get_horizontal_3d_velocity()

U_s = PISM.util.standard2dVelocityVec(grid,name="_sia_surf")
tmp = PISM.IceModelVec2S()
tmp.create(grid,'tmp',False)
u.getSurfaceValues(tmp,thickness)
U_s.set_component(0,tmp)
v.getSurfaceValues(tmp,thickness)
U_s.set_component(1,tmp)

pio = PISM.PIO(grid.com,grid.rank,"netcdf3")
pio.open(output_file,PISM.NC_WRITE,False)
pio.def_time(grid.config.get_string("time_dimension_name"),
             grid.config.get_string("calendar"), grid.time.units())
pio.append_time(grid.config.get_string("time_dimension_name"),grid.time.current())
pio.close()

# Save time & command line
PISM.util.writeProvenance(output_file)
U_s.write(output_file)
