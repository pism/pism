#! /usr/bin/env python
#

import sys, petsc4py
petsc4py.init(sys.argv)
from petsc4py import PETSc

import PISM, math

# Default constants that  may get overridden later.

Ly = 25e3   #  25 km
Lx = 50e3   #  50 km
Lz = 4000
My = 13
Mx = 23
Mz = 21

sea_level = 0;      # m sea level elevation
H0 = 60.;           # ice thickness at cliff
alpha = 0.008       # constant surface slope
Lext = 15e3          # width of strip beyond cliff
Lstream_x = 50e3
Lstream_y = 30e3

Hext = 0.    # m ice thickeness beyond the cliff

tauc_hi = 2e6       # Pa
tauc_lo = 1e4       # Pa
tauc_free_bedrock = 0  # Will get set later

enth0  = 528668.35; # Hmmm. 263.15 Kelvin at depth=0.
bed0  = 0;

def geometry(x,y):
  x0=-Lx+Lext;
  if x<x0:
    return (0,Hext)
  return (0,H0+alpha*(x-x0))

def stream_tauc(x,y):
  x0 = -Lx+Lext
  if x<x0:
    return tauc_free_bedrock
  if x<x0+Lstream_x:
    if abs(y)<Lstream_y/2:
      return tauc_lo
  return tauc_hi
      

# The main code for a run follows:
if __name__ == '__main__':
  PISM.set_abort_on_sigint(True)
  context = PISM.Context()
  for o in PISM.OptionsGroup(context.com,"","build_tiny"):
    Mx = PISM.optionsInt("-Mx","Number of grid points in x-direction",default=Mx)
    My = PISM.optionsInt("-My","Number of grid points in y-direction",default=My)
    output_filename = PISM.optionsString("-o","output file",default="tiny.nc")
    verbosity = PISM.optionsInt("-verbose","verbosity level",default=3)

  # Build the grid.
  grid = PISM.Context().newgrid()
  config = grid.config
  PISM.util.init_grid(grid,Lx,Ly,Lz,Mx,My,Mz,PISM.NOT_PERIODIC)
  vecs = PISM.model.ModelVecs();
  vecs.add( PISM.util.standardIceSurfaceVec( grid ), 'surface')
  vecs.add( PISM.util.standardIceThicknessVec( grid ), 'thickness')
  vecs.add( PISM.util.standardBedrockElevationVec(grid), 'bed')
  vecs.add( PISM.util.standardYieldStressVec( grid ), 'tauc')
  vecs.add( PISM.util.standardEnthalpyVec( grid ), 'enthalpy' )
  vecs.add( PISM.util.standardIceMask( grid ), 'ice_mask' )
  vecs.add( PISM.util.standardNoModelMask( grid ), 'no_model_mask' )
  vecs.add( PISM.util.standard2dVelocityVec( grid, name='_ssa_bc',desc='SSA Dirichlet BC') )

  # Set constant coefficients.
  vecs.enthalpy.set(enth0)

  # Build the continent
  bed = vecs.bed
  thickness = vecs.thickness

  with PISM.util.Access(comm=[bed,thickness]):
    for (i,j) in grid.points():
      x=grid.x[i]; y=grid.y[j]
      (b,t) = geometry(x,y)
      bed[i,j]=b; thickness[i,j]=t;

  # Compute mask and surface elevation from geometry variables.
  gc = PISM.GeometryCalculator(sea_level,grid.config)
  gc.compute(bed,thickness,vecs.ice_mask,vecs.surface)

  tauc = vecs.tauc; mask = vecs.ice_mask
  tauc_free_bedrock = config.get('high_tauc')
  with PISM.util.Access(comm=tauc,nocomm=mask):
    for (i,j) in grid.points():
      x=grid.x[i]; y=grid.y[j]
      tauc[i,j] = stream_tauc(x,y)

  vecs.vel_ssa_bc.set(0.)
  no_model_mask = vecs.no_model_mask
  no_model_mask.set(0)
  with PISM.util.Access(comm=[no_model_mask]):
    for (i,j) in grid.points():
      if (i==0) or (i==grid.Mx-1) or (j==0) or (j==grid.My-1):
        no_model_mask[i,j] = 1
  

  pio = PISM.PIO(grid.com, grid.rank, "netcdf3")
  pio.open(output_filename, PISM.NC_WRITE)
  pio.def_time(grid.config.get_string("time_dimension_name"),
               "365_day", "seconds since 1-1-1")
  pio.append_time(grid.config.get_string("time_dimension_name"),0.0)
  pio.close()
  vecs.writeall(output_filename)
  PISM.util.writeProvenance(output_filename)
