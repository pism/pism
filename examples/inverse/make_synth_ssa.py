#! /usr/bin/env python
#
# Copyright (C) 2011 David Maxwell and Constantine Khroulev
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


import sys, petsc4py
petsc4py.init(sys.argv)
from petsc4py import PETSc
import PISM, time, math
import pismssaforward

tauc_prior_scale = 0.2
tauc_prior_const = None


def addGroundedIceMisfitWeight(modeldata):
  grid = modeldata.grid
  weight = PISM.util.standardVelocityMisfitWeight(grid)
  mask = modeldata.vecs.ice_mask
  with PISM.util.Access(comm=weight,nocomm=mask):
    weight.set(0.)
    grounded = PISM.MASK_GROUNDED
    for (i,j) in grid.points():
      if mask[i,j] == grounded:
        weight[i,j] = 1
  modeldata.vecs.add(weight,writing=True)    

# The main code for a run follows:
if __name__ == '__main__':
  context = PISM.Context()
  com = context.com

  PISM.set_abort_on_sigint(True)

  PISM.verbosityLevelFromOptions()
  PISM.verbPrintf(2,PISM.Context().com,"SSA forward model.\n")
  PISM.stop_on_version_option()
  usage = \
"""  ssa_forward.py -i IN.nc -Mx number -My number [-o file.nc]
  or (at python prompt)
    run ssa_forward -i IN.nc -Mx number -My number [-o file.nc]
  where:
    -i      IN.nc is input file in NetCDF format: contains PISM-written model state
    -Mx     number of grid points in the x direction
    -My     number of grid points in the y direction
  notes:
    * -i is required
  """

  PISM.show_usage_check_req_opts(com,"ssa_forward",["-i"],usage)

  config = context.config
  for o in PISM.OptionsGroup(com,"","SSA Forward"):
    input_file_name = PISM.optionsString("-i","file to bootstrap from")
    output_file_name = PISM.optionsString("-o","output file",default="make_synth_ssa.nc")
    tauc_prior_scale = PISM.optionsReal("-tauc_prior_scale","initial guess for tauc to be this factor of the true value",default=tauc_prior_scale)
    tauc_prior_const = PISM.optionsReal("-tauc_prior_const","initial guess for tauc to be this constant",default=tauc_prior_const)
    noise = PISM.optionsReal("-rms_noise","pointwise rms noise to add (in m/a)",default=None)
    generate_ssa_observed = PISM.optionsFlag("-generate_ssa_observed","generate observed SSA velocities",default=False)

  config.set_string("ssa_method","fem")
  
  ssa_run = PISM.ssa.SSAFromBootFile(input_file_name)

  ssa_run.setup()

  solve_t0 = time.clock()
  ssa_run.solve()
  solve_t = time.clock()-solve_t0

  PISM.verbPrintf(2,context.com,"Solve time %g seconds.\n",solve_t)

  modeldata = ssa_run.modeldata
  grid = modeldata.grid
  vecs = modeldata.vecs

  # Add the misfit weight.
  addGroundedIceMisfitWeight(modeldata)

  # Generate a prior guess for tauc
  tauc_prior = PISM.util.standardYieldStressVec(grid,name='tauc_prior',desc="initial guess for (pseudo-plastic) basal yield stress in an inversion")
  vecs.add(tauc_prior,writing=True)
  if not tauc_prior_const is None:
    tauc_prior.set(tauc_prior_const)
  else:
    tauc_prior.copy_from(modeldata.vecs.tauc)
    tauc_prior.scale(tauc_prior_scale)
  
  tauc_true = modeldata.vecs.tauc
  tauc_true.set_name('tauc_true')
  tauc_true.set_attrs("diagnostic", "value of basal yield stress used to generate synthetic SSA velocities", "Pa", ""); 
  vecs.markForWriting(tauc_true)

  vel_ssa_observed = vecs.vel_ssa_blat
  vel_ssa_observed.rename("_ssa_observed","'observed' SSA velocities'","")
  if generate_ssa_observed:
    vecs.markForWriting(vel_ssa_observed)
    final_velocity = vel_ssa_observed
  else:
    vel_sia_observed = pismssaforward.computeSIASurfaceVelocities(modeldata)
    vel_ssa_observed.rename("_sia_observed","'observed' SIA velocities'","")
    vel_surface_observed = PISM.util.standard2dVelocityVec(grid,"_surface_observed","observed surface velocities",stencil_width=1)
    vel_surface_observed.copy_from(vel_sia_observed)
    vel_surface_observed.add(1.,vel_ssa_observed)
    vecs.markForWriting(vel_surface_observed)
    final_velocity = vel_surface_observed

  if not noise is None:
    u_noise = pismssaforward.randVectorV(grid,noise/math.sqrt(2),final_velocity.get_stencil_width())
    final_velocity.add(1./PISM.secpera,u_noise)

  pio = PISM.PISMIO(grid)
  pio.open_for_writing(output_file_name,False,True)
  pio.append_time(grid.config.get_string("time_dimension_name"),0.0)
  pio.close()

  vecs.write(output_file_name)
  
  # Save time & command line
  PISM.util.writeProvenance(output_file_name)
