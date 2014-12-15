#! /usr/bin/env python
#
# Copyright (C) 2011, 2012, 2013, 2014 David Maxwell and Constantine Khroulev
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


import PISM, sys, time, math

design_prior_scale = 0.2
design_prior_const = None

def groundedIceMisfitWeight(modeldata):
  grid = modeldata.grid
  weight = PISM.model.createVelocityMisfitWeightVec(grid)
  mask = modeldata.vecs.ice_mask
  with PISM.vec.Access(comm=weight,nocomm=mask):
    weight.set(0.)
    grounded = PISM.MASK_GROUNDED
    for (i,j) in grid.points():
      if mask[i,j] == grounded:
        weight[i,j] = 1
  return weight


def fastIceMisfitWeight(modeldata,threshold):
  grid = modeldata.grid
  weight = PISM.model.createVelocityMisfitWeightVec(grid)
  mask    = modeldata.vecs.ice_mask
  vel_ssa = modeldata.vecs.vel_ssa
  threshold = threshold*threshold
  with PISM.vec.Access(comm=weight,nocomm=[vel_ssa,mask]):
    weight.set(0.)
    grounded = PISM.MASK_GROUNDED
    for (i,j) in grid.points():
      u=vel_ssa[i,j].u; v=vel_ssa[i,j].v
      if mask[i,j] == grounded:
        if u*u+v*v > threshold:
          weight[i,j] = 1
  return weight


# The main code for a run follows:
if __name__ == '__main__':
  context = PISM.Context()
  com = context.com

  PISM.set_abort_on_sigint(True)

  PISM.verbosityLevelFromOptions()
  PISM.verbPrintf(2,PISM.Context().com,"SSA forward model.\n")
  PISM.stop_on_version_option()
  usage = \
"""  %s -i IN.nc -Mx number -My number [-o file.nc]
  or (at python prompt)
    run %s -i IN.nc -Mx number -My number [-o file.nc]
  where:
    -i      IN.nc is input file in NetCDF format: contains PISM-written model state
    -Mx     number of grid points in the x direction
    -My     number of grid points in the y direction
  notes:
    * -i is required
  """ % (sys.argv[0], sys.argv[0])

  PISM.show_usage_check_req_opts(com, sys.argv[0], ["-i"], usage)

  config = context.config
  if not PISM.optionsIsSet("-ssa_method"):
    config.set_string("ssa_method","fem")

  for o in PISM.OptionsGroup(com,"","SSA Forward"):
    input_file_name = PISM.optionsString("-i","file to bootstrap from")
    output_file_name = PISM.optionsString("-o","output file",default="make_synth_ssa.nc")
    design_prior_scale = PISM.optionsReal("-design_prior_scale","initial guess for design variable to be this factor of the true value",default=design_prior_scale)
    design_prior_const = PISM.optionsReal("-design_prior_const","initial guess for design variable to be this constant",default=design_prior_const)
    noise = PISM.optionsReal("-rms_noise","pointwise rms noise to add (in m/a)",default=None)
    misfit_weight_type = PISM.optionsList(context.com,"-misfit_type","Choice of misfit weight function",["grounded","fast"],"grounded")
    fast_ice_speed = PISM.optionsReal("-fast_ice_speed","Threshold in m/a for determining if ice is fast", 500.)
    generate_ssa_observed = PISM.optionsFlag("-generate_ssa_observed","generate observed SSA velocities",default=False)
    is_regional = PISM.optionsFlag("-regional","Compute SIA/SSA using regional model semantics",default=False)
    design_var = PISM.optionsList(context.com,"-inv_ssa","design variable for inversion", ["tauc", "hardav"], "tauc")
    
  
  ssa_run = PISM.ssa.SSAFromInputFile(input_file_name)

  ssa_run.setup()

  solve_t0 = time.clock()
  ssa_run.solve()
  solve_t = time.clock()-solve_t0

  PISM.verbPrintf(2,context.com,"Solve time %g seconds.\n",solve_t)

  modeldata = ssa_run.modeldata
  grid = modeldata.grid
  vecs = modeldata.vecs

  if design_var == 'tauc':
    # Generate a prior guess for tauc
    tauc_prior = PISM.model.createYieldStressVec(grid,name='tauc_prior',desc="initial guess for (pseudo-plastic) basal yield stress in an inversion")
    vecs.add(tauc_prior,writing=True)
    if design_prior_const is not None:
      tauc_prior.set(design_prior_const)
    else:
      tauc_prior.copy_from(modeldata.vecs.tauc)
      tauc_prior.scale(design_prior_scale)
  
    tauc_true = modeldata.vecs.tauc
    tauc_true.set_name('tauc_true')
    tauc_true.set_attrs("diagnostic", "value of basal yield stress used to generate synthetic SSA velocities", "Pa", ""); 
    vecs.markForWriting(tauc_true)
  elif design_var == 'hardav':
    # Generate a prior guess for hardav

    EC = PISM.EnthalpyConverter(config)
    ice_factory = PISM.IceFlowLawFactory(grid.com, "ssa_", config, EC);
    ice_factory.removeType(PISM.ICE_GOLDSBY_KOHLSTEDT);
    ice_factory.setType(config.get_string("ssa_flow_law")); 
    ice_factory.setFromOptions();
    flow_law = ice_factory.create()
    hardav = PISM.model.createAveragedHardnessVec(grid)
    flow_law.averaged_hardness_vec(vecs.thickness, vecs.enthalpy, hardav)
    vecs.add(hardav)

    hardav_prior = PISM.model.createAveragedHardnessVec(grid,name='hardav_prior',desc="initial guess for vertically averaged ice hardness in an inversion")
    vecs.add(hardav_prior,writing=True)
    if design_prior_const is not None:
      hardav_prior.set(design_prior_const)
    else:
      hardav_prior.copy_from(modeldata.vecs.hardav)
      hardav_prior.scale(hardav_prior_scale)
  
    hardav_true = modeldata.vecs.hardav
    hardav_true.set_name('hardav_true')
    hardav_true.set_attrs("diagnostic", "vertically averaged ice hardness used to generate synthetic SSA velocities", "Pa s^0.33333", ""); 
    vecs.markForWriting(hardav_true)

  vel_ssa_observed = vecs.vel_ssa
  vel_ssa_observed.rename("_ssa_observed","'observed' SSA velocities'","")
  if generate_ssa_observed:
    vecs.markForWriting(vel_ssa_observed)
    final_velocity = vel_ssa_observed
  else:
    sia_solver = PISM.SIAFD
    if is_regional:
      sia_solver = PISM.SIAFD_Regional
    vel_sia_observed = PISM.sia.computeSIASurfaceVelocities(modeldata,sia_solver)
    vel_ssa_observed.rename("_sia_observed","'observed' SIA velocities'","")
    vel_surface_observed = PISM.model.create2dVelocityVec(grid,"_surface_observed","observed surface velocities",stencil_width=1)
    vel_surface_observed.copy_from(vel_sia_observed)
    vel_surface_observed.add(1.,vel_ssa_observed)
    vecs.markForWriting(vel_surface_observed)
    final_velocity = vel_surface_observed

  # Add the misfit weight.
  if misfit_weight_type == "fast":
    misfit_weight = fastIceMisfitWeight(modeldata,
                                        grid.convert(fast_ice_speed, "m/year", "m/second"))
  else:
    misfit_weight = groundedIceMisfitWeight(modeldata)
  modeldata.vecs.add(misfit_weight,writing=True)    

  if not noise is None:
    u_noise = PISM.vec.randVectorV(grid,noise/math.sqrt(2),final_velocity.get_stencil_width())
    final_velocity.add(grid.convert(1.0, "m/year", "m/second"),
                       u_noise)

  pio = PISM.PIO(grid.com, "netcdf3", grid.config.get_unit_system())
  pio.open(output_file_name, PISM.PISM_READWRITE_MOVE)
  pio.def_time(grid.config.get_string("time_dimension_name"),
               grid.config.get_string("calendar"), grid.time.units_string())
  pio.append_time(grid.config.get_string("time_dimension_name"),grid.time.current())
  pio.close()

  vecs.write(output_file_name)
  
  # Save time & command line
  PISM.util.writeProvenance(output_file_name)
