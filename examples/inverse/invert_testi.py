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
import PISM.invert_ssa
import PISM.sipletools
from PISM import util

import numpy as np
import math

import siple
siple.reporting.clear_loggers()
siple.reporting.add_logger(PISM.sipletools.pism_logger)
siple.reporting.set_pause_callback(PISM.sipletools.pism_pause)


import matplotlib.pyplot as pp


class PlotListener(PISM.invert_ssa.PlotListener):
  def __call__(self,inv_solver,it,data):

    grid = self.grid
        
    zeta = self.toproczero(data.zeta)
    u    = self.toproczero(data.u)
    if inv_solver.method.startswith('tikhonov'):
      eta = data.eta
      sWeight=1
      dWeight=1/eta
      grad_zeta = self.toproczero(data.grad_zeta)
      grad_u = self.toproczero(data.grad_u)
      grad= self.toproczero(data.grad)
    else:
      r = self.toproczero(data.r)

    if zeta is not None:
      pp.figure(self.figure())
      pp.clf()
      pp.subplot(1,3,1)
      pp.plot(grid.y,zeta[:,grid.Mx/2])
    
      mag = np.max(np.abs(u))
      pp.subplot(1,3,2)
      pp.plot(grid.y,u[0,:,grid.Mx/2]/mag)

      pp.subplot(1,3,3)
      if inv_solver.method.startswith('tikhonov'):
        pp.plot(grid.y,-grad_zeta[:,grid.Mx/2]*dWeight,grid.y,grad_u[:,grid.Mx/2]*sWeight,grid.y,grad[:,grid.Mx/2])
      else:
        pp.plot(grid.y,r[0,:,grid.Mx/2])
        
      pp.ion()
      pp.show()


class LinPlotListener(PISM.invert_ssa.PlotListener):
  def __call__(self,inv_solver,it,data):

    grid = self.grid

    x  = self.toproczero(data.x)

    if x is not None:
      pp.figure(self.figure())
      pp.clf()
      mag = np.max(np.abs(x));
      if mag==0 : mag = 1
      pp.plot(grid.y,x[:,grid.Mx/2]/mag)

      pp.ion()
      pp.show()

Mx = 11 
My = 61

kFEMStencilWidth = 1

m_schoof = 2;        # (pure number)
L_schoof = 40e3;      # meters
aspect_schoof = 0.05; # (pure)
H0_schoof = aspect_schoof * L_schoof; 
                      # = 2000 m THICKNESS
B_schoof = 3.7e8;     # Pa s^{1/3}; hardness 
                      # given on p. 239 of Schoof; why so big?
p_schoof = 4.0/3.0;   # = 1 + 1/n

slope = 0.001

rms_error = 1. #m/a

right_side_weight = 1.

tauc_guess_scale = 0.3
tauc_guess_const = None

def testi_tauc(grid, tauc):
  standard_gravity = grid.config.get("standard_gravity");
  ice_density = grid.config.get("ice_density");
  f = ice_density * standard_gravity * H0_schoof * slope

  with PISM.util.Access(comm=tauc):
    for (i,j) in grid.points():
      y=grid.y[j]
      tauc[i,j] = f* (abs(y/L_schoof)**m_schoof)

class testi_run(PISM.invert_ssa.InvSSARun):
  def __init__(self,Mx,My):
    self.grid = PISM.Context().newgrid()
    self.Mx = Mx
    self.My = My

  def _initGrid(self):
    Mx=self.Mx; My=self.My
    Ly = 3*L_schoof   # 300.0 km half-width (L=40.0km in Schoof's choice of variables)
    Lx = max(60.0e3, ((Mx - 1) / 2.) * (2.0 * Ly / (My - 1)) )
    PISM.util.init_shallow_grid(self.grid,Lx,Ly,Mx,My,PISM.X_PERIODIC);

  def _initPhysics(self):
    config = self.config
    basal = PISM.IceBasalResistancePlasticLaw(config)

    # irrelevant
    enthalpyconverter = PISM.EnthalpyConverter(config);

    config.set_string("ssa_flow_law", "isothermal_glen")
    config.set("ice_softness", pow(3.7e8, -config.get("Glen_exponent")))

    self.modeldata.setPhysics(basal,enthalpyconverter)


  def _initSSACoefficients(self):
    vecs = self.modeldata.vecs; grid = self.grid
    vecs.add( util.standardIceThicknessVec( grid ), 'thickness')
    vecs.add( util.standardBedrockElevationVec(grid), 'bed')
    vecs.add( util.standardYieldStressVec( grid ), 'tauc')
    vecs.add( util.standardEnthalpyVec( grid ), 'enthalpy' )
    vecs.add( util.standardIceMask( grid ), 'ice_mask' )
    vecs.add( util.standardDrivingStressX( grid ) )
    vecs.add( util.standardDrivingStressY( grid ) )
    vecs.add( util.standardVelocityMisfitWeight(grid) )

    self._allocateBCs()

    vecs.thickness.set(H0_schoof)
    vecs.ice_mask.set(PISM.MASK_GROUNDED)
    vecs.bed.set(0.)

    testi_tauc(self.modeldata.grid, vecs.tauc)

    grid = self.grid
    standard_gravity = grid.config.get("standard_gravity");
    ice_density = grid.config.get("ice_density");
    f = ice_density * standard_gravity * H0_schoof * slope

    vecs.ssa_driving_stress_y.set(0)
    vecs.ssa_driving_stress_x.set(f)
    
    with PISM.util.Access(comm=[vecs.bc_mask,vecs.vel_bc]):
      for (i,j) in grid.points():
        if (j == 0) or (j==grid.My-1):
          vecs.bc_mask[i,j]=1
          vecs.vel_bc[i,j].u=0
          vecs.vel_bc[i,j].v=0

    misfit_weight=vecs.vel_misfit_weight
    with PISM.util.Access(comm=misfit_weight):
      for (i,j) in grid.points():
        if grid.y[j] <= 0:
          misfit_weight[i,j] = 1.;
        else:
          misfit_weight[i,j] = right_side_weight;

## Main code starts here
if __name__ == "__main__":
  context = PISM.Context()
  config = context.config
  PISM.set_abort_on_sigint(True)

  for o in PISM.OptionsGroup(context.com,"","Invert test I"):
    Mx = PISM.optionsInt("-Mx","Number of grid points in x-direction",default=Mx)
    My = PISM.optionsInt("-My","Number of grid points in y-direction",default=My)
    output_file = PISM.optionsString("-o","output file",default="invert_testi.nc")
    verbosity = PISM.optionsInt("-verbose","verbosity level",default=2)
    right_side_weight = PISM.optionsReal("-right_side_weight","L2 weight for y>0",default=right_side_weight)
    tauc_guess_scale = PISM.optionsReal("-tauc_guess_scale","initial guess for tauc to be this factor of the true value",default=tauc_guess_scale)
    tauc_guess_const = PISM.optionsReal("-tauc_guess_const","initial guess for tauc to be this constant",default=tauc_guess_const)
    do_plotting = PISM.optionsFlag("-inv_plot","perform visualization during the computation",default=False)
    do_final_plot = PISM.optionsFlag("-inv_final_plot","perform visualization at the end of the computation",default=True)
    do_pause = PISM.optionsFlag("-inv_pause","pause each iteration",default=False)
    test_adjoint = PISM.optionsFlag("-inv_test_adjoint","Test that the adjoint is working",default=False)
    monitor_adjoint = PISM.optionsFlag("-inv_monitor_adjoint","Track accuracy of the adjoint during computation",default=False)

  inv_method = config.get_string("inv_ssa_method")

  length_scale  = L_schoof
  slope = 0.001
  standard_gravity = config.get("standard_gravity");
  ice_density = config.get("ice_density");
  f0 = ice_density * standard_gravity * H0_schoof * slope
  stress_scale  = f0
  Ly = 3*L_schoof   # 300.0 km half-width (L=40.0km in Schoof's choice of variables)
  Lx = max(60.0e3, ((Mx - 1) / 2.) * (2.0 * Ly / (My - 1)) )
  area_scale    = Lx*Ly
  depth_scale   = H0_schoof

  B = B_schoof
  velocity_scale   = (f0/B)**(3.)*(length_scale/depth_scale)**(3.)*length_scale
  time_scale       = length_scale / velocity_scale 
  strainrate_scale = 1./ time_scale
  viscosity_scale  = B*(strainrate_scale**(-2./3.))
  nuH_scale        = viscosity_scale * depth_scale

  PISM.setVerbosityLevel(verbosity)
  testi = testi_run(Mx,My)
  testi.setup()
  solver = PISM.invert_ssa.InvSSASolver(testi)
  tauc_param = solver.ssarun.tauc_param

  grid = testi.grid

  # Build the true yeild stress for test I
  tauc_true = PISM.util.standardYieldStressVec(grid,name="tauc_true")
  testi_tauc(grid, tauc_true)

  # Convert tauc_true to zeta_true
  zeta_true = PISM.IceModelVec2S();
  zeta_true.create(grid,"zeta_true",PISM.kHasGhosts,kFEMStencilWidth)
  tauc_param = PISM.invert_ssa.tauc_param_factory.create(config)
  tauc_param.convertFromTauc(tauc_true,zeta_true)

  # Build the initial guess for tauc for the inversion.
  tauc = PISM.util.standardYieldStressVec(grid)
  if not tauc_guess_const is None:
    tauc.set(tauc_guess_const)
  else:
    testi_tauc(grid, tauc)
    tauc.scale(tauc_guess_scale)

  # Convert tauc guess to zeta guess
  zeta0 = PISM.IceModelVec2S();
  zeta0.create(grid, "zeta", PISM.kHasGhosts, kFEMStencilWidth)
  tauc_param.convertFromTauc(tauc,zeta0)

  if test_adjoint:
    if solver.method.startswith('tikhonov'):
      siple.reporting.msg("option -inv_test_adjoint cannot be used with inverse method %s",solver.method)
      exit(1)
    from PISM.sipletools import PISMLocalVector as PLV
    stencil_width=1
    forward_problem = solver.forward_problem
    d = PLV(PISM.sipletools.randVectorS(grid,1e5,stencil_width))
    r = PLV(PISM.sipletools.randVectorV(grid,1./PISM.secpera,stencil_width))
    (domainIP,rangeIP)=forward_problem.testTStar(PLV(zeta0),d,r,3)
    siple.reporting.msg("domainip %g rangeip %g",domainIP,rangeIP)
    exit(0)

  # Setup the output file.
  pio = PISM.PIO(grid.com,grid.rank,"netcdf3")
  pio.open(output_file,PISM.NC_WRITE,False)
  pio.def_time(grid.config.get_string("time_dimension_name"),
               grid.config.get_string("calendar"), grid.time.units())
  pio.append_time(grid.config.get_string("time_dimension_name"),grid.time.current())
  pio.close()
  zeta0.write(output_file)

  # Send the true yeild stress through the forward problem to 
  # get at true velocity field.
  u_obs = PISM.util.standard2dVelocityVec( grid, name='_ssa_true', desc='SSA velocity boundary condition',intent='intent' )
  solver.solveForward(zeta_true,out=u_obs)

  # Attach various iteration listeners to the solver as needed for:
  # progress reporting,
  if inv_method.startswith('tikhonov'):
    solver.addIterationListener(PISM.invert_ssa.printTikhonovProgress)
  # adjoint monitoring,
  if monitor_adjoint:    
    solver.addIterationListener(PISM.invert_ssa.MonitorAdjoint())
      # if inv_method=='ign':
      #   solver.addLinearIterationListener(MonitorAdjointLin())
  # Plotting
  if do_plotting:
    solver.addIterationListener(PlotListener(grid))
    if inv_method=='ign':
      solver.addLinearIterationListener(LinPlotListener(grid))
  # Pausing
  if do_pause:
    solver.addIterationListener(PISM.invert_ssa.pauseListener)            
  # Iteration saving
  solver.addXUpdateListener(PISM.invert_ssa.ZetaSaver(output_file)) 

  # Try solving
  reason = solver.solveInverse(zeta0,u_obs,zeta0)
  if reason.failed():
    PISM.verbPrintf(1,grid.com,"Inverse solve FAILURE (%s)!\n" % reason.description());
    quit()
  PISM.verbPrintf(1,grid.com,"Inverse solve success (%s)!\n" % reason.description());

  (zeta_i,u_i) = solver.inverseSolution()

  tauc_param.convertToTauc(zeta_i,tauc)

  # Write solution out to netcdf file
  testi.write(output_file)
  tauc.write(output_file)
  tauc_true.write(output_file)
  u_i.set_name("_computed",0)
  u_i.write(output_file)
  u_obs.write(output_file)

  # Draw a pretty picture
  tz = PISM.toproczero.ToProcZero(grid)
  tauc_a = tz.communicate(tauc)
  tauc_true = tz.communicate(tauc_true)
  tz2 = PISM.toproczero.ToProcZero(grid,dof=2,dim=2)
  u_i_a = tz2.communicate(u_i)
  u_obs_a = tz2.communicate(u_obs)

  if do_final_plot and (not tauc_a is None):
    from matplotlib import pyplot
    pyplot.clf()
    pyplot.subplot(1,2,1)
    pyplot.plot(grid.y,tauc_a[:,Mx/2])
    pyplot.plot(grid.y,tauc_true[:,Mx/2])

    pyplot.subplot(1,2,2)
    pyplot.plot(grid.y,u_i_a[0,:,Mx/2]*PISM.secpera)
    pyplot.plot(grid.y,u_obs_a[0,:,Mx/2]*PISM.secpera)

    pyplot.ion()
    pyplot.show()
    siple.reporting.endpause()
