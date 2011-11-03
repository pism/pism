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

import sys, petsc4py
petsc4py.init(sys.argv)
from petsc4py import PETSc
import numpy as np
import tozero
import siple
import math

import PISM

from pismssaforward import InvSSARun, SSAForwardProblem, InvertSSANLCG, InvertSSAIGN, \
tauc_params, LinearPlotListener, PlotListener, pism_print_logger, pism_pause, pauseListener, \
CaptureLogger
from linalg_pism import PISMLocalVector as PLV
import pismssaforward

siple.reporting.clear_loggers()
siple.reporting.add_logger(pism_print_logger)
siple.reporting.set_pause_callback(pism_pause)

default_ssa_l2_coeff = 1.
default_ssa_h1_coeff = 0.
default_rms_error = 100

tauc_guess_scale = 0.2
tauc_guess_const = None


class Vel2Tauc(PISM.ssa.SSAFromBootFile):

  def _constructSSA(self):
    return pismssaforward.SSAForwardSolver(self.grid,self.config)

  def setup(self):
    PISM.ssa.SSAFromBootFile.setup(self)

    # FIXME: This is a lousy name....
    self.solver.init_vars()

  def _initSSACoefficients(self):
    # Read PISM SSA related state variables
    solver = self.solver
    solver.allocateCoeffs(using_l2_weight=True)

    thickness = solver.thickness; bed = solver.bed; enthalpy = solver.enthalpy
    mask = solver.ice_mask; surface = solver.surface

    # Read in the PISM state variables that are used directly in the SSA solver
    for v in [thickness, bed, enthalpy]:
      v.regrid(self.boot_file,True)
  
    # variables mask and surface are computed from the geometry previously read
    sea_level = 0 # FIXME setFromOption?
    gc = PISM.GeometryCalculator(sea_level,self.solver.ice,self.config)
    gc.compute(bed,thickness,mask,surface)

    l2_weight = solver.range_l2_weight
    with PISM.util.Access(comm=l2_weight,nocomm=mask):
      l2_weight.set(0.)
      grounded = PISM.MASK_GROUNDED
      for (i,j) in solver.grid.points():
        # if mask[i-1,j-1]==grounded and mask[i-1,j]==grounded and mask[i-1,j+1]==grounded \
        #    and mask[i,j-1]==grounded and mask[i,j]==grounded and mask[i,j+1]==grounded \
        #    and mask[i+1,j-1]==grounded and mask[i+1,j]==grounded and mask[i+1,j+1]==grounded:
        if mask[i,j] == grounded:
          l2_weight[i,j] = 1

class Vel2TaucPlotListener(PlotListener):
  def __init__(self,grid,Vmax):
    PlotListener.__init__(self,grid)
    self.Vmax = Vmax
    
    self.figure =None
    
  def iteration(self,solver,count,x,Fx,y,d,r,*args):      
    import matplotlib.pyplot as pp

    if self.figure is None:
      self.figure = pp.figure()
    else:
      pp.figure(self.figure.number)

    l2_weight=self.tz_scalar.communicate(solver.forward_problem.solver.range_l2_weight)

    pp.clf()
    
    V = self.Vmax
    pp.subplot(1,3,1)
    rx = l2_weight*r[0,:,:]
    rx = np.maximum(rx,-V)
    rx = np.minimum(rx,V)
    pp.imshow(rx,origin='lower')
    pp.colorbar()
    pp.title('ru')
    pp.jet()

    pp.subplot(1,3,2)
    ry = l2_weight*r[1,:,:]
    ry = np.maximum(ry,-V)
    ry = np.minimum(ry,V)
    pp.imshow(ry,origin='lower')
    pp.colorbar()
    pp.title('rv')
    pp.jet()

    d *= -1
    pp.subplot(1,3,3)      
    pp.imshow(d,origin='lower')
    pp.colorbar()
    pp.jet()
    pp.title('-d')
    
    pp.ion()
    pp.show()

class Vel2TaucLinPlotListener(LinearPlotListener):
  def __init__(self,grid,Vmax):
    LinearPlotListener.__init__(self,grid)
    self.Vmax = Vmax

    self.figure =None

  def iteration(self,solver,count,x,y,d,r,*args):      
    import matplotlib.pyplot as pp

    if self.figure is None:
      self.figure = pp.figure()
    else:
      pp.figure(self.figure.number)

    l2_weight=self.tz_scalar.communicate(solver.forward_problem.solver.range_l2_weight)

    pp.clf()

    V = self.Vmax
    pp.subplot(1,3,1)
    rx = l2_weight*r[0,:,:]
    rx = np.maximum(rx,-V)
    rx = np.minimum(rx,V)
    pp.imshow(rx,origin='lower')
    pp.colorbar()
    pp.title('ru')
    pp.jet()

    pp.subplot(1,3,2)
    ry = l2_weight*r[1,:,:]
    ry = np.maximum(ry,-V)
    ry = np.minimum(ry,V)
    pp.imshow(ry,origin='lower')
    pp.colorbar()
    pp.title('rv')
    pp.jet()

    d *= -1
    pp.subplot(1,3,3)      
    pp.imshow(d,origin='lower')
    pp.colorbar()
    pp.jet()
    pp.title('-d')

    pp.ion()
    pp.show()
  
  


## Main code starts here
if __name__ == "__main__":
  context = PISM.Context()
  config = context.config()
  PISM.set_abort_on_sigint(True)

  usage = \
  """  vel2tauc.py -i IN.nc [-o file.nc]
    where:
      -i      IN.nc is input file in NetCDF format: contains PISM-written model state
    notes:
      * -i is required
    """

  PISM.show_usage_check_req_opts(context.com,"ssa_forward",["-i"],usage)
  
  for o in PISM.OptionsGroup(context.com,"","vel2tauc"):
    bootfile = PISM.optionsString("-i","input file")
    output_file = PISM.optionsString("-o","output file",default="vel2tauc_"+bootfile)
    verbosity = PISM.optionsInt("-verbose","verbosity level",default=2)
    method = PISM.optionsList(context.com,"-inv_method","Inversion algorithm",["nlcg","ign","sd"],"ign")
    rms_error = PISM.optionsReal("-rms_error","RMS velocity error",default=default_rms_error)
    tauc_param_type = PISM.optionsList(context.com,"-tauc_param","zeta->tauc parameterization",["ident","square","exp"],"ident")
    ssa_l2_coeff = PISM.optionsReal("-ssa_l2_coeff","L2 coefficient for domain inner product",default=default_ssa_l2_coeff)
    ssa_h1_coeff = PISM.optionsReal("-ssa_h1_coeff","H1 coefficient for domain inner product",default=default_ssa_h1_coeff)
    tauc_guess_scale = PISM.optionsReal("-tauc_guess_scale","initial guess for tauc to be this factor of the true value",default=tauc_guess_scale)
    tauc_guess_const = PISM.optionsReal("-tauc_guess_const","initial guess for tauc to be this constant",default=tauc_guess_const)
    do_plotting = PISM.optionsFlag("-inv_plot","perform visualization during the computation",default=False)
    do_final_plot = PISM.optionsFlag("-inv_final_plot","perform visualization at the end of the computation",default=False)
    do_pause = PISM.optionsFlag("-inv_pause","pause each iteration",default=False)
    test_adjoint = PISM.optionsFlag("-inv_test_adjoint","Test that the adjoint is working",default=False)
    ls_verbose = PISM.optionsFlag("-inv_ls_verbose","Turn on a verbose linesearch.",default=False)
    ign_theta  = PISM.optionsReal("-ign_theta","theta parameter for IGN algorithm",default=0.5)
    Vmax = PISM.optionsReal("-inv_plot_vmax","maximum velocity for plotting residuals",default=30)
    noise = PISM.optionsReal("-rms_noise","pointwise rms noise to add (in m/a)",default=None)

  config.set_string("inv_ssa_tauc_param",tauc_param_type)
  config.set("inv_ssa_domain_l2_coeff",ssa_l2_coeff)
  config.set("inv_ssa_domain_h1_coeff",ssa_h1_coeff)

  tauc_param = tauc_params[tauc_param_type]

  PISM.setVerbosityLevel(verbosity)
  vel2tauc = Vel2Tauc(bootfile)
  vel2tauc.setup()

  forward_problem = SSAForwardProblem(vel2tauc)

  grid = vel2tauc.grid

  bmr   = PISM.util.standardBasalMeltRateVec(grid)
  tillphi = PISM.util.standardTillPhiVec(grid)
  bwat = PISM.util.standardBasalWaterVec(grid)
  for v in [bmr,tillphi,bwat]:
     v.regrid(bootfile,True)
  pvars = PISM.PISMVars()
  for v in [vel2tauc.solver.thickness,vel2tauc.solver.bed,vel2tauc.solver.ice_mask,bmr,tillphi,bwat]:
     pvars.add(v)

  yieldstress = PISM.PISMDefaultYieldStress(grid,grid.config)
  yieldstress.init(pvars)
  tauc_true = PISM.util.standardYieldStressVec(grid,name="tauc_true")  
  yieldstress.basal_material_yield_stress(tauc_true)

  # Convert tauc_true to zeta_true
  if config.get_string('inv_ssa_tauc_param')=='ident':
    zeta_true = tauc_true
  else:
    zeta_true = PISM.IceModelVec2S();
    zeta_true.create(grid, "zeta_true", True, PISM.util.WIDE_STENCIL)
    with PISM.util.Access(nocomm=[tauc_true],comm=[zeta_true]):
      for (i,j) in grid.points():
        zeta_true[i,j] = tauc_param.fromTauc(tauc_true[i,j])

  # Send the true yeild stress through the forward problem to 
  # get at true velocity field.
  u_true = PISM.util.standard2dVelocityVec(grid,name="_true")
  forward_problem.F(PLV(zeta_true),out=PLV(u_true))

  if not noise is None:
    u_noise = pismssaforward.randVectorV(grid,noise/math.sqrt(2))
    unv = PLV(u_noise)
    unv += PLV(u_true)
  else:
    u_noise = u_true
  
  # Build the initial guess for tauc for the inversion.
  tauc = PISM.util.standardYieldStressVec(grid)
  if not tauc_guess_const is None:
    tauc.set(tauc_guess_const)
  else:
    tauc.copy_from(tauc_true)
    tauc.scale(tauc_guess_scale)

  # Convert tauc guess to zeta guess
  if config.get_string('inv_ssa_tauc_param')=='ident':
    zeta = tauc
  else:
    zeta = PISM.IceModelVec2S();
    zeta.create(grid, "zeta", True, PISM.util.WIDE_STENCIL)
    with PISM.util.Access(nocomm=tauc,comm=zeta):
      for (i,j) in grid.points():
        zeta[i,j] = tauc_param.fromTauc(tauc[i,j])

  if test_adjoint:
    d = PLV(pismssaforward.randVectorS(grid,1e5))
    r = PLV(pismssaforward.randVectorV(grid,1./PISM.secpera))
    (domainIP,rangeIP)=forward_problem.testTStar(PLV(zeta),d,r,3)
    siple.reporting.msg("domainip %g rangeip %g",domainIP,rangeIP)
    exit(0)

  # Determine the inversion algorithm, and set up its arguments.
  if method == "ign":
    Solver = InvertSSAIGN
  else:
    Solver = InvertSSANLCG

  # u_true = PISM.util.standard2dVelocityVec(grid,"bar_ssa")
  # u_true.regrid(bootfile,True)

  params=Solver.defaultParameters()
  if method == "sd":
    params.steepest_descent = True
    params.ITER_MAX=10000
  elif method =="ign":
    params.linearsolver.ITER_MAX=10000
    params.linearsolver.verbose = True
    params.thetaMax=ign_theta
  if ls_verbose:
    params.linesearch.verbose = True
  params.verbose   = True
  params.deriv_eps = 0.

  # Run the inversion
  logger = CaptureLogger();
  solver=Solver(forward_problem,params=params)  
  if do_plotting:
    solver.addIterationListener(Vel2TaucPlotListener(grid,Vmax))
    if method=='ign':
      solver.addLinearIterationListener(Vel2TaucLinPlotListener(grid,Vmax))
  if do_pause:
    solver.addIterationListener(pauseListener)

  rms_error /= PISM.secpera # m/s
  (zeta,u) = solver.solve(zeta,u_true,rms_error)

  # u = u_true

  # Convert back from zeta to tauc
  if config.get_string('inv_ssa_tauc_param')=='ident':
    tauc = zeta
  else:
    with PISM.util.Access(nocomm=zeta,comm=tauc):
      for (i,j) in grid.points():
        (tauc[i,j],dummy) = tauc_param.toTauc(zeta[i,j])

  # Write solution out to netcdf file
  vel2tauc.write(output_file)
  tauc.write(output_file)
  tauc_true.write(output_file)
  u.set_name("_computed",0)
  u.write(output_file)
  logger.write(output_file)
  
  # Draw a pretty picture
  tz = tozero.ToProcZero(grid)
  tauc_a = tz.communicate(tauc)
  tauc_true = tz.communicate(tauc_true)
  if do_final_plot and (not tauc_a is None):
    from matplotlib import pyplot
    pyplot.clf()
    pyplot.subplot(1,2,1)
    pyplot.imshow(tauc_a,origin='lower')
    pyplot.subplot(1,2,2)
    pyplot.imshow(tauc_true,origin='lower')
    pyplot.ion()
    pyplot.show()
    siple.reporting.endpause()
