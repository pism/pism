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
import os, math

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



class Vel2Tauc(PISM.ssa.SSAFromBootFile):
  def __init__(self,input_filename,inv_data_filename):
    PISM.ssa.SSAFromBootFile.__init__(self,input_filename)
    self.inv_data_filename = inv_data_filename

  def _constructSSA(self):
    return pismssaforward.SSAForwardSolver(self.grid,self.config)

  def setup(self):
    PISM.ssa.SSAFromBootFile.setup(self)

    # FIXME: This is a lousy name....
    self.solver.init_vars()

  def _initSSACoefficients(self):
    # Read PISM SSA related state variables
    solver = self.solver
    solver.allocateCoeffs()

    thickness = solver.thickness; bed = solver.bed; enthalpy = solver.enthalpy
    mask = solver.ice_mask; surface = solver.surface

    # Read in the PISM state variables that are used directly in the SSA solver
    for v in [thickness, bed, enthalpy]:
      v.regrid(self.boot_file,True)
  
    # variables mask and surface are computed from the geometry previously read
    sea_level = 0 # FIXME setFromOption?
    gc = PISM.GeometryCalculator(sea_level,self.solver.ice,self.config)
    gc.compute(bed,thickness,mask,surface)

    weight = solver.vel_misfit_weight
    weight.regrid(self.inv_data_filename,True)


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

    l2_weight=self.tz_scalar.communicate(solver.forward_problem.solver.vel_misfit_weight)

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

    l2_weight=self.tz_scalar.communicate(solver.forward_problem.solver.vel_misfit_weight)

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
    input_filename = PISM.optionsString("-i","input file")
    output_filename = PISM.optionsString("-o","output file",default="vel2tauc_"+os.path.basename(input_filename))
    inv_data_filename = PISM.optionsString("-inv_data","inverse data file",default=input_filename)
    verbosity = PISM.optionsInt("-verbose","verbosity level",default=2)
    method = PISM.optionsList(context.com,"-inv_method","Inversion algorithm",["nlcg","ign","sd"],"ign")
    rms_error = PISM.optionsReal("-rms_error","RMS velocity error",default=default_rms_error)
    tauc_param_type = PISM.optionsList(context.com,"-tauc_param","zeta->tauc parameterization",["ident","square","exp"],"ident")
    ssa_l2_coeff = PISM.optionsReal("-ssa_l2_coeff","L2 coefficient for domain inner product",default=default_ssa_l2_coeff)
    ssa_h1_coeff = PISM.optionsReal("-ssa_h1_coeff","H1 coefficient for domain inner product",default=default_ssa_h1_coeff)
    do_plotting = PISM.optionsFlag("-inv_plot","perform visualization during the computation",default=False)
    do_final_plot = PISM.optionsFlag("-inv_final_plot","perform visualization at the end of the computation",default=False)
    do_pause = PISM.optionsFlag("-inv_pause","pause each iteration",default=False)
    test_adjoint = PISM.optionsFlag("-inv_test_adjoint","Test that the adjoint is working",default=False)
    ls_verbose = PISM.optionsFlag("-inv_ls_verbose","Turn on a verbose linesearch.",default=False)
    ign_theta  = PISM.optionsReal("-ign_theta","theta parameter for IGN algorithm",default=0.5)
    Vmax = PISM.optionsReal("-inv_plot_vmax","maximum velocity for plotting residuals",default=30)

  config.set_string("inv_ssa_tauc_param",tauc_param_type)
  config.set("inv_ssa_domain_l2_coeff",ssa_l2_coeff)
  config.set("inv_ssa_domain_h1_coeff",ssa_h1_coeff)

  tauc_param = tauc_params[tauc_param_type]

  PISM.setVerbosityLevel(verbosity)
  vel2tauc = Vel2Tauc(input_filename,inv_data_filename)
  vel2tauc.setup()

  forward_problem = SSAForwardProblem(vel2tauc)

  grid = vel2tauc.grid

  modeldata = vel2tauc.solver
  grid = modeldata.grid

  vel_ssa_observed = PISM.util.standard2dVelocityVec(grid,'_ssa_observed')
  vel_ssa_observed.regrid(inv_data_filename,True)
  
  tauc_prior = PISM.util.standardYieldStressVec(grid,'tauc_prior')
  tauc_prior.regrid(inv_data_filename,True)

  tauc = PISM.util.standardYieldStressVec(grid)
  tauc.copy_from(tauc_prior)

  # Convert tauc guess to zeta guess
  if config.get_string('inv_ssa_tauc_param')=='ident':
    zeta = tauc_prior
  else:
    zeta = PISM.IceModelVec2S();
    zeta.create(grid, "zeta", True, PISM.util.WIDE_STENCIL)
    with PISM.util.Access(nocomm=tauc_prior,comm=zeta):
      for (i,j) in grid.points():
        zeta[i,j] = tauc_param.fromTauc(tauc_prior[i,j])

  if test_adjoint:
    d = PLV(pismssaforward.randVectorS(grid,1e5,PISM.util.WIDE_STENCIL))
    r = PLV(pismssaforward.randVectorV(grid,1./PISM.secpera,PISM.util.WIDE_STENCIL))
    (domainIP,rangeIP)=forward_problem.testTStar(PLV(zeta),d,r,3)
    siple.reporting.msg("domainip %g rangeip %g",domainIP,rangeIP)
    exit(0)

  # Determine the inversion algorithm, and set up its arguments.
  if method == "ign":
    Solver = InvertSSAIGN
  else:
    Solver = InvertSSANLCG

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
  (zeta,u) = solver.solve(zeta,vel_ssa_observed,rms_error)

  # Convert back from zeta to tauc
  if config.get_string('inv_ssa_tauc_param')=='ident':
    tauc = zeta
  else:
    with PISM.util.Access(nocomm=zeta,comm=tauc):
      for (i,j) in grid.points():
        (tauc[i,j],dummy) = tauc_param.toTauc(zeta[i,j])

  # Write solution out to netcdf file
  vel2tauc.write(output_filename)
  tauc.write(output_filename)
  tauc_prior.write(output_filename)
  
  u.set_name("_ssa_inv",0)
  u.write(output_filename)

  logger.write(output_filename)
