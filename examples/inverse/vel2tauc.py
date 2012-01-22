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

from pismssaforward import SSAForwardProblem, InvertSSANLCG, InvertSSAIGN, \
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

  def setup(self):

    PISM.ssa.SSAFromBootFile.setup(self)

    vecs = self.modeldata.vecs

    # The SSA instance will not keep a reference to pismVars; it only uses it to extract
    # its desired variables.  So it is safe to pass it pismVars and then let pismVars
    # go out of scope at the end of this method.

    self.ssa.init(vecs.asPISMVars())

    if vecs.has('vel_bc'):
      self.ssa.set_boundary_conditions(vecs.bc_mask,vecs.vel_bc)

    # FIXME: Fix this lousy name
    self.ssa.setup_vars()

    def _constructSSA(self):
      md = self.modeldata
      tauc_param_type = self.config.get_string("inv_ssa_tauc_param")
      self.tauc_param = tauc_params[tauc_param_type]
      return PISM.InvSSAForwardProblem(md.grid,md.basal,md.ice,md.enthalpyconverter,self.tauc_param,self.config)

  def _initSSACoefficients(self):
    self._allocStdSSACoefficients()
    
    # Read PISM SSA related state variables

    vecs = self.modeldata.vecs
    thickness = vecs.thickness; bed = vecs.bed; enthalpy = vecs.enthalpy
    mask = vecs.ice_mask; surface = vecs.surface

    # Read in the PISM state variables that are used directly in the SSA solver
    for v in [thickness, bed, enthalpy]:
      v.regrid(self.boot_file,True)
  
    # variables mask and surface are computed from the geometry previously read
    sea_level = 0 # FIXME setFromOption?
    gc = PISM.GeometryCalculator(sea_level,self.modeldata.ice,self.config)
    gc.compute(bed,thickness,mask,surface)

    vecs.add( PISM.util.standardVelocityMisfitWeight(self.grid) )
    weight = vecs.vel_misfit_weight
    weight.regrid(self.inv_data_filename,True)


  def _constructSSA(self):
    md = self.modeldata
    vecs  = self.modeldata.vecs
    tauc_param_type = self.config.get_string("inv_ssa_tauc_param")
    self.tauc_param = tauc_params[tauc_param_type]
    return PISM.InvSSAForwardProblem(md.grid,md.basal,md.ice,md.enthalpyconverter,self.tauc_param,self.config)


  def write(self,filename,append=False):
    if not append:
      PISM.ssa.SSAFromBootFile.write(self,filename)
    else:
      grid = self.grid
      vecs = self.modeldata.vecs

      pio = PISM.PISMIO(grid)
      pio.open_for_writing(filename,True,True) #append mode!
      
      self.modeldata.vecs.write(filename)
      pio.close()

      # Save time & command line
      PISM.util.writeProvenance(filename)


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
  config = context.config
  PISM.set_abort_on_sigint(True)

  usage = \
  """  vel2tauc.py -i IN.nc [-o file.nc]
    where:
      -i      IN.nc is input file in NetCDF format: contains PISM-written model state
    notes:
      * -i is required
    """

  # FIXME:  Required should be -i or -a
  # PISM.show_usage_check_req_opts(context.com,"ssa_forward",["-i"],usage)
  
  append_mode = False
  for o in PISM.OptionsGroup(context.com,"","vel2tauc"):
    input_filename = PISM.optionsString("-i","input file")
    append_filename = PISM.optionsString("-a","append file",default=None)
    output_filename = PISM.optionsString("-o","output file",default=None)

    if (not input_filename is None) and (not append_filename is None):
      raise RuntimeError("Only one of -i/-a is allowed.")

    if (not output_filename is None) and (not append_filename is None):
      raise RuntimeError("Only one of -a/-o is allowed.")

    if not(append_filename is None):
      input_filename = append_filename
      output_filename = append_filename
      append_mode = True

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
    use_tauc_prior = PISM.optionsFlag("-use_tauc_prior","Use tauc_prior from inverse data file as initial guess.",default=False)
    ign_theta  = PISM.optionsReal("-ign_theta","theta parameter for IGN algorithm",default=0.5)
    Vmax = PISM.optionsReal("-inv_plot_vmax","maximum velocity for plotting residuals",default=30)

  if output_filename is None:
    output_filename = "vel2tauc_"+os.path.basename(input_filename)    

  saving_inv_data = (inv_data_filename != output_filename)

  config.set_string("inv_ssa_tauc_param",tauc_param_type)
  config.set("inv_ssa_domain_l2_coeff",ssa_l2_coeff)
  config.set("inv_ssa_domain_h1_coeff",ssa_h1_coeff)

  tauc_param = tauc_params[tauc_param_type]

  PISM.setVerbosityLevel(verbosity)
  vel2tauc = Vel2Tauc(input_filename,inv_data_filename)
  vel2tauc.setup()

  forward_problem = SSAForwardProblem(vel2tauc)

  grid = vel2tauc.grid

  modeldata = vel2tauc.modeldata
  vecs = modeldata.vecs
  grid = modeldata.grid

  vel_ssa_observed = None
  vel_ssa_observed = PISM.util.standard2dVelocityVec(grid,'_ssa_observed',stencil_width=2)
  if PISM.util.fileHasVariable(inv_data_filename,"u_ssa_observed"):
    vel_ssa_observed.regrid(inv_data_filename,True)
    vecs.add(vel_ssa_observed,writing=saving_inv_data)
  else:
    if not PISM.util.fileHasVariable(inv_data_filename,"u_surface_observed"):
      PISM.verbPrintf(1,context.com,"Neither u/v_ssa_observed nor u/v_surface_observed is available in %s.\nAt least one must be specified.\n" % inv_data_filename)
      exit()
    vel_surface_observed = PISM.util.standard2dVelocityVec(grid,'_surface_observed',stencil_width=2)
    vel_surface_observed.regrid(inv_data_filename,True)
    vecs.add(vel_surface_observed,writing=saving_inv_data)
    
    vel_sia_observed = pismssaforward.computeSIASurfaceVelocities(modeldata)
    vel_sia_observed.rename('_sia_observed',"'observed' SIA velocities'","")
    vel_ssa_observed.copy_from(vel_surface_observed)
    vel_ssa_observed.add(-1,vel_sia_observed)
    vecs.add(vel_ssa_observed,writing=True)

  # Determine the prior guess for tauc: either the tauc from the input
  # file, or if -using-tauc-prior the tauc_prior from the inv_datafile
  tauc_prior = PISM.util.standardYieldStressVec(grid,'tauc_prior')
  tauc = PISM.util.standardYieldStressVec(grid)
  if use_tauc_prior:
    tauc_prior.regrid(inv_data_filename,True)
  else:
    if not PISM.util.fileHasVariable(input_filename,"tauc"):
      verbPrintf(1,com,"Initial guess for tauc is not available as 'tauc' in %s.\nYou can provide an initial guess as 'tauc_prior' using the command line option -use_tauc_prior." % input_filename)
      exit()
    tauc.regrid(inv_data_filename,True)
    tauc_prior.copy_from(tauc)
  vecs.add(tauc_prior,writing=saving_inv_data)


  # Convert tauc guess to zeta guess
  if config.get_string('inv_ssa_tauc_param')=='ident':
    tauc.copy_from(tauc_prior)
    zeta = tauc
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

  # It may be that a 'tauc' was read in earlier.  We replace it with
  # our newly generated one.
  if vecs.has('tauc'): vecs.remove('tauc')
  vecs.add(tauc,writing=True)

  u.rename("_ssa_inv","SSA velocity computed by inversion","")
  vecs.add(u,writing=True)

  # Write solution out to netcdf file
  vel2tauc.write(output_filename,append=append_mode)

  logger.write(output_filename)
