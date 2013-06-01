#! /usr/bin/env python
#
# Copyright (C) 2011, 2012, 2013 David Maxwell
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
import PISM.invert.ssa
import numpy as np
import sys, os, math

class Vel2Tauc(PISM.invert.ssa.SSAForwardRunFromInputFile):

  def write(self,filename,append=False):
    if not append:
      PISM.invert.ssa.SSAForwardRunFromInputFile.write(self,filename)
    else:
      grid = self.grid
      vecs = self.modeldata.vecs

      pio = PISM.PIO(grid.com,grid.rank,"netcdf3", grid.get_unit_system())
      pio.open(filename,PISM.NC_WRITE,True) #append mode!

      self.modeldata.vecs.write(filename)
      pio.close()

class Vel2TaucPlotListener(PISM.invert.listener.PlotListener):
  def __init__(self,grid,Vmax):
    PISM.invert.listener.PlotListener.__init__(self,grid)
    self.Vmax = Vmax
    self.l2_weight = None
    self.l2_weight_init = False

  def __call__(self,inverse_solver,count,data):

    if self.l2_weight_init == False:
      vecs = inverse_solver.ssarun.modeldata.vecs;
      self.l2_weight=self.toproczero(vecs.vel_misfit_weight)
      self.l2_weight_init = True

    method = inverse_solver.method

    r=self.toproczero(data.r)
    Td = None
    if data.has_key('Td'): Td = self.toproczero(data.Td)
    TStarR = None
    if data.has_key('TStarR'): TStarR = self.toproczero(data.TStarR)
    d = None
    if data.has_key('d'): d = self.toproczero(data.d)
    zeta = self.toproczero(data.zeta)      

    secpera = grid.convert(1.0, "year", "second")

    if self.grid.rank == 0:
      import matplotlib.pyplot as pp
    
      pp.figure(self.figure())

      l2_weight=self.l2_weight

      pp.clf()
    
      V = self.Vmax

      pp.subplot(2,3,1)
      rx = l2_weight*r[0,:,:]*secpera
      rx = np.maximum(rx,-V)
      rx = np.minimum(rx,V)
      pp.imshow(rx,origin='lower',interpolation='nearest')
      pp.colorbar()
      pp.title('r_x')
      pp.jet()

      pp.subplot(2,3,4)
      ry = l2_weight*r[1,:,:]*secpera
      ry = np.maximum(ry,-V)
      ry = np.minimum(ry,V)
      pp.imshow(ry,origin='lower',interpolation='nearest')
      pp.colorbar()
      pp.title('r_y')
      pp.jet()
    
    
      if method == 'ign':
        pp.subplot(2,3,2)
        Tdx = Td[0,:,:]*secpera
        pp.imshow(Tdx,origin='lower',interpolation='nearest')
        pp.colorbar()
        pp.title('Td_x')
        pp.jet()

        pp.subplot(2,3,5)
        Tdy = Td[1,:,:]*secpera
        pp.imshow(Tdy,origin='lower',interpolation='nearest')
        pp.colorbar()
        pp.title('Td_y')
        pp.jet()
      elif method == 'sd' or method == 'nlcg':
        pp.subplot(2,3,2)
        pp.imshow(TStarR,origin='lower',interpolation='nearest')
        pp.colorbar()
        pp.title('TStarR')
        pp.jet()

      if d is not None:
        d *= -1
        pp.subplot(2,3,3)      
        pp.imshow(d,origin='lower',interpolation='nearest')
        pp.colorbar()
        pp.jet()
        pp.title('-d')

      pp.subplot(2,3,6)
      pp.imshow(zeta,origin='lower',interpolation='nearest')
      pp.colorbar()
      pp.jet()
      pp.title('zeta')
    
      pp.ion()
      pp.show()

class Vel2TaucLinPlotListener(PISM.invert.listener.PlotListener):
  def __init__(self,grid,Vmax):
    PISM.invert.listener.PlotListener.__init__(self,grid)
    self.Vmax = Vmax
    self.l2_weight = None
    self.l2_weight_init = False

  def __call__(self,inverse_solver,count,data):
    # On the first go-around, extract the l2_weight vector onto 
    # processor zero.
    if self.l2_weight_init == False:
      vecs = inverse_solver.ssarun.modeldata.vecs;
      self.l2_weight = self.toproczero(vecs.vel_misfit_weight)
      self.l2_init = True

    l2_weight=self.l2_weight
    r = self.toproczero(data.r)
    d = self.toproczero(data.d)

    if self.grid.rank == 0:
      import matplotlib.pyplot as pp
      pp.figure(self.figure())
      pp.clf()
    
      V = self.Vmax
      pp.subplot(1,3,1)
      rx = l2_weight*r[0,:,:]
      rx = np.maximum(rx,-V)
      rx = np.minimum(rx,V)
      pp.imshow(rx,origin='lower',interpolation='nearest')
      pp.colorbar()
      pp.title('ru')
      pp.jet()
    
      pp.subplot(1,3,2)
      ry = l2_weight*r[1,:,:]
      ry = np.maximum(ry,-V)
      ry = np.minimum(ry,V)
      pp.imshow(ry,origin='lower',interpolation='nearest')
      pp.colorbar()
      pp.title('rv')
      pp.jet()
    
      d *= -1
      pp.subplot(1,3,3)      
      pp.imshow(d,origin='lower',interpolation='nearest')
      pp.colorbar()
      pp.jet()
      pp.title('-d')
    
      pp.ion()
      pp.show()

def adjustTauc(mask,tauc):
  """Where ice is floating or land is ice-free, tauc should be adjusted to have some preset default values."""

  grid = mask.get_grid()
  high_tauc = grid.config.get("high_tauc")

  with PISM.vec.Access(comm=tauc,nocomm=mask):
    mq = PISM.MaskQuery(mask)
    for (i,j) in grid.points():
      if mq.ocean(i,j):
        tauc[i,j] = 0;
      elif mq.ice_free(i,j):
        tauc[i,j] = high_tauc

## Main code starts here
if __name__ == "__main__":
  context = PISM.Context()
  config = context.config
  com = context.com
  PISM.set_abort_on_sigint(True)

  WIDE_STENCIL = 2

  usage = \
  """  vel2tauc.py -i IN.nc [-o file.nc]
    where:
      -i      IN.nc is input file in NetCDF format: contains PISM-written model state
    notes:
      * -i is required
    """

  append_mode = False
  PISM.setVerbosityLevel(1)
  for o in PISM.OptionsGroup(context.com,"","vel2tauc"):
    input_filename = PISM.optionsString("-i","input file")
    append_filename = PISM.optionsString("-a","append file",default=None)
    output_filename = PISM.optionsString("-o","output file",default=None)

    if (input_filename is None) and (append_filename is None):
      PISM.verbPrintf(1,com,"\nError: No input file specified. Use one of -i [file.nc] or -a [file.nc].\n")
      PISM.PISMEndQuiet()

    if (input_filename is not None) and (append_filename is not None):
      PISM.verbPrintf(1,com,"\nError: Only one of -i/-a is allowed.\n")
      PISM.PISMEndQuiet()

    if (output_filename is not None) and (append_filename is not None):
      PISM.verbPrintf(1,com,"\nError: Only one of -a/-o is allowed.\n")
      PISM.PISMEndQuiet()

    if append_filename is not None:
      input_filename = append_filename
      output_filename = append_filename
      append_mode = True

    inv_data_filename = PISM.optionsString("-inv_data","inverse data file",default=input_filename)
    verbosity = PISM.optionsInt("-verbose","verbosity level",default=2)

    do_plotting = PISM.optionsFlag("-inv_plot","perform visualization during the computation",default=False)
    do_final_plot = PISM.optionsFlag("-inv_final_plot","perform visualization at the end of the computation",default=False)
    do_pause = PISM.optionsFlag("-inv_pause","pause each iteration",default=False)
    test_adjoint = PISM.optionsFlag("-inv_test_adjoint","Test that the adjoint is working",default=False)
    ls_verbose = PISM.optionsFlag("-inv_ls_verbose","Turn on a verbose linesearch.",default=False)
    do_restart = PISM.optionsFlag("-inv_restart","Restart a stopped computation.",default=False)
    use_tauc_prior = PISM.optionsFlag("-inv_use_tauc_prior","Use tauc_prior from inverse data file as initial guess.",default=False)
    ign_theta  = PISM.optionsReal("-ign_theta","theta parameter for IGN algorithm",default=0.5)
    Vmax = PISM.optionsReal("-inv_plot_vmax","maximum velocity for plotting residuals",default=30)

    prep_module = PISM.optionsString("-inv_prep_module","Python module used to do final setup of inverse solver",default=None)

    is_regional = PISM.optionsFlag("-regional","Compute SIA/SSA using regional model semantics",default=False)


  inv_method = config.get_string("inv_ssa_method")
  
  
  if output_filename is None:
    output_filename = "vel2tauc_"+os.path.basename(input_filename)    

  saving_inv_data = (inv_data_filename != output_filename)

  PISM.setVerbosityLevel(verbosity)
  vel2tauc = Vel2Tauc(input_filename,inv_data_filename,'tauc')
  vel2tauc.setup()
  tauc_param = vel2tauc.design_var_param
  solver = PISM.invert.ssa.createInvSSASolver(vel2tauc)

  # if forward_type == 'classic':
  #   forward_problem = SSAForwardProblem(vel2tauc)
  # else:
  #   forward_problem = SSAForwardProblemFIXME(vel2tauc)

  modeldata = vel2tauc.modeldata
  vecs = modeldata.vecs
  grid = modeldata.grid

  # Determine the prior guess for tauc. This can be one of 
  # a) tauc from the input file (default)
  # b) tauc_prior from the inv_datafile if -inv_use_tauc_prior is set
  tauc_prior = PISM.model.createYieldStressVec(grid,'tauc_prior')
  tauc_prior.set_attrs("diagnostic", "initial guess for (pseudo-plastic) basal yield stress in an inversion", "Pa", "");
  tauc = PISM.model.createYieldStressVec(grid)
  if use_tauc_prior:
    tauc_prior.regrid(inv_data_filename,critical=True)
    vecs.add(tauc_prior,writing=saving_inv_data)
  else:
    if not PISM.util.fileHasVariable(input_filename,"tauc"):
      PISM.verbPrintf(1,com,"Initial guess for tauc is not available as 'tauc' in %s.\nYou can provide an initial guess as 'tauc_prior' using the command line option -use_tauc_prior." % input_filename)
      exit(1)
    tauc.regrid(input_filename,True)
    tauc_prior.copy_from(tauc)
    vecs.add(tauc_prior,writing=True)

  adjustTauc(vecs.ice_mask,tauc_prior)

  # Convert tauc_prior -> zeta_prior
  zeta_prior = PISM.IceModelVec2S();
  zeta_prior.create(grid, "zeta_prior", PISM.kHasGhosts, WIDE_STENCIL)
  tauc_param.convertFromDesignVariable(tauc_prior,zeta_prior)
  vecs.add(zeta_prior,writing=True)

  # If the inverse data file has a variable tauc_true, this is probably
  # a synthetic inversion.  We'll load it now so that it will get written
  # out, if needed, at the end of the computation in the output file.
  if PISM.util.fileHasVariable(inv_data_filename,"tauc_true"):
    tauc_true = PISM.model.createYieldStressVec(grid,'tauc_true')
    tauc_true.regrid(inv_data_filename,True)
    tauc_true.read_attributes(inv_data_filename)
    vecs.add(tauc_true,writing=saving_inv_data)

  # Determine the initial guess for zeta.  If we are not
  # restarting, we convert tauc_prior to zeta.  If we are restarting,
  # we load in zeta from the output file.
  zeta = PISM.IceModelVec2S();
  zeta.create(grid, "zeta_inv", PISM.kHasGhosts, WIDE_STENCIL)
  if do_restart:
    # Just to be sure, verify that we have a 'zeta_inv' in the output file.
    if not PISM.util.fileHasVariable(output_filename,'zeta_inv'):
      PISM.verbPrintf(1,com,"Unable to restart computation: file %s is missing variable 'zeta_inv'", output_filename)
      exit(1)
    zeta.regrid(output_filename,True)
  else:
    zeta.copy_from(zeta_prior)

  if test_adjoint:
    if solver.method.startswith('tikhonov') and solver.method != 'tikhonov_gn':
      PISM.logging.logMessage("option -inv_test_adjoint cannot be used with inverse method %s",solver.method)
      exit(1)

    if solver.method == 'tikhonov_gn':
      pass
    else:
      import numpy as np
      (seed,seed_set) = PISM.optionsIntWasSet("-inv_test_adjoint_seed","")
      if seed_set:
        np.random.seed(seed+PISM.Context().rank)
      d = PISM.vec.randVectorS(grid,1e5,WIDE_STENCIL)
      # If we're fixing some tauc values, we need to ensure that we don't
      # move in a direction 'd' that changes those values in this test.
      if vel2tauc.using_zeta_fixed_mask:
        zeta_fixed_mask = vecs.zeta_fixed_mask
        with PISM.vec.Access(comm=d, nocomm=zeta_fixed_mask):
          for (i,j) in grid.points():
            if zeta_fixed_mask[i,j] != 0:
              d[i,j] = 0;
      r = PISM.vec.randVectorV(grid,1./secpera,WIDE_STENCIL)
      from PISM.invert.sipletools import PISMLocalVector as PLV
      forward_problem = solver.forward_problem
      (domainIP,rangeIP)=forward_problem.testTStar(PLV(zeta),PLV(d),PLV(r),3)
      PISM.logging.logMessage("domainip %.10g rangeip %.10g\n" % (domainIP,rangeIP) )
      PISM.logging.logMessage("relative error %.10g\n" % abs((domainIP-rangeIP)/domainIP))
      exit(0)

  vel_ssa_observed = None
  vel_ssa_observed = PISM.model.create2dVelocityVec(grid,'_ssa_observed',stencil_width=2)
  if PISM.util.fileHasVariable(inv_data_filename,"u_ssa_observed"):
    vel_ssa_observed.regrid(inv_data_filename,True)
    vecs.add(vel_ssa_observed,writing=saving_inv_data)
  else:
    if not PISM.util.fileHasVariable(inv_data_filename,"u_surface_observed"):
      PISM.verbPrintf(1,context.com,"Neither u/v_ssa_observed nor u/v_surface_observed is available in %s.\nAt least one must be specified.\n" % inv_data_filename)
      exit(1)
    vel_surface_observed = PISM.model.create2dVelocityVec(grid,'_surface_observed',stencil_width=2)
    vel_surface_observed.regrid(inv_data_filename,True)
    vecs.add(vel_surface_observed,writing=saving_inv_data)
    
    sia_solver=PISM.SIAFD
    if is_regional:
      sia_solver=PISM.SIAFD_Regional
    vel_sia_observed = PISM.sia.computeSIASurfaceVelocities(modeldata,sia_solver)
    vel_sia_observed.rename('_sia_observed',"'observed' SIA velocities'","")
    vel_ssa_observed.copy_from(vel_surface_observed)
    vel_ssa_observed.add(-1,vel_sia_observed)
    vecs.add(vel_ssa_observed,writing=True)

  # Establish a logger which will save logging messages to the output file.  
  logger = PISM.logging.CaptureLogger(output_filename,'vel2tauc_log');
  PISM.logging.add_logger(logger)
  if append_mode or do_restart:
    logger.readOldLog()
  
  # Prep the output file from the grid so that we can save zeta to it during the runs.
  if not append_mode:
    pio = PISM.PIO(grid.com,grid.rank,"netcdf3", grid.get_unit_system())
    pio.open(output_filename,PISM.NC_WRITE,False)
    pio.def_time(grid.config.get_string("time_dimension_name"),
                 grid.config.get_string("calendar"), grid.time.units_string())
    pio.append_time(grid.config.get_string("time_dimension_name"),grid.time.current())
    pio.close()
  zeta.write(output_filename)


  # Log the command line to the output file now so that we have a record of
  # what was attempted
  PISM.util.writeProvenance(output_filename)    

  # Attach various iteration listeners to the solver as needed for:
  # Plotting
  if do_plotting:
    solver.addIterationListener(Vel2TaucPlotListener(grid,Vmax))
    if solver.method=='ign':
      solver.addLinearIterationListener(Vel2TaucLinPlotListener(grid,Vmax))
  # Pausing
  if do_pause:
    solver.addIterationListener(PISM.invert.listener.pauseListener)
  # Progress reporting
  progress_reporter = None
  if inv_method.startswith('tik'):
    progress_reporter = PISM.invert.ssa.PrintTikhonovProgress()
  else:
    progress_reporter = PISM.invert.ssa.PrintRMSMisfit()
  if progress_reporter is not None:
    solver.addIterationListener(progress_reporter)

  # Saving the current iteration
  solver.addDesignUpdateListener(PISM.invert.ssa.ZetaSaver(output_filename)) 

  # Solver is set up.  Give the user's prep module a chance to do any final
  # setup.
  
  if prep_module is not None:
    exec "import %s as user_prep_module" % prep_module
    user_prep_module.prep_solver(solver)

  # Run the inverse solver!
  if do_restart:
    PISM.logging.logMessage('************** Restarting inversion. ****************\n')
  else:
    PISM.logging.logMessage('============== Starting inversion. ==================\n')  

  # Try solving
  reason = solver.solveInverse(zeta_prior,vel_ssa_observed,zeta);
  if reason.failed():
    PISM.logging.logError("Inverse solve FAILURE:\n%s\n" % reason.nested_description(1));
    quit()

  PISM.logging.logMessage("Inverse solve success (%s)!\n" % reason.description());

  (zeta,u) = solver.inverseSolution()

  # Convert back from zeta to tauc
  tauc_param.convertToDesignVariable(zeta,tauc)

  # It may be that a 'tauc' was read in earlier.  We replace it with
  # our newly generated one.
  if vecs.has('tauc'): vecs.remove('tauc')
  vecs.add(tauc,writing=True)

  vecs.add(zeta,writing=True)

  u.rename("_ssa_inv","SSA velocity computed by inversion","")
  vecs.add(u,writing=True)

  residual = PISM.model.create2dVelocityVec(grid,name='_inv_ssa_residual')
  residual.copy_from(u)
  residual.add(-1,vel_ssa_observed);
  
  r_mag = PISM.IceModelVec2S();
  r_mag.create(grid,"inv_ssa_residual", PISM.kNoGhosts,0);
  
  r_mag.set_attrs("diagnostic","magnitude of mismatch between observed surface velocities and their reconstrution by inversion",
            "m s-1", "inv_ssa_residual", 0);
  r_mag.set_attr("_FillValue", grid.convert(-0.01,'m/year','m/s'));
  r_mag.set_attr("valid_min", 0.0);
  r_mag.set_glaciological_units("m year-1")
  r_mag.write_in_glaciological_units = True

  residual.magnitude(r_mag)
  r_mag.mask_by(vecs.thickness)
  
  vecs.add(residual,writing=True)
  vecs.add(r_mag,writing=True)

  # Write solution out to netcdf file
  vel2tauc.write(output_filename,append=append_mode)
  # If we're not in append mode, the previous command just nuked
  # the output file.  So we rewrite the siple log.
  if not append_mode:
    logger.write(output_filename)

  # Save the misfit history
  if progress_reporter is not None:
    progress_reporter.write(output_filename)
