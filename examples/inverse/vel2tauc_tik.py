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

import sys, petsc4py
petsc4py.init(sys.argv)
from petsc4py import PETSc
import numpy as np
import os, math

import PISM
from PISM.invert_ssa import SSAForwardProblem, InvertSSANLCG, InvertSSAIGN, \
tauc_param_factory, LinearPlotListener, PlotListener 
from PISM.sipletools import pism_print_logger, pism_pause, pauseListener, \
CaptureLogger, CarefulCaptureLogger
from PISM.sipletools import PISMLocalVector as PLV

import siple

siple.reporting.clear_loggers()
siple.reporting.add_logger(pism_print_logger)
siple.reporting.set_pause_callback(pism_pause)

default_ssa_l2_coeff = 1.
default_ssa_h1_coeff = 0.
default_eta = 1

class Listener(PISM.PythonTikhonovSVListener):
  def __init__(self):
    PISM.PythonTikhonovSVListener.__init__(self)
  def iteration(self,it,eta,objVal,penaltyVal,d,diff_d,grad_d,u,diff_u,grad_u,grad):
    print "----------------------------------------------------------"
    print "Iteration %d" % it
    print "RMS misfit: %g" % math.sqrt(penaltyVal)
    print "sqrt(design objective) %g; weighted %g" % (math.sqrt(objVal),math.sqrt(objVal/eta)) 
    print "gradient: design %g state %g sum %g " % (grad_d.norm(PETSc.NormType.NORM_2)/eta,grad_u.norm(PETSc.NormType.NORM_2),grad.norm(PETSc.NormType.NORM_2))
    print "tikhonov functional: %g" % (objVal/eta + penaltyVal)
#    print "doing iteration %d: gradients %g %g %g" % (it,grad_d.norm(PETSc.NormType.NORM_2),eta*grad_u.norm(PETSc.NormType.NORM_2),grad.norm(PETSc.NormType.NORM_2))

class Vel2TaucTikhonov(PISM.ssa.SSAFromInputFile):
  def __init__(self,input_filename,inv_data_filename,tauc_param):
    PISM.ssa.SSAFromInputFile.__init__(self,input_filename)
    self.tauc_param = tauc_param;
    self.inv_data_filename = inv_data_filename

  def _setFromOptions(self):
    PISM.ssa.SSAFromInputFile._setFromOptions(self)
    for o in PISM.OptionsGroup(PISM.Context().com,"","Vel2Tauc"):
      self.using_zeta_fixed_mask = PISM.optionsFlag("-use_zeta_fixed_mask","Keep tauc constant except where grounded ice is present",default=False)

  def _initGrid(self):
    # The implementation in PISM.ssa.SSAFromInputFile uses a non-periodic
    # grid only if the run is regional and "ssa_method=fem" in the config
    # file.  For inversions, we always use an FEM type method, so for
    # regional inversions, we always use a non-periodic grid.
    periodicity = PISM.XY_PERIODIC
    if self.is_regional:
      periodicity=PISM.NOT_PERIODIC
    PISM.util.init_grid_from_file(self.grid,self.boot_file,periodicity);

  def setup(self):

    PISM.ssa.SSAFromInputFile.setup(self)

    vecs = self.modeldata.vecs

    if vecs.has('vel_bc'):
      self.ssa.set_boundary_conditions(vecs.bc_mask,vecs.vel_bc)

    if vecs.has('zeta_fixed_mask') and self.using_zeta_fixed_mask:
      self.ssa.set_tauc_fixed_locations(vecs.zeta_fixed_mask)

    self.ssa.init(vecs.asPISMVars())

    # Cache the values of the coefficients at quadrature points once here.
    # Subsequent solves will then not need to cache these values.
    self.ssa.cacheQuadPtValues();
    self.ssa.set_functionals()


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
    gc = PISM.GeometryCalculator(sea_level, self.config)
    gc.compute(bed,thickness,mask,surface)

    if self.is_regional:
      vecs.add( PISM.util.standardNoModelMask(self.grid), 'no_model_mask' )
      vecs.no_model_mask.regrid(self.boot_file,True)
      vecs.add( vecs.surface, 'usurfstore')
      vecs.setPISMVarsName('usurfstore','usurfstore')

    if self.config.get_flag('ssa_dirichlet_bc'):
      vecs.add( PISM.util.standard2dVelocityVec( self.grid, name='_ssa_bc', desc='SSA velocity boundary condition',intent='intent' ), "vel_ssa_bc" )
      has_u_ssa_bc = PISM.util.fileHasVariable(self.boot_file,'u_ssa_bc');
      has_v_ssa_bc = PISM.util.fileHasVariable(self.boot_file,'v_ssa_bc');
      if (not has_u_ssa_bc) or (not has_v_ssa_bc):
        PISM.verbPrintf(2,grid.com, "Input file '%s' missing Dirichlet boundary data u/v_ssa_bc; using zero default instead." % self.boot_file)
        vecs.vel_ssa_bc.set(0.)
      else:
        vecs.vel_ssa_bc.regrid(self.boot_file,True)

      if self.is_regional:
        vecs.add( vecs.no_model_mask, 'bc_mask')
      else:
        vecs.add( PISM.util.standardBCMask( self.grid ), 'bc_mask' )
        bc_mask_name = vecs.bc_mask.string_attr("name")
        if PISM.util.fileHasVariable(self.boot_file,bc_mask_name):
          vecs.bc_mask.regrid(self.boot_file,True)          
        else:
          PISM.verbPrintf(2,grid.com,"Input file '%s' missing Dirichlet location mask '%s'.  Default to no Dirichlet locations." %(self.boot_file,bc_mask_name))
          vecs.bc_mask.set(0)
      # We call this variable 'bc_mask' in the python code, it is called
      # 'bcflag' when passed between pism components, and it has yet
      # another name when written out to a file.  Anyway, we flag its
      # export to PISMVars name here.
      vecs.setPISMVarsName('bc_mask','bcflag')

    vecs.add( PISM.util.standardVelocityMisfitWeight(self.grid) )
    weight = vecs.vel_misfit_weight
    weight.regrid(self.inv_data_filename,True)

    if PISM.util.fileHasVariable(self.inv_data_filename,'misfit_element_mask'):
      vecs.add( PISM.util.standardMisfitElementMask(self.grid) )
      vecs.misfit_element_mask.regrid(self.inv_data_filename,True)
    else:
      raise Exception()

    zeta_fixed_mask = PISM.IceModelVec2Int()
    zeta_fixed_mask.create(self.grid, 'zeta_fixed_mask', True, self.grid.max_stencil_width);
    zeta_fixed_mask.set_attrs("model_state", "tauc_unchanging integer mask", "", "");
    mask_values=[0,1]
    zeta_fixed_mask.set_attr("flag_values", mask_values);
    zeta_fixed_mask.set_attr("flag_meanings","tauc_changable tauc_unchangeable");
    zeta_fixed_mask.output_data_type = PISM.PISM_BYTE;
    
    zeta_fixed_mask.set(1);
    with PISM.util.Access(comm=zeta_fixed_mask,nocomm=mask):
      mq = PISM.MaskQuery(mask)
      for (i,j) in self.grid.points():
        if mq.grounded_ice(i,j):
          zeta_fixed_mask[i,j] = 0;
    vecs.add(zeta_fixed_mask)

  def _constructSSA(self):
    md = self.modeldata
    self.tauc_param = tauc_param_factory.create(self.config)
    return PISM.InvSSATikhonov(md.grid,md.basal,md.enthalpyconverter,self.tauc_param,self.config)


  def write(self,filename,append=False):
    if not append:
      PISM.ssa.SSAFromInputFile.write(self,filename)
    else:
      grid = self.grid
      vecs = self.modeldata.vecs

      pio = PISM.PIO(grid.com,grid.rank,"netcdf3")
      pio.open(filename,PISM.NC_WRITE,True) #append mode!
      
      self.modeldata.vecs.write(filename)
      pio.close()

# class ZetaSaver:
#   """Iteration listener used to save a copy of the current value
#   of zeta (i.e. parameterized tauc) at each iteration during an inversion."""
#   def __init__(self,output_filename):
#     self.output_filename = output_filename
# 
#   def __call__(self,inverse_solver,count,x,Fx,y,r):
#     zeta = x.core()
#     # The solver doesn't care what the name of zeta is, and we
#     # want it called 'zeta_inv' in the output file, so we rename it.
#     zeta.rename('zeta_inv', 'last iteration of parameterized basal yeild stress computed by inversion','')
#     zeta.write(self.output_filename)

# class Vel2TaucPlotListener(PlotListener):
#   def __init__(self,grid,Vmax,method):
#     PlotListener.__init__(self,grid)
#     self.Vmax = Vmax
#     self.rank = grid.rank
#     self.l2_weight = None
#     self.l2_weight_init = False
#     self.figure =None
#     self.method = method
# 
#   def __call__(self,inverse_solver,count,x,y,d,r,*args):
#     if self.l2_weight_init == False:
#       vecs = inverse_solver.forward_problem.ssarun.modeldata.vecs;
#       self.l2_weight=self.tz_scalar.communicate(vecs.vel_misfit_weight)
#       self.l2_weight_init = True
#     PlotListener.__call__(self,inverse_solver,count,x,y,d,r,*args)
# 
#   def iteration(self,inverse_solver,count,x,Fx,y,d,r,arg_extra,*args):      
# 
#     import matplotlib.pyplot as pp
# 
#     if self.figure is None:
#       self.figure = pp.figure()
#     else:
#       pp.figure(self.figure.number)
# 
#     l2_weight=self.l2_weight
# 
#     pp.clf()
#     
#     V = self.Vmax
# 
#     pp.subplot(2,3,1)
#     rx = l2_weight*r[0,:,:]
#     rx = np.maximum(rx,-V)
#     rx = np.minimum(rx,V)
#     pp.imshow(rx,origin='lower',interpolation='nearest')
#     pp.colorbar()
#     pp.title('r_x')
#     pp.jet()
# 
#     pp.subplot(2,3,4)
#     ry = l2_weight*r[1,:,:]
#     ry = np.maximum(ry,-V)
#     ry = np.minimum(ry,V)
#     pp.imshow(ry,origin='lower',interpolation='nearest')
#     pp.colorbar()
#     pp.title('r_y')
#     pp.jet()
#     
#     
#     if self.method == 'ign':
#       Td = arg_extra
#       pp.subplot(2,3,2)
#       Tdx = Td[0,:,:]
#       pp.imshow(Tdx,origin='lower',interpolation='nearest')
#       pp.colorbar()
#       pp.title('Td_x')
#       pp.jet()
# 
#       pp.subplot(2,3,5)
#       Tdy = Td[1,:,:]
#       pp.imshow(Tdy,origin='lower',interpolation='nearest')
#       pp.colorbar()
#       pp.title('Td_y')
#       pp.jet()
#     else:
#       TStarR = arg_extra
#       pp.subplot(2,3,2)
#       pp.imshow(TStarR,origin='lower',interpolation='nearest')
#       pp.colorbar()
#       pp.title('TStarR')
#       pp.jet()
# 
#     d *= -1
#     pp.subplot(2,3,3)      
#     pp.imshow(d,origin='lower',interpolation='nearest')
#     pp.colorbar()
#     pp.jet()
#     pp.title('-d')
# 
#     pp.subplot(2,3,6)      
#     pp.imshow(x,origin='lower',interpolation='nearest')
#     pp.colorbar()
#     pp.jet()
#     pp.title('zeta')
#     
#     pp.ion()
#     pp.show()

# class MonitorAdjoint:
#   def __init__(self):
#     self.Td = None
#     self.TStarR = None
# 
#   def __call__(self,inverse_solver,count,x,Fx,y,d,r,*args):
#     fp = inverse_solver.forward_problem
#     self.Td=fp.T(d,self.Td)
#     self.TStarR=forward_problem.TStar(r,out=self.TStarR)
#     ip1 = fp.domainIP(d,self.TStarR)
#     ip2 = fp.rangeIP(self.Td,r)
#     siple.reporting.msg("adjoint test: <Td,r>=%g <d,T^*r>=%g (percent error %g)",ip1,ip2,(abs(ip1-ip2))/max(abs(ip1),abs(ip2)))
# 
# class MonitorAdjointLin:
#   def __init__(self):
#     self.Td = None
#     self.TStarR = None
# 
#   def __call__(self,inverse_solver,count,x,y,d,r,*args):
#     fp = inverse_solver.forward_problem
#     self.Td=fp.T(d,self.Td)
#     self.TStarR=forward_problem.TStar(r,out=self.TStarR)
#     ip1 = fp.domainIP(d,self.TStarR)
#     ip2 = fp.rangeIP(self.Td,r)
#     siple.reporting.msg("adjoint test: <Td,r>=%g <d,T^*r>=%g (percent error %g)",ip1,ip2,(abs(ip1-ip2))/max(abs(ip1),abs(ip2)))


# class Vel2TaucLinPlotListener(LinearPlotListener):
#   def __init__(self,grid,Vmax):
#     LinearPlotListener.__init__(self,grid)
#     self.Vmax = Vmax
#     self.l2_weight = None
#     self.l2_weight_init = False
#     self.figure =None
# 
#   def __call__(self,inverse_solver,count,x,y,d,r,*args):
#     # On the first go-around, extract the l2_weight vector onto 
#     # processor zero.
#     if self.l2_weight_init == False:
#       vecs = inverse_solver.forward_problem.ssarun.modeldata.vecs;
#       self.l2_weight=self.tz_scalar.communicate(vecs.vel_misfit_weight)
#       self.l2_init = True
#     LinearPlotListener.__call__(self,solver,count,x,y,d,r,*args)
# 
#   def iteration(self,inverse_solver,count,x,y,d,r,*args):      
#     import matplotlib.pyplot as pp
# 
#     if self.figure is None:
#       self.figure = pp.figure()
#     else:
#       pp.figure(self.figure.number)
# 
#     l2_weight=self.l2_weight
# 
#     pp.clf()
#     
#     V = self.Vmax
#     pp.subplot(1,3,1)
#     rx = l2_weight*r[0,:,:]
#     rx = np.maximum(rx,-V)
#     rx = np.minimum(rx,V)
#     pp.imshow(rx,origin='lower',interpolation='nearest')
#     pp.colorbar()
#     pp.title('ru')
#     pp.jet()
#     
#     pp.subplot(1,3,2)
#     ry = l2_weight*r[1,:,:]
#     ry = np.maximum(ry,-V)
#     ry = np.minimum(ry,V)
#     pp.imshow(ry,origin='lower',interpolation='nearest')
#     pp.colorbar()
#     pp.title('rv')
#     pp.jet()
#     
#     d *= -1
#     pp.subplot(1,3,3)      
#     pp.imshow(d,origin='lower',interpolation='nearest')
#     pp.colorbar()
#     pp.jet()
#     pp.title('-d')
#     
#     pp.ion()
#     pp.show()


## Main code starts here
if __name__ == "__main__":
  context = PISM.Context()
  config = context.config
  com = context.com
  PISM.set_abort_on_sigint(True)

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
    method = PISM.optionsList(context.com,"-inv_method","Inversion algorithm",["nlcg","ign","sd"],"ign")
    # rms_error = PISM.optionsReal("-rms_error","RMS velocity error",default=default_rms_error)
    eta = PISM.optionsReal("-eta","tikhonov penalty parameter",default=default_eta)
    inv_ssa_cL2 = PISM.optionsReal("-inv_ssa_cL2","L2 coefficient for domain inner product",default=default_ssa_l2_coeff)
    inv_ssa_cH1 = PISM.optionsReal("-inv_ssa_cH1","H1 coefficient for domain inner product",default=default_ssa_h1_coeff)
    do_plotting = PISM.optionsFlag("-inv_plot","perform visualization during the computation",default=False)
    do_final_plot = PISM.optionsFlag("-inv_final_plot","perform visualization at the end of the computation",default=False)
    do_pause = PISM.optionsFlag("-inv_pause","pause each iteration",default=False)
    test_adjoint = PISM.optionsFlag("-inv_test_adjoint","Test that the adjoint is working",default=False)
    ls_verbose = PISM.optionsFlag("-inv_ls_verbose","Turn on a verbose linesearch.",default=False)
    do_restart = PISM.optionsFlag("-inv_restart","Restart a stopped computation.",default=False)
    use_tauc_prior = PISM.optionsFlag("-use_tauc_prior","Use tauc_prior from inverse data file as initial guess.",default=False)
    ign_theta  = PISM.optionsReal("-ign_theta","theta parameter for IGN algorithm",default=0.5)
    Vmax = PISM.optionsReal("-inv_plot_vmax","maximum velocity for plotting residuals",default=30)
    monitor_adjoint = PISM.optionsFlag("-inv_monitor_adjoint","Track accuracy of the adjoint during computation",default=False)
    is_regional = PISM.optionsFlag("-regional","Compute SIA/SSA using regional model semantics",default=False)
    prep_module = PISM.optionsString("-inv_prep_module","Python module used to do final setup of inverse solver",default=None)
    
  if output_filename is None:
    output_filename = "vel2tauc_"+os.path.basename(input_filename)    

  saving_inv_data = (inv_data_filename != output_filename)

  length_scale   = 100e3 # 100 km
  area_scale     = length_scale*length_scale
  velocity_scale = 100 / PISM.secpera # 100m/a
  inv_ssa_cL2 /= area_scale 
  #eta /= velocity_scale**2

  config.set_string("inv_ssa_tauc_param","ident")
  config.set("inv_ssa_cL2",inv_ssa_cL2)
  config.set("inv_ssa_cH1",inv_ssa_cH1)

  stress_scale = 50000 # Pa
  config.set("tauc_param_trunc_tauc0",.1*stress_scale)
  config.set("tauc_param_tauc_eps",.001*stress_scale)
  config.set("tauc_param_tauc_scale",stress_scale)
  config.set("inv_ssa_velocity_scale",1) # m/a

  tauc_param = tauc_param_factory.create(config)

  PISM.setVerbosityLevel(verbosity)
  vel2tauc = Vel2TaucTikhonov(input_filename,inv_data_filename,tauc_param)
  vel2tauc.setup()

  modeldata = vel2tauc.modeldata
  vecs = modeldata.vecs
  grid = modeldata.grid

  # Determine the prior guess for tauc. This can be one of 
  # a) tauc from the input file (default)
  # b) tauc_prior from the inv_datafile if -use_tauc_prior is set
  tauc_prior = PISM.util.standardYieldStressVec(grid,'tauc_prior')
  tauc_prior.set_attrs("diagnostic", "initial guess for (pseudo-plastic) basal yield stress in an inversion", "Pa", "");
  tauc = PISM.util.standardYieldStressVec(grid)
  if use_tauc_prior:
    tauc_prior.regrid(inv_data_filename,True)
  else:
    if not PISM.util.fileHasVariable(input_filename,"tauc"):
      PISM.verbPrintf(1,com,"Initial guess for tauc is not available as 'tauc' in %s.\nYou can provide an initial guess as 'tauc_prior' using the command line option -use_tauc_prior." % input_filename)
      exit(1)
    tauc.regrid(input_filename,True)
    tauc_prior.copy_from(tauc)
  vecs.add(tauc_prior,writing=saving_inv_data)

  # If the inverse data file has a variable tauc_true, this is probably
  # a synthetic inversion.  We'll load it now so that it will get written
  # out, if needed, at the end of the computation in the output file.
  if PISM.util.fileHasVariable(inv_data_filename,"tauc_true"):
    tauc_true = PISM.util.standardYieldStressVec(grid,'tauc_true')
    tauc_true.regrid(inv_data_filename,True)
    tauc_true.read_attributes(inv_data_filename)
    vecs.add(tauc_true,writing=saving_inv_data)


  # Determine the initial guess for zeta.  If we are not
  # restarting, we convert tauc_prior to zeta.  If we are restarting,
  # we load in zeta from the output file.
  zeta = PISM.IceModelVec2S();
  zeta.create(grid, "zeta_inv", True, PISM.util.WIDE_STENCIL)
  if do_restart:
    # Just to be sure, verify that we have a 'zeta_inv' in the output file.
    if not PISM.util.fileHasVariable(output_filename,'zeta_inv'):
      PISM.verbPrintf(1,com,"Unable to restart computation: file %s is missing variable 'zeta_inv'", output_filename)
      exit(1)
    zeta.regrid(output_filename,True)
  else:
    tauc_param.convertFromTauc(tauc_prior,zeta)
  vecs.add(zeta,writing=True) # Ensure that the last value of zeta will
                              # be written out

  vel_ssa_observed = None
  vel_ssa_observed = PISM.util.standard2dVelocityVec(grid,'_ssa_observed',stencil_width=2)
  if PISM.util.fileHasVariable(inv_data_filename,"u_ssa_observed"):
    vel_ssa_observed.regrid(inv_data_filename,True)
    vecs.add(vel_ssa_observed,writing=saving_inv_data)
  else:
    if not PISM.util.fileHasVariable(inv_data_filename,"u_surface_observed"):
      PISM.verbPrintf(1,context.com,"Neither u/v_ssa_observed nor u/v_surface_observed is available in %s.\nAt least one must be specified.\n" % inv_data_filename)
      exit(1)
    vel_surface_observed = PISM.util.standard2dVelocityVec(grid,'_surface_observed',stencil_width=2)
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

  # # Prep the output file from the grid so that we can save zeta to it during the runs.
  # if not append_mode:
  #   if do_restart:
  #     logger = CarefulCaptureLogger(output_filename);
  #     
  #   pio = PISM.PIO(grid.com,grid.rank,"netcdf3")
  #   pio.open(output_filename,PISM.NC_WRITE,False)
  #   pio.def_time(grid.config.get_string("time_dimension_name"),
  #                grid.config.get_string("calendar"), grid.time.units())
  #   pio.append_time(grid.config.get_string("time_dimension_name"),grid.time.current())
  #   pio.close()
  zeta.write(output_filename)

  # Log the command line to the output file now so that we have a record of
  # what was attempted
  PISM.util.writeProvenance(output_filename)    

  # # If we haven't set up the aforementioned logger yet, it's safe to do so now.
  # if logger is None: logger = CarefulCaptureLogger(output_filename);

  ip     = PISM.InvSSATikhonovProblem(vel2tauc.ssa,zeta,vel_ssa_observed,eta)
  l = Listener()
  ip.addListener(l)
  solver = PISM.InvSSATikhonovSolver(grid.com,"tao_lmvm",ip)
  # # Attach various iteration listeners to the solver as needed for:
  # # Plotting
  # if do_plotting:
  #   solver.addIterationListener(Vel2TaucPlotListener(grid,Vmax,method))
  #   if method=='ign':
  #     solver.addLinearIterationListener(Vel2TaucLinPlotListener(grid,Vmax))
  # if monitor_adjoint:
  #   solver.addIterationListener(MonitorAdjoint())
  #   if method=='ign':
  #     solver.addLinearIterationListener(MonitorAdjointLin())
  # # Pausing
  # if do_pause:
  #   solver.addIterationListener(pauseListener)            
  # # Saving the current iteration
  # solver.addXUpdateListener(ZetaSaver(output_filename)) 

  # Solver is set up.  Give the user's prep module a chance to do any final
  # setup.
  
  # if prep_module is not None:
  #   exec "import %s as user_prep_module" % prep_module
  #   user_prep_module.prep_solver(solver,method)

  # Run the inverse solver!
  if do_restart:
    siple.reporting.msg('************** Restarting inversion. ****************')
  else:
    siple.reporting.msg('============== Starting inversion. ==================')  

  # rms_error /= PISM.secpera # m/s
  # (zeta,u) = solver.solve(zeta,vel_ssa_observed,rms_error)

  # Try solving
  if not solver.solve():
    PISM.verbPrintf(1,grid.com,"Inverse solve FAILURE (%s)!\n" % solver.reasonDescription());
    PISM.PISMEnd()
  PISM.verbPrintf(1,grid.com,"Inverse solve converged (%s)!\n" % solver.reasonDescription());

  vel_ssa_inv      = ip.stateSolution();
  zeta_inv         = ip.designSolution();

  vel_ssa_misfit = PISM.util.standard2dVelocityVec(grid,'_ssa_misfit',stencil_width=2)
  vel_ssa_misfit.copy_from(vel_ssa_inv)
  vel_ssa_misfit.add(-1,vel_ssa_observed)
  vel_ssa_misfit.beginGhostComm();  vel_ssa_misfit.endGhostComm()
  misfit = math.sqrt(vel2tauc.ssa.evalPenalty(vel_ssa_misfit))
  print "Velocity misfit %g m/a" % (misfit * PISM.secpera)

  # Convert back from zeta to tauc
  tauc_param.convertToTauc(zeta,tauc)

  # It may be that a 'tauc' was read in earlier.  We replace it with
  # our newly generated one.
  if vecs.has('tauc'): vecs.remove('tauc')
  vecs.add(tauc,writing=True)

  vel_ssa_inv.rename("_ssa_inv","SSA velocity computed by inversion","")
  vecs.add(vel_ssa_inv,writing=True)

  # Write solution out to netcdf file
  vel2tauc.write(output_filename,append=append_mode)

  # # If we're not in append mode, the previous command just nuked
  # # the output file.  So we rewrite the siple log.
  # if not append_mode:
  #   logger.write(output_filename)
