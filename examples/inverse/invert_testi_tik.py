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
import math

import PISM
from PISM import util

from PISM.invert_ssa import InvSSARun, SSAForwardProblem, InvertSSANLCG, InvertSSAIGN, \
tauc_param_factory, PlotListener 

import matplotlib.pyplot as pp

class InvSSATikRun(PISM.ssa.SSARun):

  def setup(self):

    PISM.ssa.SSARun.setup(self)

    vecs = self.modeldata.vecs

    # The SSA instance will not keep a reference to pismVars; it only uses it to extract
    # its desired variables.  So it is safe to pass it pismVars and then let pismVars
    # go out of scope at the end of this method.

    self.ssa.init(vecs.asPISMVars())

    if vecs.has('vel_bc'):
      self.ssa.set_boundary_conditions(vecs.bc_mask,vecs.vel_bc)

    # Cache the values of the coefficients at quadrature points once here.
    # Subsequent solves will then not need to cache these values.
    self.ssa.cacheQuadPtValues();
    self.ssa.set_functionals()

  def _constructSSA(self):
    md = self.modeldata
    vecs  = self.modeldata.vecs
    self.tauc_param = tauc_param_factory.create(self.config)
    return PISM.InvSSATikhonov(md.grid,md.basal,md.enthalpyconverter,self.tauc_param,self.config)


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

ssa_l2_coeff = 1.
ssa_h1_coeff = 0.

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

class testi_run(InvSSATikRun):
  def __init__(self,Mx,My):
    self.grid = PISM.Context().newgrid()
    self.Mx = Mx
    self.My = My

  def _initGrid(self):
    Mx=self.Mx; My=self.My
    Ly = 3*L_schoof   # 300.0 km half-width (L=40.0km in Schoof's choice of variables)
    Lx = max(60.0e3, ((Mx - 1) / 2.) * (2.0 * Ly / (My - 1)) )
    PISM.util.init_shallow_grid(self.grid,Lx,Ly,Mx,My,PISM.NOT_PERIODIC);

  def _initPhysics(self):
    config = self.config
    basal = PISM.IceBasalResistancePlasticLaw(
         config.get("plastic_regularization") / PISM.secpera,
         config.get_flag("do_pseudo_plastic_till"),
         config.get("pseudo_plastic_q"),
         config.get("pseudo_plastic_uthreshold") / PISM.secpera);

    basal.printInfo(1,PETSc.COMM_WORLD)
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
    method = PISM.optionsList(context.com,"-inv_method","Inversion algorithm",["nlcg","ign","sd"],"ign")
    eta = PISM.optionsReal("-eta","penalty weight",default=1)
    right_side_weight = PISM.optionsReal("-right_side_weight","L2 weight for y>0",default=right_side_weight)
    inv_ssa_cL2 = PISM.optionsReal("-inv_ssa_cL2","L2 coefficient for domain inner product",default=1)
    inv_ssa_cH1 = PISM.optionsReal("-inv_ssa_cH1","H1 coefficient for domain inner product",default=1)
    tauc_guess_scale = PISM.optionsReal("-tauc_guess_scale","initial guess for tauc to be this factor of the true value",default=tauc_guess_scale)
    tauc_guess_const = PISM.optionsReal("-tauc_guess_const","initial guess for tauc to be this constant",default=tauc_guess_const)
    do_plotting = PISM.optionsFlag("-inv_plot","perform visualization during the computation",default=False)
    do_final_plot = PISM.optionsFlag("-inv_final_plot","perform visualization at the end of the computation",default=True)
    do_pause = PISM.optionsFlag("-inv_pause","pause each iteration",default=False)
    test_adjoint = PISM.optionsFlag("-inv_test_adjoint","Test that the adjoint is working",default=False)
    ls_verbose = PISM.optionsFlag("-inv_ls_verbose","Turn on a verbose linesearch.",default=False)

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


  print stress_scale, velocity_scale, time_scale, strainrate_scale, viscosity_scale, nuH_scale

  # area_scale = 1000e2*1000e2  # 100 km^2
  # velocity_scale = 100 / PISM.secpera  # 100 m/a
  # stress_scale = 500000 # Pa

  inv_ssa_cL2 /= area_scale 
  eta /= velocity_scale**2

  config.set("inv_ssa_cL2",inv_ssa_cL2)
  config.set("inv_ssa_cH1",inv_ssa_cH1)

  config.set("tauc_param_trunc_tauc0",.1*stress_scale)
  config.set("tauc_param_tauc_eps",.001*stress_scale)
  config.set("tauc_param_tauc_scale",stress_scale)

  # Default to be overridden by command line perhaps
  config.set_string("inv_ssa_tauc_param","ident")

  PISM.setVerbosityLevel(verbosity)
  testi = testi_run(Mx,My)
  testi.setup()

  forward_problem = SSAForwardProblem(testi)

  grid = testi.grid

  # Build the true yeild stress for test I
  tauc_true = PISM.util.standardYieldStressVec(grid,name="tauc_true")
  testi_tauc(grid, tauc_true)
  tauc_true.scale(1.)
  # tauc_true.set(stress_scale)
  # Convert tauc_true to zeta_true
  zeta_true = PISM.IceModelVec2S();
  zeta_true.create(grid,"zeta_true",PISM.kHasGhosts,kFEMStencilWidth)
  tauc_param = PISM.invert_ssa.tauc_param_factory.create(config)
  tauc_param.convertFromTauc(tauc_true,zeta_true)

  # Send the true yeild stress through the forward problem to 
  # get at true velocity field.
  if not testi.ssa.linearizeAt(zeta_true):
    PISM.verbPrintf(1,grid.com,"Forward solve failed (%s)!\n" % testi.ssa.reasonDescription());
    exit(1)
  u_obs = PISM.util.standard2dVelocityVec( grid, name='_ssa_true', desc='SSA velocity boundary condition',intent='intent' )
  u_obs.copy_from(testi.ssa.solution())


  tauc_guess_const = stress_scale
  # Build the initial guess for tauc for the inversion.
  tauc = PISM.util.standardYieldStressVec(grid)
  if not tauc_guess_const is None:
    tauc.set(tauc_guess_const)
  else:
    testi_tauc(grid, tauc)
    tauc.scale(tauc_guess_scale)

  designFunc =PISM.H1NormFunctional2S(grid,inv_ssa_cL2,inv_ssa_cH1)
  stateFunc = PISM.MeanSquareObservationFunctional2V(grid)
  print "cL2, cH1", inv_ssa_cL2,inv_ssa_cH1
  grid.printInfo(1)
  print "design functional value is: %g" % designFunc.valueAt(zeta_true)
  print "tauc design functional value is: %g" % designFunc.valueAt(tauc_true)
  print "state functional value is: %g; eta=%g" % (eta*stateFunc.valueAt(u_obs),eta)
  
  maxU = 0 
  maxZeta = 0
  maxTauc = 0
  with PISM.util.Access(nocomm=[u_obs,zeta_true,tauc_true]):
    for (i,j) in grid.points():
      U=u_obs[i,j].magnitude()
      if U>maxU:
        maxU = U
      Zeta=zeta_true[i,j]
      if Zeta>maxZeta:
        maxZeta = Zeta
      TAUC=tauc_true[i,j]
      if TAUC>maxTauc:
        maxTauc = TAUC
  print "Maximum velocity is : %g m/s = %g m/a.  Relative: %g" % (maxU,maxU*PISM.secpera,maxU/velocity_scale)
  print "Maximum zeta is : %g Relative: %g" % (maxZeta,maxZeta)
  print "Maximum tauc is : %g Relative: %g" % (maxTauc,maxTauc/stress_scale)

  # tozero = PISM.toproczero.ToProcZero(grid,dof=2,dim=2)
  # u_obs_a = tozero.communicate(u_obs)
  # pp.plot(grid.y,u_obs_a[0,:,Mx/2])
  # pp.show()
  # quit()


  # Convert tauc guess to zeta guess
  zeta = PISM.IceModelVec2S();
  zeta.create(grid, "zeta", PISM.kHasGhosts, kFEMStencilWidth)
  tauc_param.convertFromTauc(tauc,zeta)


  # Build the inverse solver.  This first step feels redundant.
  ip = PISM.InvSSATikhonovProblem(testi.ssa,zeta,u_obs,eta)
  solver = PISM.InvSSATikhonovSolver(grid.com,"tao_lmvm",ip)

  # Try solving
  if not solver.solve():
    PISM.verbPrintf(1,grid.com,"Inverse solve FAILURE (%s)!\n" % solver.reasonDescription());
    
  PISM.verbPrintf(1,grid.com,"Inverse solve success (%s)!\n" % solver.reasonDescription());

  u_i    = ip.stateSolution();
  zeta_i = ip.designSolution();

  du = PISM.IceModelVec2V()
  du.create(grid, "du", PISM.kHasGhosts, kFEMStencilWidth );
  du.copy_from(u_i)
  du.add(-1,u_obs)
  du.beginGhostComm(); du.endGhostComm()
  misfit_functional = PISM.MeanSquareObservationFunctional2V(grid,testi.modeldata.vecs.vel_misfit_weight);
  misfit_functional.normalize()
  misfit = math.sqrt(misfit_functional.valueAt(du))
  misfit_norm = math.sqrt(misfit_functional.valueAt(u_i))  
  PISM.verbPrintf(1,grid.com,"RMS Misfit: %g (relative %g)\n",misfit,misfit/misfit_norm)


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
  if do_final_plot and (not tauc_a is None):
    from matplotlib import pyplot
    pyplot.clf()
    pyplot.plot(grid.y,tauc_a[:,Mx/2])
    pyplot.plot(grid.y,tauc_true[:,Mx/2])
    pyplot.show()
