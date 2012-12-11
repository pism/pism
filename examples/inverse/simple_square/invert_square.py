import sys, petsc4py
petsc4py.init(sys.argv)
from petsc4py import PETSc
import numpy as np
import siple

import PISM
from PISM import util

from PISM.invert_ssa import SSAForwardProblem, InvSSARun, tauc_params, InvertSSANLCG, InvertSSAIGN
from PISM.sipletools import PISMLocalVector as PLV

# siple.reporting.clear_loggers()
# siple.reporting.add_logger(pism_print_logger)
# siple.reporting.set_pause_callback(pism_pause)

p_schoof = 4.0/3.0   # = 1 + 1/n

rms_error = 2.754712829589844 #m/a # invParams.pointwise_error_size from conv_rate.py run when relative_noise_level is set to 1%
#rms_error /= PISM.secpera # m/s

ssa_l2_coeff = 1.
ssa_h1_coeff = 0.

tauc_guess_scale = 1.
tauc_guess_const = None # will be changed to avg value of tauc_true below

# PISM grids go from -Lx to +Lx and -Ly to +Ly
Lx = 40.e3 # meters
Ly = 80.e3 # meters
Mx = 61
My = 121  #int(Mx/Lx * Ly)

B = 7. * (PISM.secpera)**(1./3) * 1.e5 # units have to be kg, m, sec in PISM

def testi_tauc(grid,ice,tauc):
  xzero = -20.e3
  yzero = 30.e3
  xsig = 1.e8 # = 2*\sigma^2 where \sigma is the std
  ysig = xsig
  with PISM.util.Access(comm=tauc):
    for (i,j) in grid.points():
      x=grid.x[i]
      y=grid.y[j]
      tauc[i,j] = 5.e-4 + 5.e-3*np.exp(-(x-xzero)**2/xsig - (y-yzero)**2/ysig)
      tauc[i,j] = tauc[i,j] * PISM.secpera * 1.e5

class testi_run(InvSSARun):
  def __init__(self,Mx,My):
    self.grid = PISM.Context().newgrid()
    self.Mx = Mx
    self.My = My

  def _initGrid(self):
    Mx=self.Mx; My=self.My
    PISM.util.init_shallow_grid(self.grid,Lx,Ly,Mx,My,PISM.NONE); # NONE makes a non-periodic grid

  def _initPhysics(self):
    config = self.config
    config.set_flag("do_pseudo_plastic_till", True)
    config.set("pseudo_plastic_q", 1.0);
    config.set("pseudo_plastic_uthreshold", 1.0/PISM.secpera) # so that tau_b = tauc * u

    basal = PISM.IceBasalResistancePseudoPlasticLaw(config)

    # irrelevant
    enthalpyconverter = PISM.EnthalpyConverter(config);

    ice = PISM.CustomGlenIce(self.grid.com, "", config, enthalpyconverter);
    ice.setHardness(B)

    self.modeldata.setPhysics(ice,basal,enthalpyconverter)

  def _initSSACoefficients(self):
    vecs = self.modeldata.vecs; grid = self.grid
    vecs.add( util.standardIceThicknessVec( grid ), 'thickness')
    vecs.add( util.standardBedrockElevationVec(grid), 'bed')
    vecs.add( util.standardYieldStressVec( grid ), 'tauc')
    vecs.add( util.standardEnthalpyVec( grid ), 'enthalpy' )
    vecs.add( util.standardIceMask( grid ), 'ice_mask' )
    # vecs.add( util.standardDrivingStressX( grid ), 'driving_x' )
    # vecs.add( util.standardDrivingStressY( grid ), 'driving_y' )
    vecs.add( util.standardIceSurfaceVec( grid ), 'surface_altitude' )
    vecs.add( util.standardVelocityMisfitWeight(grid) )

    self._allocateBCs()

    vecs.ice_mask.set(PISM.MASK_GROUNDED)

    # why is tauc_true in here???
    testi_tauc(self.modeldata.grid,self.modeldata.ice,vecs.tauc)

    grid = self.grid

    vars = [vecs.bed,vecs.thickness,vecs.surface_altitude]
    with PISM.util.Access(comm=vars):
      for (i,j) in grid.points():
        x=grid.x[i]; y=grid.y[j]
        vecs.thickness[i,j] = 1060. + 2.e-3*y # 900. + 2.e-3*y: different zero crossing because of different origin
        vecs.bed[i,j] = 0.
        vecs.surface_altitude[i,j] = vecs.bed[i,j] + vecs.thickness[i,j]


    # standard_gravity = grid.config.get("standard_gravity");
    # with PISM.util.Access(comm=[vecs.driving_y,vecs.driving_x]):
    #   for (i,j) in grid.points():
    #     vecs.driving_y[i,j] = 0.
    #     vecs.driving_x[i,j] = self.modeldata.ice.rho*standard_gravity*vecs.thickness[i,j]*slope

    # Dirichlet bc of 0 everywhere except for the lower boundary
    with PISM.util.Access(comm=[vecs.bc_mask,vecs.vel_bc]):
      for (i,j) in grid.points():
        edge = ( (j == grid.My-1) ) or ( (i==0) or (i==grid.Mx-1) ); # (j == 0) or
        if (edge):
          vecs.bc_mask[i,j] = 1
          vecs.vel_bc[i,j].u = 0.
          vecs.vel_bc[i,j].v = 0.

    # Set misfit_weight to one everywhere
    misfit_weight=vecs.vel_misfit_weight
    with PISM.util.Access(comm=misfit_weight):
      for (i,j) in grid.points():
        misfit_weight[i,j] = 1.


## Main code starts here
context = PISM.Context()
config = context.config
PISM.set_abort_on_sigint(True)

for o in PISM.OptionsGroup(context.com,"","Invert test of simple square"):
  Mx = PISM.optionsInt("-Mx","Number of grid points in x-direction",default=Mx)
  My = PISM.optionsInt("-My","Number of grid points in y-direction",default=My)
  output_file = PISM.optionsString("-o","output file",default="simple_square.nc")
  verbosity = PISM.optionsInt("-verbose","verbosity level",default=2)
  method = PISM.optionsList(context.com,"-inv_method","Inversion algorithm",["nlcg","ign","sd"],"ign")
  rms_error = PISM.optionsReal("-rms_error","RMS velocity error",default=rms_error)
  tauc_param_type = PISM.optionsList(context.com,"-tauc_param","zeta->tauc parameterization",["ident","square","exp"],"ident")
  ssa_l2_coeff = PISM.optionsReal("-ssa_l2_coeff","L2 coefficient for domain inner product",default=ssa_l2_coeff)
  ssa_h1_coeff = PISM.optionsReal("-ssa_h1_coeff","H1 coefficient for domain inner product",default=ssa_h1_coeff)
  tauc_guess_scale = PISM.optionsReal("-tauc_guess_scale","initial guess for tauc to be this factor of the true value",default=tauc_guess_scale)
  tauc_guess_const = PISM.optionsReal("-tauc_guess_const","initial guess for tauc to be this constant",default=tauc_guess_const)
  do_plotting = PISM.optionsFlag("-inv_plot","perform visualization during the computation",default=False)
  do_final_plot = PISM.optionsFlag("-inv_final_plot","perform visualization at the end of the computation",default=True)
  do_pause = PISM.optionsFlag("-inv_pause","pause each iteration",default=False)
  test_adjoint = PISM.optionsFlag("-inv_test_adjoint","Test that the adjoint is working",default=False)
  ls_verbose = PISM.optionsFlag("-inv_ls_verbose","Turn on a verbose linesearch.",default=False)

config.set_string("inv_ssa_tauc_param",tauc_param_type)
config.set("inv_ssa_domain_l2_coeff",ssa_l2_coeff)
config.set("inv_ssa_domain_h1_coeff",ssa_h1_coeff)

tauc_param = tauc_params[tauc_param_type]

PISM.setVerbosityLevel(verbosity)
testi = testi_run(Mx,My)
testi.setup()

forward_problem = SSAForwardProblem(testi)

grid = testi.grid

# Build the true yield stress for test I
tauc_true = PISM.util.standardYieldStressVec(grid,name="tauc_true")
testi_tauc(grid,testi.modeldata.ice,tauc_true)

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

# Build the initial guess for tauc for the inversion.
tauc = PISM.util.standardYieldStressVec(grid)
testi_tauc(grid,testi.modeldata.ice,tauc)
tauc_guess_const = tauc.sum() / (grid.xm * grid.ym) # mean of tauc_true
if not tauc_guess_const is None:
  tauc.set(tauc_guess_const)
else:
  testi_tauc(grid,testi.modeldata.ice,tauc)
  tauc.scale(tauc_guess_scale)

# Convert tauc guess to zeta guess
if config.get_string('inv_ssa_tauc_param')=='ident':
  zeta = tauc
else:
  zeta = PISM.IceModelVec2S();
  zeta.create(grid, "zeta", True, PISM.util.WIDE_STENCIL)
  with PISM.util.Access(nocomm=[tauc],comm=[zeta]):
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

params=Solver.defaultParameters()
if method == "sd":
  params.steepest_descent = True
  params.ITER_MAX=10000
elif method =="ign":
  params.linearsolver.ITER_MAX=10000
  params.linearsolver.verbose = True
if ls_verbose:
  params.linesearch.verbose = True
params.verbose   = True
params.deriv_eps = 0.


# Run the inversion
solver=Solver(forward_problem,params=params)
if do_plotting:
  solver.addIterationListener(TestIPlotListener(grid))
if do_pause:
  solver.addIterationListener(pauseListener)

rms_error /= PISM.secpera # m/s
(zeta,u) = solver.solve(zeta,u_true,rms_error)

# Convert back from zeta to tauc
if config.get_string('inv_ssa_tauc_param')=='ident':
  tauc = zeta
else:
  with PISM.util.Access(nocomm=[zeta],comm=[tauc]):
    for (i,j) in grid.points():
      (tauc[i,j],dummy) = tauc_param.toTauc(zeta[i,j])

# Write solution out to netcdf file
testi.write(output_file)
tauc.write(output_file)
tauc_true.write(output_file)
u.set_name("_computed",0)
u.write_in_glaciological_units=True
u.write(output_file)
u_true.write_in_glaciological_units=True
u_true.write(output_file)

# forward_problem.solver.vel_bc.write(output_file)
