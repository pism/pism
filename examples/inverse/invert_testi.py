import sys, petsc4py
petsc4py.init(sys.argv)
from petsc4py import PETSc
import numpy as np
import tozero
import siple

import PISM

from pismssaforward import PISMSSAForwardProblem, InvertSSANLCG, InvertSSAIGN, tauc_params
from linalg_pism import PISMLocalVector

def pism_print_logger(message,severity):
  verb = severity
  com = PISM.Context().com
  PISM.verbPrintf(verb,com, "%s\n" % message)
siple.reporting.clear_loggers()
siple.reporting.add_logger(pism_print_logger)

Mx = 11 
My = 61

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

tauc_guess_scale = 0.2
tauc_guess_const = None

def testi_tauc(grid,ice,tauc):
  standard_gravity = grid.config.get("standard_gravity");
  f = ice.rho*standard_gravity*H0_schoof*slope
  with PISM.util.Access(comm=tauc):
    for (i,j) in grid.points():
      y=grid.y[j]
      tauc[i,j] = f* (abs(y/L_schoof)**m_schoof)

class testi(PISMSSAForwardProblem):
  def __init__(self,Mx,My):
    PISMSSAForwardProblem.__init__(self)
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
    basal = PISM.IceBasalResistancePlasticLaw(
         config.get("plastic_regularization") / PISM.secpera,
         config.get_flag("do_pseudo_plastic_till"),
         config.get("pseudo_plastic_q"),
         config.get("pseudo_plastic_uthreshold") / PISM.secpera);
 
    # irrelevant
    enthalpyconverter = PISM.EnthalpyConverter(config);

    ice = PISM.CustomGlenIce(self.grid.com, "", config, enthalpyconverter);
    ice.setHardness(B_schoof)
    
    self.solver.setPhysics(ice,basal,enthalpyconverter)


  def _initSSACoefficients(self):
    solver = self.solver
    solver.allocateCoeffs(using_l2_weight=True,using_explicit_driving_stress=True)
    solver.allocateBCs()
    solver.thickness.set(H0_schoof)
    solver.ice_mask.set(PISM.MASK_GROUNDED)
    solver.bed.set(0.)

    testi_tauc(solver.grid,solver.ice,solver.tauc)

    grid = self.grid
    standard_gravity = grid.config.get("standard_gravity");
    f = self.solver.ice.rho*standard_gravity*H0_schoof*slope
    driving_stress = solver.drivingstress
    with PISM.util.Access(comm=[driving_stress,solver.bc_mask,solver.vel_bc]):
      for (i,j) in grid.points():
        driving_stress[i,j].u = f
        driving_stress[i,j].v = 0

        if (j == 0) or (j==grid.My-1):
          solver.bc_mask[i,j]=1
          solver.vel_bc[i,j].u=0
          solver.vel_bc[i,j].v=0

    l2_weight=solver.range_l2_weight
    with PISM.util.Access(comm=l2_weight):
      for (i,j) in grid.points():
        if grid.y[j] < 0:
          l2_weight[i,j] = 1.;
        else:
          l2_weight[i,j] = right_side_weight;

## Main code starts here

context = PISM.Context()
config = context.config()
PISM.set_abort_on_sigint(True)

for o in PISM.OptionsGroup(context.com,"","Invert test I"):
  Mx = PISM.optionsInt("-Mx","Number of grid points in x-direction",default=Mx)
  My = PISM.optionsInt("-My","Number of grid points in y-direction",default=My)
  output_file = PISM.optionsString("-o","output file",default="invert_testi.nc")
  verbosity = PISM.optionsInt("-verbose","verbosity level",default=2)
  method = PISM.optionsList(context.com,"-inv_method","Inversion algorithm",["nlcg","ign","sd"],"ign")
  rms_error = PISM.optionsReal("-rms_error","RMS velocity error",default=rms_error)
  right_side_weight = PISM.optionsReal("-right_side_weight","L2 weight for y>0",default=right_side_weight)
  tauc_param_type = PISM.optionsList(context.com,"-tauc_param","zeta->tauc parameterization",["ident","square","exp"],"ident")
  ssa_l2_coeff = PISM.optionsReal("-ssa_l2_coeff","L2 coefficient for domain inner product",default=ssa_l2_coeff)
  ssa_h1_coeff = PISM.optionsReal("-ssa_h1_coeff","H1 coefficient for domain inner product",default=ssa_h1_coeff)
  tauc_guess_scale = PISM.optionsReal("-tauc_guess_scale","initial guess for tauc to be this factor of the true value",default=tauc_guess_scale)
  tauc_guess_const = PISM.optionsReal("-tauc_guess_const","initial guess for tauc to be this constant",default=tauc_guess_const)


config.set_string("inv_ssa_tauc_param",tauc_param_type)
config.set("inv_ssa_domain_l2_coeff",ssa_l2_coeff)
config.set("inv_ssa_domain_h1_coeff",ssa_h1_coeff)

tauc_param = tauc_params[tauc_param_type]

PISM.setVerbosityLevel(verbosity)
forward_problem = testi(Mx,My)
forward_problem.setup()

grid = forward_problem.grid

# Build the true yeild stress for test I
tauc_true = PISM.util.standardYieldStressVec(grid,name="tauc_true")
testi_tauc(grid,forward_problem.solver.ice,tauc_true)

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
forward_problem.F(PISMLocalVector(zeta_true),out=PISMLocalVector(u_true))

# Build the initial guess for tauc for the inversion.
tauc = PISM.util.standardYieldStressVec(grid)
if not tauc_guess_const is None:
  tauc.set(tauc_guess_const)
else:
  testi_tauc(grid,forward_problem.solver.ice,tauc)
  tauc.scale(tauc_guess_scale)

# Convert over to parameterized tauc
if config.get_string('inv_ssa_tauc_param')=='ident':
  zeta = tauc
else:
  zeta = PISM.IceModelVec2S();
  zeta.create(grid, "zeta", True, PISM.util.WIDE_STENCIL)
  with PISM.util.Access(nocomm=[tauc],comm=[zeta]):
    for (i,j) in grid.points():
      zeta[i,j] = tauc_param.fromTauc(tauc[i,j])

# Build the inversion algorithm, and set up its arguments.
if method == "ign":
  Solver = InvertSSAIGN
else:
  Solver = InvertSSANLCG

params=Solver.defaultParameters()
if method == "sd":
  params.steepest_descent = True
  params.ITER_MAX=10000
elif method =="ign":
  params.linearsolver.ITER_MAX=1000
  params.linearsolver.verbose = True
params.verbose   = True
params.deriv_eps = 0.

solver=Solver(forward_problem,params=params)  

# solver.params.linesearch.verbose=True
rms_error /= PISM.secpera # m/s
(zeta,u) = solver.solve(zeta,u_true,rms_error)

if config.get_string('inv_ssa_tauc_param')=='ident':
  tauc = zeta
else:
  with PISM.util.Access(nocomm=[zeta],comm=[tauc]):
    for (i,j) in grid.points():
      (tauc[i,j],dummy) = tauc_param.toTauc(zeta[i,j])

forward_problem.write(output_file)
tauc.write(output_file)
tauc_true.write(output_file)
u.set_name("_computed",0)
u.write(output_file)

tz = tozero.ToProcZero(grid)
tauc_a = tz.communicate(tauc)
tauc_true = tz.communicate(tauc_true)
if not tauc_a is None:
  from matplotlib import pyplot
  pyplot.plot(grid.y,tauc_a[:,Mx/2])
  pyplot.plot(grid.y,tauc_true[:,Mx/2])
  pyplot.draw()

context.com.barrier()
siple.reporting.endpause()
