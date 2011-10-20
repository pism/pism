import sys, petsc4py
petsc4py.init(sys.argv)
from petsc4py import PETSc
import numpy as np
from dmutils import msg
import tozero

import PISM

from pismssaforward import PISMSSAForwardProblem, InvertSSANLCG, InvertSSAIGN
from linalg_pism import PISMLocalVector

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

rms_error = 1. #m/a

def testi_tauc(grid,ice,tauc):
  standard_gravity = PISM.global_config().get("standard_gravity");
  slope = 0.001
  f = ice.rho*standard_gravity*H0_schoof*slope
  with PISM.util.Access(comm=tauc):
    for (i,j) in grid.points():
      y=grid.y[j]
      tauc[i,j] = f* (abs(y/L_schoof)**m_schoof)

class testi(PISMSSAForwardProblem):

  def initGrid(self,Mx,My):
    Ly = 3*L_schoof   # 300.0 km half-width (L=40.0km in Schoof's choice of variables)
    Lx = max(60.0e3, ((Mx - 1) / 2.) * (2.0 * Ly / (My - 1)) )
    PISM.util.init_shallow_grid(self.grid,Lx,Ly,Mx,My,PISM.XY_PERIODIC);

  def initPhysics(self):
    config = self.config
    self.basal = PISM.IceBasalResistancePlasticLaw(
         config.get("plastic_regularization") / PISM.secpera,
         config.get_flag("do_pseudo_plastic_till"),
         config.get("pseudo_plastic_q"),
         config.get("pseudo_plastic_uthreshold") / PISM.secpera);
 
    # irrelevant
    self.enthalpyconverter = PISM.EnthalpyConverter(config);

    self.ice = PISM.CustomGlenIce(self.grid.com, "", config,self.enthalpyconverter);
    self.ice.setHardness(B_schoof)


  def initSSACoefficients(self):
    solver = self.solver
    solver.allocateBCs()
    solver.bc_mask.set(0)
    solver.thickness.set(H0_schoof)
    solver.ice_mask.set(PISM.MASK_GROUNDED)

    # The finite difference code uses the following flag to treat 
    # the non-periodic grid correctly.
    self.config.set_flag("compute_surf_grad_inward_ssa", True);
    self.config.set("epsilon_ssafd", 0.0);  # don't use this lower bound
    
    testi_tauc(self.grid,self.ice,self.solver.tauc)
    # solver.tauc.set(0.)

    bc_mask = solver.bc_mask
    vel_bc  = self.solver.vel_bc
    surface = solver.surface
    bed     = solver.bed

    grid = self.grid
    vars = [surface,bed,vel_bc,bc_mask]
    with PISM.util.Access(comm=vars):
      for (i,j) in grid.points():
        x=grid.x[i]; y=grid.y[j]
        (bed_ij,junk,u,v) = PISM.exactI(m_schoof,x,y)
        bed[i,j] = bed_ij
        surface[i,j] = bed_ij + H0_schoof
    
        edge = ( (j == 0) or (j == grid.My - 1) ) or ( (i==0) or (i==grid.Mx-1) );
        if (edge):
          bc_mask[i,j] = 1;
          vel_bc[i,j].u = u;
          vel_bc[i,j].v = v;

## Main code starts here

context = PISM.Context()
config = context.config()
PISM.set_abort_on_sigint(True)

for o in PISM.OptionsGroup(context.com,"","Test J"):
  Mx = PISM.optionsInt("-Mx","Number of grid points in x-direction",default=Mx)
  My = PISM.optionsInt("-My","Number of grid points in y-direction",default=My)
  output_file = PISM.optionsString("-o","output file",default="invert_testi.nc")
  verbosity = PISM.optionsInt("-verbose","verbosity level",default=2)
  method = PISM.optionsList(context.com,"-inv_method","Inversion algorithm",["nlcg","ign","sd"],"ign")
  rms_error = PISM.optionsReal("-rms_error","RMS velocity error",default=rms_error)

PISM.setVerbosityLevel(verbosity)
forward_problem = testi(Mx,My)
grid = forward_problem.grid

tauc_true = PISM.util.standardYieldStressVec(grid,name="tauc_true")
testi_tauc(grid,forward_problem.ice,tauc_true)

u_true = PISM.util.standard2dVelocityVec(grid,name="_true")
forward_problem.F(PISMLocalVector(tauc_true),out=PISMLocalVector(u_true))

tauc = PISM.util.standardYieldStressVec(grid)
testi_tauc(grid,forward_problem.ice,tauc)
tauc.scale(0.01)
# tauc0 = 1e6 # Pa
# tauc.set(tauc0)

if method == "ign":
  solver = InvertSSAIGN(forward_problem)
else:
  solver = InvertSSANLCG(forward_problem)  
  if method == "sd":
    solver.params.steepest_descent = True
solver.params.verbose   = True
solver.params.deriv_eps = 0.
# solver.params.linesearch.verbose=True
rms_error /= PISM.secpera # m/s
(tauc,u) = solver.solve(tauc,u_true,rms_error)

tz = tozero.ToProcZero(grid)
tauc_a = tz.communicate(tauc)
if not tauc_a is None:
  from matplotlib import pyplot
  pyplot.imshow(tauc_a)
  pyplot.show()
context.com.barrier()

forward_problem.write(output_file)
tauc.write(output_file)
tauc_true.write(output_file)
u.set_name("_computed",0)
u.write(output_file)
