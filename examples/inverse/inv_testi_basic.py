import sys, petsc4py
petsc4py.init(sys.argv)
from petsc4py import PETSc
import PISM
import PISM.invert_ssa
import math
import numpy as np
import matplotlib.pyplot as pp

context = PISM.Context()
com = context.com
config = context.config

PISM.set_abort_on_sigint(True)

PISM.verbosityLevelFromOptions()
PISM.stop_on_version_option()

tauc_guess_scale = 0.3
tauc_guess_const = None

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


for o in PISM.OptionsGroup(com,"",sys.argv[0]):
  Mx = PISM.optionsInt("-Mx","problem size",default=11)
  My = PISM.optionsInt("-My","problem size",default=61)
  eta = PISM.optionsReal("-eta","regularization paramter",default=5e4)
  cL2 = PISM.optionsReal("-inv_ssa_cL2","regularization paramter",default=1)
  cH1 = PISM.optionsReal("-inv_ssa_cH1","regularization paramter",default=0)
  tauc_guess_scale = PISM.optionsReal("-tauc_guess_scale","initial guess for tauc to be this factor of the true value",default=tauc_guess_scale)
  tauc_guess_const = PISM.optionsReal("-tauc_guess_const","initial guess for tauc to be this constant",default=tauc_guess_const)

  p = PISM.optionsReal("-p","glen exponent",default=4./3.)
  q = PISM.optionsReal("-q","pseudo yeild stress exponent",default=0)
  hasFixedDesignLocs = PISM.optionsFlag("-use_fixed_design_locs","test having fixed design variables",default=False)
  hasMisfitWeight = PISM.optionsFlag("-use_misfit_weight","test misfit weight paramter",default=False)
  isDimensionless = PISM.optionsInt("-is_dimensionless","use dimensionless units",default=1)
  isDimensionless = bool(isDimensionless)
x0 = 0.;

L_schoof = 40e3;      # meters
Ly_dim = 3*L_schoof;      # meters
Lx_dim = max(60.0e3, ((Mx - 1) / 2.) * (2.0 * Ly_dim / (My - 1)) )
# # aspect_schoof = 0.05; # (pure)
# # H0_schoof = aspect_schoof * L_schoof; 
#                       # = 2000 m THICKNESS
H0_dim = 2000
B_dim = 3.7e8;     # Pa s^{1/3}; hardness 

slope = 0.001
standard_gravity = config.get("standard_gravity");
ice_density = config.get("ice_density");
f0_dim = ice_density * standard_gravity * H0_dim * slope

velocity_scale_dim   = (f0_dim/B_dim)**(3.)*(L_schoof/H0_dim)**(3.)*L_schoof
time_scale_dim       = L_schoof / velocity_scale_dim
strainrate_scale_dim = 1./ time_scale_dim
viscosity_scale_dim  = B_dim*(strainrate_scale_dim**(-2./3.))
nuH_scale_dim        = viscosity_scale_dim * H0_dim

eps_velocity = (config.get("plastic_regularization") / PISM.secpera)/velocity_scale_dim
pp_threshold = (config.get("pseudo_plastic_uthreshold") / PISM.secpera)/velocity_scale_dim
schoofLen      = config.get("Schoof_regularizing_length", "km", "m"); #// convert to meters
schoofVel      = config.get("Schoof_regularizing_velocity", "m/year", "m/s"); #// convert to m/s
eps_strainrate = abs(schoofVel/schoofLen)/strainrate_scale_dim
eps_nuH        = config.get("epsilon_ssa")/nuH_scale_dim

print "Epsilons: vel %g strainrate %g nuH %g threshold %g" % (eps_velocity,eps_strainrate,eps_nuH,pp_threshold)

if isDimensionless:
  Ly  = Ly_dim/L_schoof
  Lx =  Lx_dim/L_schoof
  L_schoof = 1.
  length_scale = 1
  H  = 1.
  f0 = 1.
  B  = 1.
else:
  length_scale = L_schoof
  Ly  = Ly_dim
  Lx  = Lx_dim
  H  = H0_dim
  f0 = f0_dim
  B  = B_dim

stress_scale  = f0
area_scale    = Ly*Lx
depth_scale   = H

velocity_scale   = (f0/B)**(3.)*(length_scale/depth_scale)**(3.)*length_scale
time_scale       = length_scale / velocity_scale 
strainrate_scale = 1./ time_scale
viscosity_scale  = B*(strainrate_scale**(-2./3.))
nuH_scale        = viscosity_scale * depth_scale

cL2 /= area_scale 
# cH1 /= (stress_scale**2)

config.set("inv_ssa_cL2",cL2)
config.set("inv_ssa_cH1",cH1)

config.set("tauc_param_trunc_tauc0",.1*stress_scale)
config.set("tauc_param_tauc_eps",.001*stress_scale)
config.set("tauc_param_tauc_scale",stress_scale)

m=2

one_m_per_a = (1./PISM.secpera)*(velocity_scale/velocity_scale_dim)
config.set("inv_ssa_velocity_scale",one_m_per_a) # m/a

config.set_string("inv_ssa_tauc_param","ident")

# Build the grid.
grid = PISM.Context().newgrid()
PISM.util.init_shallow_grid(grid,Lx,Ly,Mx,My,PISM.X_PERIODIC)
grid.printInfo(1)

# Coefficient data
cf = lambda x,y:  f0*(abs(y/L_schoof))**m
ff = lambda x,y:  PISM.PISMVector2(f0,0)

# one neighbouring ghost only for all variables with ghosts.
stencil_width = 1

# Main model parameters: f,c and the associated utrue
f = PISM.IceModelVec2V();
f.create(grid,"f",PISM.kHasGhosts,stencil_width)
c = PISM.IceModelVec2S();
c.create(grid,"c",PISM.kHasGhosts,stencil_width)
# Initial guess for c.
c0 = PISM.IceModelVec2S();
c0.create(grid,"c",PISM.kHasGhosts,stencil_width)
with PISM.util.Access(comm=[c,f,c0]):
  for (i,j) in grid.points():
    x=grid.x[i]; y=grid.y[j]
    f[i,j] = ff(x,y);
    c[i,j] = cf(x,y)

if tauc_guess_const is not None:
  c0.set(tauc_guess_const)
else:
  c0.copy_from(c)
  c0.scale(tauc_guess_scale)

# Convert c/c0 to zeta/zeta0
tauc_param = PISM.invert_ssa.tauc_param_factory.create(config)
zeta = PISM.IceModelVec2S();
zeta.create(grid,"zeta",PISM.kHasGhosts,stencil_width)
zeta0 = PISM.IceModelVec2S();
zeta0.create(grid,"zeta0",PISM.kHasGhosts,stencil_width)
tauc_param.convertFromTauc(c,zeta)
tauc_param.convertFromTauc(c0,zeta0)

# Create Dirichlet data and apply it
dirichletIndices = PISM.PISM.IceModelVec2Int()
dirichletIndices.create(grid, "bc_mask", PISM.kHasGhosts, stencil_width);
dirichletValues = PISM.PISM.IceModelVec2V()
dirichletValues.create(grid, "bc_vals", PISM.kHasGhosts, stencil_width);
dirichletValues.set(0)
with PISM.util.Access(comm=[dirichletIndices,dirichletValues]):
  for (i,j) in grid.points():
    x=grid.x[i]; y=grid.y[j]
    if  j==0 or j==(My-1):
      dirichletIndices[i,j]=1;
invProblem = PISM.InvSSABasicTikhonov(grid,f,p,q,tauc_param)
invProblem.setDirichletData(dirichletIndices,dirichletValues)
invProblem.setBH(B,H)


invProblem.setEpsilon(eps_velocity*velocity_scale,eps_strainrate*strainrate_scale,eps_nuH*nuH_scale)
invProblem.setPseudoPlasticThreshold(pp_threshold*velocity_scale)

# If testing zeta_fixed mask (locations where zeta is not allowed to change)
# set it here.
fixedDesignLocs = PISM.PISM.IceModelVec2Int()
if hasFixedDesignLocs:
  fixedDesignLocs.create(grid, "zeta_fixed_mask", PISM.kHasGhosts, stencil_width);
  fixedDesignLocs.set(0);
  with PISM.util.Access(comm=[fixedDesignLocs]):
    for (i,j) in grid.points():
      if(j>My*2./3.):
        fixedDesignLocs[i,j]=1;
  invProblem.setFixedDesignLocations(fixedDesignLocs)  

# If testing the misfit weight function set it here
misfitWeight = PISM.PISM.IceModelVec2S()
if hasMisfitWeight:
  misfitWeight.create(grid, "vel_misfit_weight", PISM.kNoGhosts, stencil_width);
  misfitWeight.set(1);
  with PISM.util.Access(nocomm=[misfitWeight]):
    for (i,j) in grid.points():
      if(j<My/2):
        misfitWeight[i,j]=0;
  invProblem.setObservationWeights(misfitWeight)

# Solve the forward problem for the true value of c/zeta and obtain an
# observed solution 
invProblem.setZeta(zeta)
if not invProblem.solve():
  PISM.verbPrintf(1,grid.com,"Forward solve failed (%s)!\n" % invProblem.reasonDescription());
  exit(1)
u_obs = PISM.PISM.IceModelVec2V()
u_obs.create(grid, "true value", PISM.kHasGhosts, stencil_width);
u_obs.copy_from(invProblem.solution())

invProblem.setZeta(zeta0)
if not invProblem.solve():
  PISM.verbPrintf(1,grid.com,"Forward solve failed (%s)!\n" % invProblem.reasonDescription());
  exit(1)
u0 = PISM.PISM.IceModelVec2V()
u0.create(grid, "initial value", PISM.kHasGhosts, stencil_width);
u0.copy_from(invProblem.solution())
du = PISM.PISM.IceModelVec2V()
du.create(grid, "du", PISM.kHasGhosts, stencil_width);
du.copy_from(u0)
du.add(-1,u_obs)


designFunc =PISM.H1NormFunctional2S(grid,cL2,cH1)
stateFunc = PISM.MeanSquareFunctional2V(grid)
stateFunc.normalize( grid.config.get("inv_ssa_velocity_scale")  )
print "cL2, cH1", cL2,cH1
grid.printInfo(1)
print "design functional value is: %g" % designFunc.valueAt(zeta0)
print "state functional value is: %g; eta=%g" % (stateFunc.valueAt(du),eta)

maxU = 0 
maxZeta = 0
maxTauc = 0
with PISM.util.Access(nocomm=[u_obs,zeta,c]):
  for (i,j) in grid.points():
    U=u_obs[i,j].magnitude()
    if U>maxU:
      maxU = U
    Zeta=zeta[i,j]
    if Zeta>maxZeta:
      maxZeta = Zeta
    TAUC=c[i,j]
    if TAUC>maxTauc:
      maxTauc = TAUC
print "Maximum velocity is : %g m/s = %g m/a.  Relative: %g" % (maxU,maxU*PISM.secpera,maxU/velocity_scale)
print "Maximum zeta is : %g Relative: %g" % (maxZeta,maxZeta)
print "Maximum tauc is : %g Relative: %g" % (maxTauc,maxTauc/stress_scale)

# tozero = PISM.toproczero.ToProcZero(grid,dof=2,dim=2)
# import schoof
# import matplotlib.pyplot as pp
# standard_gravity = grid.config.get("standard_gravity");
# ice_density = grid.config.get("ice_density");
# s = schoof.SchoofSSA1dExact(L_schoof,m,B=B,h=H,f=f0)
# u_schoof_a = [s.eval(y) for y in grid.y]
# u_true_a = tozero.communicate(u_obs)
# pp.plot(grid.y,u_true_a[0,:,Mx/2],grid.y,u_schoof_a)
# pp.show()
# quit()

# Build the inverse solver.  This first step feels redundant.
ip = PISM.InvSSABasicTikhonovProblem(invProblem,zeta0,u_obs,eta)
l=Listener()
ip.addListener(l)
solver = PISM.InvSSABasicTikhonovSolver(grid.com,"tao_lmvm",ip)

# Try solving
if not solver.solve():
  PISM.verbPrintf(1,grid.com,"Inverse solve FAILED (%s)!\n" % solver.reasonDescription());
else:
  PISM.verbPrintf(1,grid.com,"Inverse solve success (%s)!\n" % solver.reasonDescription());

u_i    = ip.stateSolution();
zeta_i = ip.designSolution();

c_i = PISM.IceModelVec2S();
c_i.create(grid,"c_i",PISM.kHasGhosts,stencil_width)
tauc_param.convertToTauc(zeta_i,c_i);

du = PISM.IceModelVec2V()
du.create(grid, "du", PISM.kHasGhosts, stencil_width);
du.copy_from(u_i)
du.add(-1,u_obs)
if hasMisfitWeight:
  misfit_functional = PISM.MeanSquareFunctional2V(grid,misfitWeight);
else:
  misfit_functional = PISM.MeanSquareFunctional2V(grid);
misfit_functional.normalize( config.get("inv_ssa_velocity_scale")  )
misfit = math.sqrt(misfit_functional.valueAt(du))
misfit_norm = math.sqrt(misfit_functional.valueAt(u_i))  
PISM.verbPrintf(1,grid.com,"RMS Misfit: %g (relative %g)\n",misfit,misfit/misfit_norm)


tozero1 = PISM.toproczero.ToProcZero(grid,dof=1,dim=2)
ci_a = tozero1.communicate(c_i);
c_a = tozero1.communicate(c);
c0_a = tozero1.communicate(c0);


tozero2 = PISM.toproczero.ToProcZero(grid,dof=2,dim=2)
ui_a = tozero2.communicate(u_i);
u_a = tozero2.communicate(u_obs);
if context.rank == 0:
  cmin=np.min(c_a); cmax = np.max(c_a)

  pp.imshow(ci_a,vmin=cmin,vmax=cmax)
  pp.colorbar()
  pp.title("c from inversion")
  pp.draw()
  
  pp.figure()
  pp.imshow(c_a,vmin=cmin,vmax=cmax)
  pp.title("c true value")
  pp.colorbar()
  pp.draw()

  du_a = ui_a-u_a
  du_norm = np.sqrt(du_a[0,:,:]*du_a[0,:,:]+du_a[1,:,:]*du_a[1,:,:])
  pp.figure()
  pp.subplot(1,2,1)
  pp.imshow(u_a[0,:,:])
  pp.title("u_x")
  pp.colorbar()

  pp.subplot(1,2,2)
  pp.imshow(du_norm)
  pp.title("u-u_obs")
  pp.colorbar()

  pp.draw()

  pp.figure()
  pp.plot(grid.y,ci_a[:,Mx/2])
  pp.plot(grid.y,c_a[:,Mx/2])
  pp.ion()
  pp.show()

  pause_time = PISM.optionsReal("-final_draw_pause","",default=10)
  import time
  time.sleep(pause_time)
