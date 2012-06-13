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


for o in PISM.OptionsGroup(com,"",sys.argv[0]):
  M = PISM.optionsInt("-M","problem size",default=30)
  eta = PISM.optionsReal("-eta","regularization paramter",default=5e4)
  p = PISM.optionsReal("-p","glen exponent",default=2)
  q = PISM.optionsReal("-q","pseudo yeild stress exponent",default=1)
  hasFixedDesignLocs = PISM.optionsFlag("-fixed_design_locs","test having fixed design variables",default=False)
  hasFixedDesignLocs = PISM.optionsFlag("-fixed_design_locs","test having fixed design variables",default=False)
  hasMisfitWeight = PISM.optionsFlag("-use_misfit_weight","test misfit weight paramter",default=False)
  inv_ssa_cL2 = PISM.optionsReal("-inv_ssa_cL2","L2 coefficient for domain inner product",default=1)
  inv_ssa_cH1 = PISM.optionsReal("-inv_ssa_cH1","H1 coefficient for domain inner product",default=0)

x0 = 0.;
L = math.pi;
m=4

config.set_string("inv_ssa_tauc_param","ident")
config.set('tauc_param_tauc_scale',1)
config.set('inv_ssa_velocity_scale',1)
config.set('inv_ssa_cL2',inv_ssa_cL2)
config.set('inv_ssa_cH1',inv_ssa_cH1)

# Build the grid.
grid = PISM.Context().newgrid()
PISM.util.init_shallow_grid(grid,L,L,M,M,PISM.NOT_PERIODIC)

# Coefficient data
cf = lambda x,y: 1.5+0.4*np.cos(2*x)*np.sin(1.5*y)
ff = lambda x,y:  PISM.PISMVector2(2*np.sin(x)*np.cos(y),2*np.sin(x)*np.cos(y)) 

# one neighbouring ghost only for all variables with ghosts.
stencil_width = 1

# Main model parameters: f,c and the associated utrue
f = PISM.IceModelVec2V();
f.create(grid,"f",PISM.kHasGhosts,stencil_width)
c = PISM.IceModelVec2S();
c.create(grid,"c",PISM.kHasGhosts,stencil_width)
with PISM.util.Access(comm=[c,f]):
  for (i,j) in grid.points():
    x=grid.x[i]; y=grid.y[j]
    f[i,j] = ff(x,y);
    c[i,j] = cf(x,y)

# Initial guess for c.
c0 = PISM.IceModelVec2S();
c0.create(grid,"c",PISM.kHasGhosts,stencil_width)
c0.set(1)

# Convert c/c0 to zeta/zeta0
tauc_param = PISM.invert_ssa.tauc_param_factory.create(config)
zeta0 = PISM.IceModelVec2S();
zeta0.create(grid,"zeta0",PISM.kHasGhosts,stencil_width)
zeta = PISM.IceModelVec2S();
zeta.create(grid,"zeta",PISM.kHasGhosts,stencil_width)
with PISM.util.Access(comm=[zeta0,zeta],nocomm=[c0,c]):
  for (i,j) in grid.points():
    zeta0[i,j] = tauc_param.fromTauc(c0[i,j])
    zeta[i,j] = tauc_param.fromTauc(c[i,j])


# Create Dirichlet data and apply it
dirichletIndices = PISM.PISM.IceModelVec2Int()
dirichletIndices.create(grid, "bc_mask", PISM.kHasGhosts, stencil_width);
dirichletValues = PISM.PISM.IceModelVec2V()
dirichletValues.create(grid, "bc_vals", PISM.kHasGhosts, stencil_width);
dirichletValues.set(0)
with PISM.util.Access(comm=[dirichletIndices,dirichletValues]):
  for (i,j) in grid.points():
    x=grid.x[i]; y=grid.y[j]
    if i==0 or i==(M-1) or j==0 or j==(M-1):
      dirichletIndices[i,j]=1;
invProblem = PISM.InvSSABasicTikhonov(grid,f,p,q,tauc_param)
invProblem.setZeta(zeta0)
invProblem.setDirichletData(dirichletIndices,dirichletValues)

# If testing zeta_fixed mask (locations where zeta is not allowed to change)
# set it here.
fixedDesignLocs = PISM.PISM.IceModelVec2Int()
if hasFixedDesignLocs:
  fixedDesignLocs.create(grid, "zeta_fixed_mask", PISM.kHasGhosts, stencil_width);
  fixedDesignLocs.set(0);
  with PISM.util.Access(comm=[fixedDesignLocs]):
    for (i,j) in grid.points():
      if(i>M/2) and (j>M/2):
        fixedDesignLocs[i,j]=1;
  invProblem.setFixedDesignLocations(fixedDesignLocs)  

# If testing the misfit weight function set it here
misfitWeight = PISM.PISM.IceModelVec2S()
if hasMisfitWeight:
  misfitWeight.create(grid, "vel_misfit_weight", PISM.kNoGhosts, stencil_width);
  misfitWeight.set(1);
  with PISM.util.Access(nocomm=[misfitWeight]):
    for (i,j) in grid.points():
      if(j<M/2):
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

# Build the inverse solver.  This first step feels redundant.
ip = PISM.InvSSABasicTikhonovProblem(invProblem,zeta0,u_obs,eta)
solver = PISM.InvSSABasicTikhonovSolver(grid.com,"tao_lmvm",ip)

# Try solving
if solver.solve():
  PISM.verbPrintf(1,grid.com,"Inverse solve success (%s)!\n" % solver.reasonDescription());
  u_i = ip.stateSolution();
  zeta_i = ip.designSolution();
  
  du = PISM.IceModelVec2V()
  du.create(grid, "du", PISM.kHasGhosts, stencil_width);
  du.copy_from(u_i)
  du.add(-1,u_obs)
  if hasMisfitWeight:
    misfit_functional = PISM.MeanSquareFunctional2V(grid,misfitWeight);
  else:
    misfit_functional = PISM.MeanSquareFunctional2V(grid);
  misfit_functional.normalize(1)
  misfit = math.sqrt(misfit_functional.valueAt(du))
  misfit_norm = math.sqrt(misfit_functional.valueAt(u_i))  
  PISM.verbPrintf(1,grid.com,"RMS Misfit: %g (relative %g)\n",misfit,misfit/misfit_norm)

  c_i = PISM.IceModelVec2S();
  c_i.create(grid,"c_i",PISM.kHasGhosts,stencil_width)
  with PISM.util.Access(comm=c_i,nocomm=zeta_i):
    for (i,j) in grid.points():
      (c_i[i,j],dummy) = tauc_param.toTauc(zeta_i[i,j])
  
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
    pp.imshow(u_a[0,:,:])
    pp.title("observation error")
    pp.colorbar()
    pp.draw()
  
    pause_time = PISM.optionsReal("-final_draw_pause","",default=10)
    import time
    time.sleep(pause_time)
else:
  PISM.verbPrintf(1,grid.com,"Inverse solve failure (%s)!\n" % solver.reasonDescription());
