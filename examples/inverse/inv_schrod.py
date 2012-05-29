import sys, petsc4py
petsc4py.init(sys.argv)
from petsc4py import PETSc
import PISM
import math
import numpy as np

context = PISM.Context()
com = context.com

PISM.set_abort_on_sigint(True)

PISM.verbosityLevelFromOptions()
PISM.stop_on_version_option()

config = context.config

for o in PISM.OptionsGroup(com,"",sys.argv[0]):
  M = PISM.optionsInt("-M","problem size",default=30)
  eta = PISM.optionsReal("-eta","regularization paramter",default=300.)

x0 = 0.;
L  = math.pi;

# Build the grid.
grid = PISM.Context().newgrid()
PISM.util.init_shallow_grid(grid,L,L,M,M,PISM.NOT_PERIODIC)


cf = lambda x,y: 1.5+0.4*np.cos(2*x)*np.sin(1.5*y)
utruef = lambda x,y: PISM.PISMVector2(1+np.sin(x)*np.cos(y),np.sin(x)*np.cos(y))
ff = lambda x,y:  PISM.PISMVector2(2*np.sin(x)*np.cos(y),2*np.sin(x)*np.cos(y)) + PISM.PISMVector2(utruef(x,y).u*cf(x,y),utruef(x,y).v*cf(x,y))


stencil_width = 1

f = PISM.IceModelVec2V();
f.create(grid,"f",PISM.kHasGhosts,stencil_width)

utrue = PISM.IceModelVec2V();
utrue.create(grid,"utrue",PISM.kHasGhosts,stencil_width)

c = PISM.IceModelVec2S();
c.create(grid,"c",PISM.kHasGhosts,stencil_width)
c0 = PISM.IceModelVec2S();
c0.create(grid,"c",PISM.kHasGhosts,stencil_width)
c0.set(1)

with PISM.util.Access(comm=[c,f,utrue]):
  for (i,j) in grid.points():
    x=grid.x[i]; y=grid.y[j]
    f[i,j] = ff(x,y);
    c[i,j] = cf(x,y)
    utrue[i,j] = utruef(x,y)

dirichletIndices = PISM.PISM.IceModelVec2Int()
dirichletIndices.create(grid, "bc_mask", PISM.kHasGhosts, stencil_width);

dirichletValues = PISM.PISM.IceModelVec2V()
dirichletValues.create(grid, "bc_vals", PISM.kHasGhosts, stencil_width);

with PISM.util.Access(comm=[dirichletIndices,dirichletValues]):
  for (i,j) in grid.points():
    x=grid.x[i]; y=grid.y[j]
    dirichletValues[i,j] = utruef(x,y)
    if i==0 or j==0 or (i==M-1) or (j==M-1):
      dirichletIndices[i,j]=1;


invProblem = PISM.InvSchrodTikhonov(grid,f)
invProblem.set_c(c)
invProblem.setDirichletData(dirichletIndices,dirichletValues)

if not invProblem.solve():
  PISM.verbPrintf(1,grid.com,"Forward solve failed (%s)!\n" % invProblem.reasonDescription());
  exit(1)

u_obs = PISM.PISM.IceModelVec2V()
u_obs.create(grid, "true value", PISM.kHasGhosts, stencil_width);
u_obs.copy_from(invProblem.solution())


ip = PISM.InvSchrodTikhonovProblem(invProblem,c0,u_obs,eta)
solver = PISM.InvSchrodTikhonovSolver(grid.com,"tao_lmvm",ip)

if solver.solve():
  PISM.verbPrintf(1,grid.com,"Inverse solve success (%s)!\n" % solver.reasonDescription());
  u_i = ip.stateSolution();
  c_i = ip.designSolution();
  
  tozero1 = PISM.toproczero.ToProcZero(grid,dof=1,dim=2)
  ci_a = tozero1.communicate(c_i);
  c_a = tozero1.communicate(c);
  c0_a = tozero1.communicate(c0);
  
  
  tozero2 = PISM.toproczero.ToProcZero(grid,dof=2,dim=2)
  ui_a = tozero2.communicate(u_i);
  u_a = tozero2.communicate(u_obs);
  if context.rank == 0:
    import matplotlib.pyplot as pp
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
    pp.imshow(du_norm)
    pp.title("observation error")
    pp.colorbar()
    pp.draw()
  
  
    pause_time = PISM.optionsReal("-final_draw_pause","",default=10)
    import time
    time.sleep(pause_time)
else:
  PISM.verbPrintf(1,grid.com,"Inverse solve failure (%s)!\n" % solver.reasonDescription());
