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
  N = PISM.optionsInt("-N","problem size",default=10)

x0 = 0.;
L  = math.pi;

# Build the grid.
grid = PISM.Context().newgrid()
PISM.util.init_shallow_grid(grid,L,L,N,N,PISM.NOT_PERIODIC)


cf = lambda x,y: 1.5+0.4*np.cos(4*x)*np.sin(3*y)
utruef = lambda x,y: PISM.PISMVector2(1+np.sin(x)*np.cos(y),np.sin(x)*np.cos(y))
ff = lambda x,y:  PISM.PISMVector2(2*np.sin(x)*np.cos(y),2*np.sin(x)*np.cos(y)) + PISM.PISMVector2(utruef(x,y).u*cf(x,y),utruef(x,y).v*cf(x,y))

# cf = lambda x,y: 1.
# utruef = lambda x,y: PISM.PISMVector2(np.sin(x),0)
# ff = lambda x,y:  PISM.PISMVector2(2*np.sin(x),0)



stencil_width = 1

f = PISM.IceModelVec2V();
f.create(grid,"f",PISM.kHasGhosts,stencil_width)

utrue = PISM.IceModelVec2V();
utrue.create(grid,"utrue",PISM.kHasGhosts,stencil_width)

c = PISM.IceModelVec2S();
c.create(grid,"c",PISM.kHasGhosts,stencil_width)

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
    if i==0 or j==0 or (i==N-1) or (j==N-1):
      dirichletIndices[i,j]=1;


invProblem = PISM.InvSchrodTikhonov(grid,f)
invProblem.set_c(c)
invProblem.setDirichletData(dirichletIndices,dirichletValues)

if invProblem.solve():
  PISM.verbPrintf(1,grid.com,"Success (%s)!\n" % invProblem.reasonDescription());
  u = invProblem.solution();
  tozero = PISM.toproczero.ToProcZero(grid,dof=2,dim=2)
  u_a = tozero.communicate(u);
  ut_a = tozero.communicate(utrue);
  if context.rank == 0:
    import matplotlib.pyplot as pp
    U = u_a[0,:,:]
    V = u_a[1,:,:]
    Ut = ut_a[0,:,:]
    Vt = ut_a[1,:,:]
    
    pp.imshow(U-Ut)
    pp.colorbar()
    pp.draw()

    pp.figure()
    pp.imshow(V-Vt)
    pp.colorbar()
    pp.draw()

    pause_time = PISM.optionsReal("-draw_pause","",default=1)
    import time
    time.sleep(pause_time)

    # pp.show()
else:
  PISM.verbPrintf(1,grid.com,"Failure (%s)!\n" % invProblem.reasonDescription());
