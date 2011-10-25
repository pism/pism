from petsc4py import PETSc

class ToProcZero:
  
  def __init__(self, grid,dof=1,dim=2):
    self.grid = grid
    self.dof = dof
    self.dim = dim

    if dim != 2:
      raise NotImplementedError()

    if dof == 1:
      self.owns_da = False
      self.da = grid.da2
    elif dof == 2:
      da = grid.da2
      self.da = PETSc.DA().create(dim=da.dim,dof=2,sizes=da.sizes,
      proc_sizes=da.proc_sizes,periodic_type=da.periodic_type,
      stencil_width =da.stencil_width,comm=grid.com)
      self.owns_da = True
    else:
      raise NotImplementedError()

    
    self.tmp_U         = self.da.createGlobalVector()
    self.tmp_U_natural = self.da.createNaturalVector()
    self.scatter, self.U0 = PETSc.Scatter.toZero(self.tmp_U_natural)

  def __del__(self):
    self.tmp_U.destroy()
    self.tmp_U_natural.destroy()
    if self.owns_da:
      self.da.destroy()

  def communicate(self,u):    
    comm = self.grid.da2.getComm()
    rank = comm.getRank()

    u.copy_to(self.tmp_U)
    self.da.globalToNatural(self.tmp_U,self.tmp_U_natural)
    self.scatter.scatter(self.tmp_U_natural, self.U0, False, PETSc.Scatter.Mode.FORWARD)

    rv = None
    if rank == 0:
      if self.dof == 1:
        rv = self.U0[...].reshape(self.da.sizes, order='f').copy()
      else:
        s=self.da.sizes
        rv = self.U0[...].reshape((2,s[0],s[1]), order='f').copy()
      
    comm.barrier()

    return rv
