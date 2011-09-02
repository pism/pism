from petsc4py import PETSc

class ToProcZero:
  
  def __init__(self, grid,dof=1,dim=2):
    self.grid = grid
    self.dof = dof
    self.dim = dim

    if dim != 2:
      raise NotImplementedError()
    
    if dof != 1:
      raise NotImplementedError()

    self.owns_da = False
    self.da = grid.da2
    
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
      rv = self.U0[...].reshape(self.da.sizes, order='f').copy()
    comm.barrier()

    return rv
