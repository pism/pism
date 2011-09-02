import PISM
from siple.linalg import AbstractVector
from petsc4py import PETSc

class PISMLocalVector(AbstractVector):
  """Implements the invtools.linalg.AbstractVector interface for a dolfin generic vector."""
  def __init__(self,u):
    if isinstance(u,PISM.IceModelVec2S):
      self._core = u
      self.dof = 1
    elif isinstance(u,PISM.IceModelVec2V):
      self._core = u
      self.dof = 2
    else:
      raise ValueError("An PISMLocalVector wraps PISM IceModelVec2S or IceModelVec2V: found a %s" % u)
    self.grid = u.get_grid()
  
  def set(self,rhs):
    self._core.copy_from(rhs._core)
  
  def acc(self,rhs):
    self._core.add(1.,rhs._core)
  
  def scale(self,t):
    self._core.scale(t)
  
  def axpy(self,t,v):
    self._core.add(t,v._core)
  
  def copy(self):
    c = self.vector_like()
    self._core.copy_to(c._core)
    return c
  
  def vector_like(self):
    if self.dof == 1:
      c = PISM.IceModelVec2S()
    else:
      c = PISM.IceModelVec2V()    
    c.create(self.grid,"",True,2)
    return PISMLocalVector(c)
    
  def zero_like(self):
    z = self.vector_like()
    z._core.set(0.)
    return z
  
  def dim(self):
    grid = self._core.get_grid()
    return self.dof*grid.Mx*grid.My
  
  def core(self):
    return self._core

  def norm(self,name):
    if name == 'linf':
      return self._core.norm(PETSc.NormType.NORM_INFINITY)
    if name == 'l2':
      return self._core.norm(PETSc.NormType.NORM_2)
    if name == 'l1':
      return self._core.norm(PETSc.NormType.NORM_1)
    
    raise ValueError()
