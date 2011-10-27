import PISM

from siple.gradient.forward import NonlinearForwardProblem
from siple.gradient.nonlinear import InvertNLCG, InvertIGN
from siple.reporting import msg, std_pause


from linalg_pism import PISMLocalVector
from math import sqrt
import numpy as np

def pism_print_logger(message,severity):
  verb = severity
  com = PISM.Context().com
  PISM.verbPrintf(verb,com, "%s\n" % message)

def pism_pause(message_in=None,message_out=None):
  import sys, os
  fd = sys.stdin.fileno()
  if os.isatty(fd):
    return std_pause(message_in,message_out)
  if not message_in is None:
    PISM.verbPrintf(1,com,message_in+"\n")
  import sys
  ch = sys.stdin.read(1)
  if not message_out is None:
    PISM.verbPrintf(1,com,message_out+"\n")

def randVectorS(grid,scale):
  rv = PISM.IceModelVec2S();
  rv.create(grid, 'rand vec', True, PISM.util.WIDE_STENCIL)
  shape=(grid.xm,grid.ym)
  r = np.random.normal(scale=scale,size=shape)
  with PISM.util.Access(comm=rv):
    for (i,j) in grid.points():
      rv[i,j] = r[i,j]
  return rv

def randVectorV(grid,scale):
  rv = PISM.IceModelVec2V();
  rv.create(grid, 'rand vec', True, PISM.util.WIDE_STENCIL)
  shape=(grid.xm,grid.ym)
  r_u = np.random.normal(scale=scale,size=shape)
  r_v = np.random.normal(scale=scale,size=shape)
  with PISM.util.Access(comm=rv):
    for (i,j) in grid.points():
      rv[i,j].u = r_u[i,j]
      rv[i,j].v = r_v[i,j]
  return rv

tauc_params = {"ident":PISM.cvar.g_InvTaucParamIdent, 
               "square":PISM.cvar.g_InvTaucParamSquare,
               "exp":PISM.cvar.g_InvTaucParamExp }

class SSAForwardSolver(PISM.ssa.SSASolver):
  def __init__(self,grid,config):
    PISM.ssa.SSASolver.__init__(self,grid,config)
    tauc_param_type = config.get_string("inv_ssa_tauc_param")
    self.tauc_param = tauc_params[tauc_param_type]
    self.range_l2_weight=None

  def allocateCoeffs(self,using_l2_weight=False,**kwargs):
    PISM.ssa.SSASolver.allocateCoeffs(self,**kwargs)

    if using_l2_weight:      
      self.range_l2_weight = PISM.IceModelVec2S();
      self.range_l2_weight.create(self.grid, 'range_l2_weight', True, PISM.util.WIDE_STENCIL)
      self.range_l2_weight.set_attrs("diagnostic", "range l2 weight", "", "");
    
  
  def init_vars(self):
    pismVars = PISM.PISMVars()
    for var in [self.thickness,self.bed,self.tauc,
                self.enthalpy,self.ice_mask]:
        pismVars.add(var)
    
    if not self.surface is None:
      pismVars.add(self.surface)
    if not self.drivingstress is None:
      pismVars.add(self.drivingstress,'ssa_driving_stress')
    
    
    if not self.range_l2_weight is None:
      pismVars.add(self.range_l2_weight)
      
    
    # The SSA instance will not keep a reference to pismVars; it only uses it to extract
    # its desired variables.  So it is safe to pass it pismVars and then let pismVars
    # go out of scope at the end of this method.
    self.ssa.init(pismVars)

    if not self.vel_bc is None:
      self.ssa.set_boundary_conditions(self.bc_mask,self.vel_bc)

    self.ssa.setup_vars()

    self.ssa_init = True
    
  def buildSSA(self):
    self.ssa = PISM.InvSSAForwardProblem(self.grid,self.basal,self.ice,self.enthalpyconverter,self.tauc_param,self.config)

  def write(self,filename):
    PISM.ssa.SSASolver.write(self,filename)
    if not self.range_l2_weight is None:
      self.range_l2_weight.write(filename)

WIDE_STENCIL = 2

class InvSSARun(PISM.ssa.SSARun):
  def setup(self):
    PISM.ssa.SSARun.setup(self)

    # FIXME: This is a lousy name....
    self.solver.init_vars()

  def _constructSSA(self):
    return SSAForwardSolver(self.grid,self.config)

class SSAForwardProblem(NonlinearForwardProblem):
  
  def __init__(self,ssarun):
    self.ssarun = ssarun
    self.solver = ssarun.solver
    self.grid = ssarun.grid

  def F(self, x,out=None,guess=None):
    """
    Returns the value of the forward problem at x.

    Nonlinear problems often make use of an initial guess; this can be provided in 'guess'.

    Storage in 'out', if given, is used for the return value.
    """
    # if not guess is None:
    #   self.solver.set_initial_velocity_guess(guess)
    if out is None:
      out = self.rangeVector()
    self.solver.ssa.set_tauc(x.core())
    self.solver.ssa.solveF(out.core())
    return out

  def T(self,d,out=None):
    """
    Returns the value of the linearization, T, of F, at the point x specified previously in linearizeAt, 
    in the direction d.

    Storage in 'out', if given, is used for the return value.
    """
    if out is None:
      out = self.rangeVector()
    self.solver.ssa.solveT(d.core(),out.core())
    return out

  def TStar(self,r,out=None):
    """
    Let T be the linearization of F at the point x (at the point x specified previously in linearizeAt).  
    Its adjoint is T*.  This method returns the value of T* in the direction r.

    Storage in 'out', if given, is used for the return value.
    """
    if out is None:
      out = self.domainVector()
    self.solver.ssa.solveTStar(r.core(),out.core())
    return out

  def linearizeAt(self,x,guess=None):
    """
    Instructs the class that subsequent calls to T and TStar will be conducted for the given value of x.

    Nonlinear problems often make use of an initial guess; this can be provided in 'guess'.
    """
    self.solver.ssa.set_tauc(x.core())

  def evalFandLinearize(self,x,out=None,guess=None):
    """
    Computes the value of F(x) and locks in a linearization.  Sometimes there are efficiencies that
    can be acheived this way.

    Default implementation simply calls F, then linearizeAt.
    """
    if out is None:
      out = self.rangeVector()
    self.linearizeAt(x)
    self.solver.ssa.solveF(out.core())
    return out
  
  def rangeIP(self,a,b):
    """
    Computes the inner product of two vectors in the range.
    """
    return self.solver.ssa.rangeIP(a.core(),b.core())

  def domainIP(self,a,b):
    """
    Computes the inner product of two vectors in the domain.
    """
    return self.solver.ssa.domainIP(a.core(),b.core())

  def rangeVector(self):
    """Constructs a brand new vector from the range vector space"""
    v = PISM.IceModelVec2V()
    v.create(self.grid,"",True,WIDE_STENCIL)
    return PISMLocalVector(v)

  def domainVector(self):
    """Constructs a brand new vector from the domain vector space"""
    v = PISM.IceModelVec2S()
    v.create(self.grid,"",True,WIDE_STENCIL)
    return PISMLocalVector(v)

class InvertSSANLCG(InvertNLCG):
  """
  Inversion of the map

    F: gamma |-> u

  where u is the solution of the PDE

    -Laplacian u + gamma u = f

  F is a map from L^2 to L^2.
  """

  @staticmethod
  def defaultParameters():
    params = InvertNLCG.defaultParameters()
    return params

  def __init__(self,forward_problem,params=None):
    InvertNLCG.__init__(self,params)
    self.forward_problem = forward_problem

  def forwardProblem(self):
    return self.forward_problem

  def stopConditionMet(self,count,x,Fx,y,r):
    """
    Determines if minimization should be halted (based, e.g. on a Morozov discrepancy principle)

    In: count: current iteration count
        x:     point in domain of potential minimizer.
        Fx:    value of nonlinear function at x
        r:     current residual, i.e. y-F(x)    
    """

    J = sqrt(abs(self.forward_problem.rangeIP(r,r)));

    msg('(%d) J=%g goal=%g',count,J,self.Jgoal)

    if( J < self.Jgoal ):
      msg('Stop condition met')
      return True
    return False

  def initialize(self,x,y,deltaLInf):
    """
    Hook called at the start of solve.  This gives the class a chance to massage the input.
    For example, x and y might be dolfin (finite element) Functions; this method should return 
    the associated dolfin GenericVectors.

    The remaining arguments are passed directly from solve, and can be used for determining the
    final stopping criterion.

    Returns dolfin vectors corresponding to the initial value of x and the desired value of y=F(x).    
    """
    xv = PISMLocalVector(x)
    yv = PISMLocalVector(y)

    self.Jgoal = self.params.mu*deltaLInf

    return (xv,yv)

  def finalize(self,x,y):
    """
    Hook called at the end of 'solve'.  Gives the chance to massage the return values.
    """
    tauc = x.core()
    u = y.core()
    return (tauc,u)

  def solve(self,x,y,deltaLInf):
    """
    Solve the ill posed problem F(x)=y where y is know to an L infinity error deltaLInf
    """
    return InvertNLCG.solve(self,x,y,deltaLInf)


class InvertSSAIGN(InvertIGN):
  """
  Inversion of the map

    F: gamma |-> u

  where u is the solution of the PDE

    -Laplacian u + gamma u = f

  F is a map from L^2 to L^2.
  """

  @staticmethod
  def defaultParameters():
    params = InvertIGN.defaultParameters()
    return params

  def __init__(self,forward_problem,params=None):
    InvertIGN.__init__(self,params)
    self.forward_problem = forward_problem

  def forwardProblem(self):
    return self.forward_problem

  def temper_d(self, x,d,y,r):
    dnorm = d.norm('linf');  xnorm = x.norm('linf')
    if dnorm > 2*xnorm:
      msg('wild change predicted by linear step. scaling')
      d.scale(2*xnorm/dnorm)

  def initialize(self,x,y,deltaLInf):
    """
    Hook called at the start of solve.  This gives the class a chance to massage the input.
    For example, x and y might be dolfin (finite element) Functions; this method should return 
    the associated dolfin GenericVectors.

    The remaining arguments are passed directly from solve, and can be used for determining the
    final stopping criterion.

    Returns dolfin vectors corresponding to the initial value of x and the desired value of y=F(x).    
    """
    xv = PISMLocalVector(x)
    yv = PISMLocalVector(y)

    Jgoal = deltaLInf

    return (xv,yv,Jgoal)

  def finalize(self,x,y):
    """
    Hook called at the end of 'solve'.  Gives the chance to massage the return values.
    """

    tauc = x.core()
    u = y.core()
    return (tauc,u)

  def solve(self,x,y,deltaLInf):
    """
    Solve the ill posed problem F(x)=y where y is know to an L infinity error deltaLInf
    """
    return InvertIGN.solve(self,x,y,deltaLInf)

class PlotListener:  
  def __init__(self,grid):
    import tozero
    self.tz_scalar = tozero.ToProcZero(grid,dof=1)
    self.tz_vector = tozero.ToProcZero(grid,dof=2)
  
  def __call__(self,solver,count,x,Fx,y,d,r,*args):
    from matplotlib import pyplot as pp
    import siple
    d = self.tz_scalar.communicate(d.core())
    x = self.tz_scalar.communicate(x.core())
    r = self.tz_vector.communicate(r.core())
    y = self.tz_vector.communicate(y.core())
    Fx = self.tz_vector.communicate(Fx.core())

    if not d is None:
      r *= PISM.secpera
      y *= PISM.secpera

      self.iteration(solver,count,x,Fx,y,d,r,*args)

  def iteration(self,solver,count,x,Fx,y,d,r,*args):      
    import matplotlib.pyplot as pp
    pp.clf()
    pp.subplot(2,3,1)
    pp.imshow(y[0,:,:],origin='lower')
    pp.colorbar()
    pp.title('yu')
    pp.jet()

    pp.subplot(2,3,4)
    pp.imshow(y[1,:,:],origin='lower')
    pp.colorbar()
    pp.title('yv')
    pp.jet()

    
    pp.subplot(2,3,2)
    pp.imshow(r[0,:,:],origin='lower')
    pp.colorbar()
    pp.title('ru')
    pp.jet()

    pp.subplot(2,3,5)
    pp.imshow(r[1,:,:],origin='lower')
    pp.colorbar()
    pp.title('rv')
    pp.jet()

    d *= -1
    pp.subplot(2,3,3)      
    pp.imshow(d,origin='lower')
    pp.colorbar()
    pp.jet()
    pp.title('-d')
    
    pp.subplot(2,3,6)      
    pp.imshow(x,origin='lower')
    pp.colorbar()
    pp.title('zeta')
    pp.jet()

    pp.ion()
    pp.show()

def pauseListener(*args):
    import siple
    siple.reporting.pause()
