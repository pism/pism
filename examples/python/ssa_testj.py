import sys, petsc4py
petsc4py.init(sys.argv)
from petsc4py import PETSc

import PISM

import ssa

class testj(ssa.SSAExactTestCase):
  def initGrid(self,Mx,My):
    halfWidth = 300.0e3
    Lx = halfWidth; Ly = halfWidth
    ssa.init_shallow_grid(self.grid,Lx,Ly,Mx,My,PISM.XY_PERIODIC);

  def initPhysics(self):
    config = self.config
    self.basal = PISM.IceBasalResistancePlasticLaw(
           config.get("plastic_regularization") / PISM.secpera,
           config.get_flag("do_pseudo_plastic_till"),
           config.get("pseudo_plastic_q"),
           config.get("pseudo_plastic_uthreshold") / PISM.secpera);

    self.ice = PISM.CustomGlenIce(self.grid.com, "", config)
    self.enthalpyconverter = PISM.EnthalpyConverter(config)

  def initSSACoefficients(self):
    solver = self.solver
    solver.allocateBCs()

    solver.tauc.set(0.0) # irrelevant for test J
    solver.bed.set(0.0) # assures shelf is floating
    solver.bc_mask.set(0) # No dirichlet data.

    solver.enthalpy.set(528668.35); # arbitrary; corresponds to 263.15 Kelvin at depth=0.

    ocean_rho = self.config.get("sea_water_density");
    
    nu0 = 30.0 * 1.0e6 * PISM.secpera; # 9.45e14 Pa s 
    H0 = 500.0;                        # 500 m typical thickness

    solver.thickness.begin_access();
    solver.surface.begin_access();
    solver.bc_mask.begin_access();
    solver.vel_bc.begin_access();

    grid = self.grid
    for i in xrange(grid.xs,grid.xs+grid.xm):
      for j in xrange (grid.ys,grid.ys+grid.ym):
        x = grid.x[i]; y = grid.y[j]
        (H,junk,u,v) = PISM.exactJ(x,y);
        solver.thickness[i,j] = H;
        solver.surface[i,j] = (1.0 - self.ice.rho / ocean_rho) * H; #// FIXME task #7297
    
        # // special case at center point (Dirichlet BC)
        if (i == (grid.Mx)/2) and (j == (grid.My)/2):
          solver.bc_mask[i,j] = 1;
          solver.vel_bc[i,j] = [u,v]

    solver.surface.end_access();
    solver.thickness.end_access();
    solver.bc_mask.end_access();
    solver.vel_bc.end_access();

    # // communicate what we have set
    solver.surface.beginGhostComm(); solver.surface.endGhostComm();
    solver.thickness.beginGhostComm(); solver.thickness.endGhostComm();
    solver.bc_mask.beginGhostComm(); solver.bc_mask.endGhostComm();
    solver.vel_bc.beginGhostComm(); solver.vel_bc.endGhostComm();

    # Test J has a viscosity that is independent of velocity.  So we force a 
    # constant viscosity by settting the strength_extension
    # thickness larger than the given ice thickness. (max = 770m).

    solver.ssa.strength_extension.set_notional_strength(nu0 * H0);
    solver.ssa.strength_extension.set_min_thickness(800.);


  def exactSolution(self,i,j,x,y):
    (j1,j2,u,v) = PISM.exactJ(x,y);
    return [u,v]


# The main code for a run follows:

comm = PETSc.COMM_WORLD
rank = PETSc.Comm.getRank(comm)
size = PETSc.Comm.getSize(comm)

config = PISM.NCConfigVariable(); overrides = PISM.NCConfigVariable();

# Handle the input arguments.  I don't see a good way yet to call PetscOptionsBegin/End
# and integrate better into the -help mechanism.
optDB = PETSc.Options()
Mx = optDB.getInt("-Mx",default=61)
My = optDB.getInt("-My",default=61)
output_file = optDB.getString("-o",default="testj.nc")
verbosity = optDB.getString("-verbose",default=3)

PISM.init_config(comm, rank, config, overrides)
PISM.setVerbosityLevel(verbosity)
tc = testj(comm,rank,size,config,Mx,My)
tc.solve()
tc.report()
tc.write(output_file)
