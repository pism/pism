import PISM
from PISM import IceModelVec2S, IceModelVec2V, PISMVars, secpera
import math

def init_shallow_grid(grid,Lx,Ly,Mx,My,p):
  grid.Lx = Lx;
  grid.Ly = Ly;
  grid.periodicity = p;
  grid.start_year = grid.year = 0.0;
  grid.Mx = Mx; grid.My=My; grid.Mz = 3;

  grid.compute_nprocs();
  grid.compute_ownership_ranges();
  grid.compute_vertical_levels()
  grid.compute_horizontal_spacing();
  grid.createDA();


class SSATestCase:
  def __init__(self,comm,rank,size,config):
    self.config = config
    self.grid = PISM.IceGrid(comm,rank,size,config)
    
    self.vars = PISMVars();

    self.surface = IceModelVec2S()
    self.thickness = IceModelVec2S()
    self.bed = IceModelVec2S()
    self.tauc = IceModelVec2S()

    self.enthalpy = PISM.IceModelVec3()

    self.vel_bc = IceModelVec2V()
    self.ice_mask =   PISM.IceModelVec2Mask()
    self.bc_mask =   PISM.IceModelVec2Mask()

    self.basal = None
    self.ice = None
    self.enthalpyconverter = None
    
  def buildSSACoefficients(self):
    WIDE_STENCIL=2;
    grid = self.grid

    # ice surface elevation
    self.surface.create(grid, "usurf", True, WIDE_STENCIL)
    self.surface.set_attrs("diagnostic", "ice upper surface elevation", "m", "surface_altitude");
    self.vars.add(self.surface);

    # land ice thickness
    self.thickness.create(grid, "thk", True, WIDE_STENCIL);
    self.thickness.set_attrs("model_state", "land ice thickness", "m", "land_ice_thickness");
    self.thickness.set_attr("valid_min", 0.0);
    self.vars.add(self.thickness);
 
    # bedrock surface elevation
    self.bed.create(grid, "topg", True, WIDE_STENCIL);
    self.bed.set_attrs("model_state", "bedrock surface elevation", "m", "bedrock_altitude");
    self.vars.add(self.bed)
 
    # yield stress for basal till (plastic or pseudo-plastic model)
    self.tauc.create(grid, "tauc", True, WIDE_STENCIL);
    self.tauc.set_attrs("diagnostic", "yield stress for basal till (plastic or pseudo-plastic model)", "Pa", ""); 
    self.vars.add(self.tauc)

    # enthalpy
    self.enthalpy.create(grid, "enthalpy", True, WIDE_STENCIL);
    self.enthalpy.set_attrs("model_state", "ice enthalpy (includes sensible heat, latent heat, pressure)", "J kg-1", "");
    self.vars.add(self.enthalpy);

    # dirichlet boundary condition
    self.vel_bc.create(grid, "_bc", True, WIDE_STENCIL);
    self.vel_bc.set_attrs("intent", "X-component of the SSA velocity boundary conditions", "m s-1", "", 0);
    self.vel_bc.set_attrs("intent", "Y-component of the SSA velocity boundary conditions", "m s-1", "", 1);
    self.vel_bc.set_glaciological_units("m year-1");
    huge_vel = 1e6/PISM.secpera;
    attrs = [ ("valid_min", -huge_vel), ("valid_max", huge_vel), ("_FillValue", 2*huge_vel) ]
    for a in attrs: 
      for component in range(2):
        self.vel_bc.set_attr(a[0],a[1],component)
    self.vel_bc.write_in_glaciological_units = True;
    self.vel_bc.set(2*huge_vel);

    # // grounded_dragging_floating integer mask
    self.ice_mask.create(grid, "mask", True, WIDE_STENCIL);
    self.ice_mask.set_attrs("model_state", "grounded_dragging_floating integer mask", "", "");
    mask_values=[PISM.MASK_ICE_FREE_BEDROCK, PISM.MASK_GROUNDED, PISM.MASK_FLOATING,
                  PISM.MASK_ICE_FREE_OCEAN, PISM.MASK_OCEAN_AT_TIME_0 ]
    self.ice_mask.set_attr("flag_values", mask_values);
    self.ice_mask.set_attr("flag_meanings","ice_free_bedrock dragging_sheet floating ice_free_ocean ocean_at_time_zero");
    self.ice_mask.output_data_type = PISM.NC_BYTE;
    self.vars.add(self.ice_mask);

    # // Dirichlet B.C. mask
    self.bc_mask.create(grid, "bc_mask", True, WIDE_STENCIL);
    self.bc_mask.set_attrs("model_state", "grounded_dragging_floating integer mask", "", "");
    mask_values = [0,1]
    self.bc_mask.set_attr("flag_values", mask_values)
    self.bc_mask.set_attr("flag_meanings","no_data dirichlet_bc_location");
    self.bc_mask.output_data_type = PISM.NC_BYTE;
    self.vars.add(self.bc_mask);


  def init(self,Mx,My):
      #// FIXME (DAM 2/17/11):  These are currently only looked at by the finite difference solver.
      self.config.scalar_from_option("ssa_eps",  "epsilon_ssa");
      self.config.scalar_from_option("ssa_maxi", "max_iterations_ssa");
      self.config.scalar_from_option("ssa_rtol", "ssa_relative_convergence");

      #// Subclass builds grid.
      self.initializeGrid(Mx,My);  

      #// We setup storage for the coefficients.
      self.buildSSACoefficients();

      #// Subclass builds ice flow law, basal resistance, etc.
      self.initializeSSAModel();

      # #// Allocate the actual SSA solver
      self.ssa = PISM.SSAFEM(self.grid, self.basal, self.ice, self.enthalpyconverter, self.config);
      #ssafactory(self.grid, self.basal, self.ice, self.enthalpyconverter, self.config);
      self.ssa.init(self.vars); 

      # // Allow the subclass to setup the coefficients.
      ierr = self.initializeSSACoefficients();

  #//! Solve the SSA
  def run(self):
    PISM.verbPrintf(2,self.grid.com,"* Solving the SSA stress balance ...\n");
    fast = False;
    self.ssa.update(fast);

  def report(self):
      grid = self.grid

      ssa_stdout = self.ssa.stdout_report()
      PISM.verbPrintf(3,grid.com,ssa_stdout)

      maxvecerr = 0.0; avvecerr = 0.0; 
      avuerr = 0.0; avverr = 0.0;
      maxuerr = 0.0; maxverr = 0.0;

      if(self.config.get_flag("do_pseudo_plastic_till")):
        PISM.verbPrintf(1,grid.com, "WARNING: numerical errors not valid for pseudo-plastic till\n")
      PISM.verbPrintf(1,grid.com, "NUMERICAL ERRORS in velocity relative to exact solution:\n")

      vel_ssa = self.ssa.get_advective_2d_velocity()
      
      vel_ssa.begin_access()

      exactvelmax = 0; gexactvelmax = 0;
      for i in range(self.grid.xs,self.grid.xs+self.grid.xm):
        for j in range(self.grid.ys,self.grid.ys+self.grid.ym):
          x=grid.x[i]; y=grid.y[j]
          (uexact,vexact) = self.exactSolution(i,j,x,y);
          exactnormsq=math.sqrt(uexact*uexact+vexact*vexact);
          exactvelmax = max(exactnormsq,exactvelmax);
          solution = vel_ssa[i,j]
          uerr = abs(solution.u-uexact)
          verr = abs(solution.v-vexact)
          avuerr += uerr;
          avverr += verr;
          maxuerr = max(maxuerr,uerr);
          maxverr = max(maxverr,verr)
          vecerr = math.sqrt(uerr * uerr + verr * verr);
          maxvecerr = max(maxvecerr,vecerr);
          avvecerr = avvecerr + vecerr;

      vel_ssa.end_access();
      
      gexactvelmax = PISM.globalMax(exactvelmax,grid.com);      
      gmaxuerr     = PISM.globalMax(maxuerr,grid.com);
      gmaxverr     = PISM.globalMax(maxverr,grid.com);
      gavuerr      = PISM.globalSum(avuerr,grid.com) / (grid.Mx*grid.My)
      gavverr      = PISM.globalSum(avverr,grid.com) / (grid.Mx*grid.My)
      gmaxvecerr   = PISM.globalMax(maxvecerr,grid.com)
      gavvecerr    = PISM.globalMax(avvecerr,grid.com) / (grid.Mx*grid.My)

      report_velocity_scale = PISM.secpera
      PISM.verbPrintf(1,grid.com, "velocity  :  maxvector   prcntavvec      maxu      maxv       avu       avv\n");
      #FIXME: variable arguments to verbPrintf are not working.  For now, do the string formatting on the python side.  Maybe
      #this is the best approach.
      PISM.verbPrintf(1,grid.com, "           %11.4f%13.5f%10.4f%10.4f%10.4f%10.4f\n" %
                            (gmaxvecerr*report_velocity_scale, (gavvecerr/gexactvelmax)*100.0,
                             gmaxuerr*report_velocity_scale, gmaxverr*report_velocity_scale, gavuerr*report_velocity_scale, 
                             gavverr*report_velocity_scale) );
      PISM.verbPrintf(1,grid.com, "NUM ERRORS DONE\n");


  def write(self,filename):
    grid = self.grid
    
    pio = PISM.PISMIO(grid)
    pio.open_for_writing(filename,False,True)
    pio.append_time(0.0)
    pio.close()
    for var in [self.surface,self.thickness,self.bc_mask,self.tauc,self.bed,self.enthalpy,self.vel_bc]:
      var.write(filename)
    
    vel_ssa = self.ssa.get_advective_2d_velocity()
    vel_ssa.write(filename)
    
    exact = IceModelVec2V();
    exact.create(grid, "_exact", False)
    exact.set_attrs("diagnostic", "X-component of the SSA exact solution", "m s-1", "", 0);
    exact.set_attrs("diagnostic", "U-component of the SSA exact solution", "m s-1", "", 1);
    exact.set_glaciological_units("m year-1");

    exact.begin_access()
    for i in xrange(grid.xs,grid.xs+grid.xm):
      for j in xrange(grid.ys,grid.ys+grid.ym):
        exact[i,j] = self.exactSolution(i,j,grid.x[i],grid.y[j])
    exact.end_access();
    exact.write(filename);

  def initializeGrid(self,Mx,My):
    raise NotImplementedError()
    
  def initializeSSAModel(self):
    raise NotImplementedError()
  
  def initializeSSACoefficients(self):
    raise NotImplementedError()

  def exactSolution(self,i,j,x,y):
    raise NotImplementedError()
    
