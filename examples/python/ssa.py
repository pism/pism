# Copyright (C) 2011 David Maxwell
# 
# This file is part of PISM.
# 
# PISM is free software; you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation; either version 2 of the License, or (at your option) any later
# version.
# 
# PISM is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.
# 
# You should have received a copy of the GNU General Public License
# along with PISM; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

import PISM
from PISM import IceModelVec2S, IceModelVec2V, PISMVars, secpera
from petsc4py import PETSc
import math

WIDE_STENCIL=2

def standardIceSurfaceVec(grid,name='usurf'):
  surface = IceModelVec2S();
  surface.create(grid, name, True, WIDE_STENCIL)
  surface.set_attrs("diagnostic", "ice upper surface elevation", "m", "surface_altitude");
  return surface;

def standardIceThicknessVec(grid,name='thk'):
  thickness = IceModelVec2S();
  thickness.create(grid, name, True, WIDE_STENCIL);
  thickness.set_attrs("model_state", "land ice thickness", "m", "land_ice_thickness");
  thickness.set_attr("valid_min", 0.0);
  return thickness

def standardBedrockElevationVec(grid,name='topg'):
  bed = IceModelVec2S()
  bed.create(grid, name, True, WIDE_STENCIL);
  bed.set_attrs("model_state", "bedrock surface elevation", "m", "bedrock_altitude");
  return bed

def standardYieldStressVec(grid,name='tauc'):
  # yield stress for basal till (plastic or pseudo-plastic model)
  tauc = IceModelVec2S()
  tauc.create(grid, name, True, WIDE_STENCIL);
  tauc.set_attrs("diagnostic", "yield stress for basal till (plastic or pseudo-plastic model)", "Pa", ""); 
  return tauc;

def standardEnthalpyVec(grid,name='enthalpy'):
  enthalpy = PISM.IceModelVec3()
  enthalpy.create(grid, name, True, WIDE_STENCIL);
  enthalpy.set_attrs("model_state", "ice enthalpy (includes sensible heat, latent heat, pressure)", "J kg-1", "");
  return enthalpy;

def standard2dVelocityVec(grid,name="",desc="",intent=""):
  if name is None: name="" # FIXME
  vel = IceModelVec2V();
  vel.create(grid,name,True,WIDE_STENCIL)
  vel.set_attrs(intent, "%s%s" %("X-component of the ",desc), "m s-1", "", 0);
  vel.set_attrs(intent, "%s%s" %("Y-component of the ",desc), "m s-1", "", 1);
  vel.set_glaciological_units("m year-1");
  vel.write_in_glaciological_units = True
  huge_vel = 1e6/PISM.secpera;
  attrs = [ ("valid_min", -huge_vel), ("valid_max", huge_vel), ("_FillValue", 2*huge_vel) ]
  for a in attrs: 
    for component in range(2):
      vel.set_attr(a[0],a[1],component)
  vel.set(2*huge_vel)
  return vel

def standardIceMask(grid,name='mask'):
  ice_mask = PISM.IceModelVec2Mask()
  ice_mask.create(grid, name, True, WIDE_STENCIL);
  ice_mask.set_attrs("model_state", "grounded_dragging_floating integer mask", "", "");
  mask_values=[PISM.MASK_ICE_FREE_BEDROCK, PISM.MASK_GROUNDED, PISM.MASK_FLOATING,
               PISM.MASK_ICE_FREE_OCEAN, PISM.MASK_OCEAN_AT_TIME_0 ]
  ice_mask.set_attr("flag_values", mask_values);
  ice_mask.set_attr("flag_meanings","ice_free_bedrock dragging_sheet floating ice_free_ocean ocean_at_time_zero");
  ice_mask.output_data_type = PISM.NC_BYTE;
  return ice_mask

def standardBCMask(grid,name='bc_mask'):
  bc_mask = PISM.IceModelVec2Mask()
  bc_mask.create(grid, name, True, WIDE_STENCIL);
  bc_mask.set_attrs("model_state", "grounded_dragging_floating integer mask", "", "");
  mask_values=[0,1]
  bc_mask.set_attr("flag_values", mask_values);
  bc_mask.set_attr("flag_meanings","no_data dirichlet_bc_location");
  bc_mask.output_data_type = PISM.NC_BYTE;
  return bc_mask


def standardLongitudeVec(grid,name="lon"):
  longitude = PISM.IceModelVec2S()
  longitude.create(grid, name, True)
  longitude.set_attrs("mapping", "longitude", "degree_east", "longitude")
  longitude.time_independent = True
  longitude.set_attr("coordinates", "")
  longitude.set_attr("grid_mapping", "")
  return longitude

def standardLatitudeVec(grid,name="lat"):
  latitude = PISM.IceModelVec2S()
  latitude.create(grid, name, True)
  latitude.set_attrs("mapping", "latitude", "degree_east", "latitude")
  latitude.time_independent = True
  latitude.set_attr("coordinates", "")
  latitude.set_attr("grid_mapping", "")
  return latitude


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

def grid2d_init_from_file(grid,bootfile,Mx=None,My=None):
  pio = PISM.PISMIO(grid)
  pio.open_for_reading(bootfile)
  grid_info = PISM.grid_info(); pio.get_grid_info_2d( grid_info )
  pio.close()

  grid.Mx = grid_info.x_len;
  grid.My = grid_info.y_len;
  grid.Mz = 2;
  grid.x0 = grid_info.x0;
  grid.y0 = grid_info.y0;
  grid.Lx = grid_info.Lx;
  grid.Ly = grid_info.Ly;

  if not Mx is None:
    grid.Mx = Mx
  if not My is None:
    grid.My = My

  grid.compute_nprocs();
  grid.compute_ownership_ranges();
  grid.compute_horizontal_spacing();
  grid.compute_vertical_levels();
  grid.createDA()


SSAAlgorithms = {"fem":PISM.SSAFEM, "fd":PISM.SSAFD}

class SSASolver:
  def __init__(self,grid,ice,basal,enthalpyconverter,config=None):
    self.grid = grid
    if config is None: config = grid.config
    self.config = config
  
    self.surface    = standardIceSurfaceVec( grid )
    self.thickness  = standardIceThicknessVec( grid )
    self.bed        = standardBedrockElevationVec( grid )
    self.tauc       = standardYieldStressVec( grid )
    self.enthalpy   = standardEnthalpyVec( grid )
    self.ice_mask   = standardIceMask( grid )
    
    self.vel_bc     = None
    self.bc_mask    = None
  
    self.ice   = ice
    self.basal = basal
    self.enthalpyconverter = enthalpyconverter

    self.setFromOptions()
    self.ssa = SSAAlgorithms[config.get_string("ssa_method")](self.grid,self.basal,self.ice,self.enthalpyconverter,self.config)

    self.vel_ssa = None

    self.ssa_init = False

  def setFromOptions(self):
    #// FIXME (DAM 2/17/11):  These are currently only looked at by the finite difference solver.
    self.config.scalar_from_option("ssa_eps",  "epsilon_ssa");
    self.config.scalar_from_option("ssa_maxi", "max_iterations_ssa");
    self.config.scalar_from_option("ssa_rtol", "ssa_relative_convergence");

    optDB = PETSc.Options()
    if optDB.hasName("-ssa_method"):
      # I'd really like to call PISMOptionsList here.
      ssa_method = optDB.getString("-ssa_method")
      if SSAAlgorithms.has_key(ssa_method):
        self.config.set_string("ssa_method",ssa_method)
      else:
        raise ValueError("-ssa_method should be one of %s, found %s" %(SSAAlgorithms.keys(),ssa_method))

  def allocateBCs(self,velname='_bc',maskname='bc_mask'):
    self.vel_bc     = standard2dVelocityVec( self.grid, name=velname, desc='SSA velocity boundary condition',intent='intent' )
    self.bc_mask    = standardBCMask( self.grid,name=maskname )

  def readCoeffsFromFile(self,filename,omit=None):
    vars = set([self.surface,self.thickness,self.bed,self.tauc,self.enthalpy,self.ice_mask])
    if not self.vel_bc is None:
      vars.add(self.vel_bc); vars.add(self.bc_mask)
    if not omit is None:
      vars.difference_update(omit)
    
    for v in vars:
      v.regrid(filename,True)

  def solve(self):
    if not self.ssa_init:
      pismVars = PISM.PISMVars()
      for var in [self.surface,self.thickness,self.bed,self.tauc,self.enthalpy,self.ice_mask]:
        pismVars.add(var)
        
      self.ssa.init(pismVars)
      if not self.vel_bc is None:
        self.ssa.set_boundary_conditions(self.bc_mask,self.vel_bc)
      self.ssa_init = True
    
    fast = False;
    self.ssa.update(fast);
    self.vel_ssa = self.ssa.get_advective_2d_velocity()

  def solution(self):
    return self.vel_ssa

  def write(self,filename):
    vars = [self.surface,self.thickness,self.tauc,self.bed,self.enthalpy]
    if not self.bc_mask is None:
      vars.append(self.bc_mask)
      vars.append(self.vel_bc)

    for var in vars:
      var.write(filename)

    self.vel_ssa.write(filename)


class SSATestCase:
  def __init__(self,comm,rank,size,config,Mx,My):
    self.grid = PISM.IceGrid(comm,rank,size,config)
    self.config = config

    self.initGrid(Mx,My)

    self.ice = None; self.basal = None; self.enthalpyconverter = None
    self.initPhysics()
    # FIXME: Check that the subclass did its job.

    self.solver = SSASolver(self.grid,self.ice,self.basal,self.enthalpyconverter)

    self.initSSACoefficients()

  #//! Solve the SSA
  def solve(self):
    PISM.verbPrintf(2,self.grid.com,"* Solving the SSA stress balance ...\n");
    self.solver.solve()

  def write(self,filename):
    grid = self.grid
    
    pio = PISM.PISMIO(grid)
    pio.open_for_writing(filename,False,True)
    pio.append_time(0.0)
    pio.close()
    
    self.solver.write(filename)

  def initGrid(self,Mx,My):
    raise NotImplementedError()
    
  def initPhysics(self):
    raise NotImplementedError()
  
  def initSSACoefficients(self):
    raise NotImplementedError()

  def report(self):
    raise NotImplementedError()


class SSAExactTestCase(SSATestCase):
  def report(self):
      grid = self.grid

      ssa_stdout = self.solver.ssa.stdout_report()
      PISM.verbPrintf(3,grid.com,ssa_stdout)

      maxvecerr = 0.0; avvecerr = 0.0; 
      avuerr = 0.0; avverr = 0.0;
      maxuerr = 0.0; maxverr = 0.0;

      if(self.config.get_flag("do_pseudo_plastic_till")):
        PISM.verbPrintf(1,grid.com, "WARNING: numerical errors not valid for pseudo-plastic till\n")
      PISM.verbPrintf(1,grid.com, "NUMERICAL ERRORS in velocity relative to exact solution:\n")

      vel_ssa = self.solver.ssa.get_advective_2d_velocity()
      
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
    SSATestCase.write(self,filename)
    
    grid=self.grid
    exact = standard2dVelocityVec(grid,name="_exact",desc="SSA exact solution",intent="diagnostic")
    exact.begin_access()
    for i in xrange(grid.xs,grid.xs+grid.xm):
      for j in xrange(grid.ys,grid.ys+grid.ym):
        exact[i,j] = self.exactSolution(i,j,grid.x[i],grid.y[j])
    exact.end_access();
    exact.write(filename);
  

class Reporter:
  
  def __init__(self,com,verbosity=1):
    self.com = com
    self.verbosity=verbosity
  
  def println(self,msg,*args):
    PISM.verbPrintf(self.verbosity,self.com,"%s\n" % (msg % args))
