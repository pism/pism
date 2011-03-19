# // Copyright (C) 2007--2011 Ed Bueler and Constantine Khroulev
# //
# // This file is part of PISM.
# //
# // PISM is free software; you can redistribute it and/or modify it under the
# // terms of the GNU General Public License as published by the Free Software
# // Foundation; either version 2 of the License, or (at your option) any later
# // version.
# //
# // PISM is distributed in the hope that it will be useful, but WITHOUT ANY
# // WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# // FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# // details.
# //
# // You should have received a copy of the GNU General Public License
# // along with PISM; if not, write to the Free Software
# // Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
# 
# #include <petsc.h>
# #include "SSAFD.hh"
# #include "SSAFEM.hh"
# #include "PISMIO.hh"
# #include "Timeseries.hh"
# 
# static char help[] =
#   "Driver for EISMINT-Ross diagnostic velocity computation in ice shelf.\n"
#   "Illustrates use of SSA stress balance plus I/O plus time series,\n"
#   "without the time-stepping mass continuity and conservation of energy\n"
#   "components of PISM.\n\n";
# 
# PetscErrorCode read_riggs_and_compare(IceGrid &grid, PISMVars &vars, IceModelVec2V &vel_ssa) {
#   PetscErrorCode  ierr;
#   PetscTruth      riggsSet;
#   char            riggsfile[PETSC_MAX_PATH_LEN];
# 
#   ierr = PetscOptionsGetString(PETSC_NULL, "-riggs", riggsfile,
#                                PETSC_MAX_PATH_LEN, &riggsSet); CHKERRQ(ierr);
#   if (riggsSet == PETSC_FALSE)
#     return 0;
# 
#   IceModelVec2S *longitude, *latitude;
#   IceModelVec2Mask *mask;
# 
#   longitude = dynamic_cast<IceModelVec2S*>(vars.get("longitude"));
#   if (longitude == NULL) SETERRQ(1, "longitude is not available");
# 
#   latitude = dynamic_cast<IceModelVec2S*>(vars.get("latitude"));
#   if (latitude == NULL) SETERRQ(1, "latitude is not available");
# 
#   mask = dynamic_cast<IceModelVec2Mask*>(vars.get("mask"));
#   if (mask == NULL) SETERRQ(1, "mask is not available");
# 
#   ierr = verbPrintf(2,grid.com,"comparing to RIGGS data in %s ...\n",
#                     riggsfile); CHKERRQ(ierr);
# 
#   Timeseries latdata(&grid, "riggslat", "count"),
#     londata(&grid, "riggslon", "count"),
#     magdata(&grid, "riggsmag", "count"),
#     udata(&grid, "riggsu", "count"),
#     vdata(&grid, "riggsv", "count");
#   PetscInt    len;
#   PetscScalar **clat, **clon;
# 
#   ierr = longitude->get_array(clon); CHKERRQ(ierr);
#   ierr =  latitude->get_array(clat); CHKERRQ(ierr);
#   ierr =    vel_ssa.begin_access(); CHKERRQ(ierr);
#   ierr =      mask->begin_access();  CHKERRQ(ierr);
# 
#   ierr = magdata.set_units("m year-1", ""); CHKERRQ(ierr);
#   ierr =   udata.set_units("m year-1", ""); CHKERRQ(ierr);
#   ierr =   vdata.set_units("m year-1", ""); CHKERRQ(ierr);
# 
#   ierr = latdata.read(riggsfile); CHKERRQ(ierr);
#   ierr = londata.read(riggsfile); CHKERRQ(ierr);
#   ierr = magdata.read(riggsfile); CHKERRQ(ierr);
#   ierr =   udata.read(riggsfile); CHKERRQ(ierr);
#   ierr =   vdata.read(riggsfile); CHKERRQ(ierr);
#       
#   // same length for all vars here
#   len = latdata.length();
#   PetscScalar  goodptcount = 0.0, ChiSqr = 0.0;
#   for (PetscInt k = 0; k<len; k++) {
#     PetscScalar lat, lon, mag, u, v;
#     lat = latdata[k];
#     lon = londata[k];
#     mag = magdata[k];
#     u   = udata[k];
#     v   = vdata[k];
#     ierr = verbPrintf(4,grid.com,
#                       " RIGGS[%3d]: lat = %7.3f, lon = %7.3f, mag = %7.2f, u = %7.2f, v = %7.2f\n",
#                       k,lat,lon,mag,u,v); CHKERRQ(ierr); 
#     const PetscScalar origdlat = (-5.42445 - (-12.3325)) / 110.0;
#     const PetscScalar lowlat = -12.3325 - origdlat * 46.0;
#     const PetscScalar dlat = (-5.42445 - lowlat) / (float) (grid.My - 1);        
#     const PetscScalar lowlon = -5.26168;
#     const PetscScalar dlon = (3.72207 - lowlon) / (float) (grid.Mx - 1);
#     const int         cj = (int) floor((lat - lowlat) / dlat);
#     const int         ci = (int) floor((lon - lowlon) / dlon);
#     if ((ci >= grid.xs) && (ci < grid.xs+grid.xm) && (cj >= grid.ys) && (cj < grid.ys+grid.ym)) {
#       const PetscScalar cu = secpera * vel_ssa(ci,cj).u;
#       const PetscScalar cv = secpera * vel_ssa(ci,cj).v;
#       const PetscScalar cmag = sqrt(PetscSqr(cu)+PetscSqr(cv));
#       ierr = verbPrintf(4,PETSC_COMM_SELF,
#                         " PISM%d[%3d]: lat = %7.3f, lon = %7.3f, mag = %7.2f, u = %7.2f, v = %7.2f\n",
#                         grid.rank,k,clat[ci][cj],clon[ci][cj],cmag,cu,cv); CHKERRQ(ierr); 
#       if (mask->value(ci,cj) == MASK_FLOATING) {
#         goodptcount += 1.0;
#         ChiSqr += PetscSqr(u-cu)+PetscSqr(v-cv);
#       }
#     }
#   }
#   ChiSqr = ChiSqr / PetscSqr(30.0); // see page 48 of MacAyeal et al
#   PetscScalar g_goodptcount, g_ChiSqr;
#   ierr = PetscGlobalSum(&goodptcount, &g_goodptcount, grid.com); CHKERRQ(ierr);
#   ierr = PetscGlobalSum(&ChiSqr, &g_ChiSqr, grid.com); CHKERRQ(ierr);
#   ierr = verbPrintf(4,grid.com,"number of RIGGS data points = %d\n"
#                     "number of RIGGS points in computed ice shelf region = %8.2f\n",
#                     len, g_goodptcount); CHKERRQ(ierr);
#   ierr = verbPrintf(2,grid.com,"Chi^2 statistic for computed results compared to RIGGS is %10.3f\n",
#                     g_ChiSqr * (156.0 / g_goodptcount)); CHKERRQ(ierr);
# 
#   ierr = longitude->end_access(); CHKERRQ(ierr);
#   ierr =  latitude->end_access(); CHKERRQ(ierr);
#   ierr =      mask->end_access(); CHKERRQ(ierr);
#   ierr =    vel_ssa.end_access(); CHKERRQ(ierr);
# 
#   return 0;
# }
# 

  



import sys, petsc4py
petsc4py.init(sys.argv)
from petsc4py import PETSc

import PISM, ssa, math


DEFAULT_MIN_THICKNESS = 5.0 # meters
DEFAULT_CONSTANT_HARDNESS_FOR_SSA = 1.9e8  # Pa s^{1/3}; see p. 49 of MacAyeal et al 1996
DEFAULT_TYPICAL_STRAIN_RATE = (100.0 / PISM.secpera) / (100.0 * 1.0e3)  # typical strain rate is 100 m/yr per 
DEFAULT_nuH = DEFAULT_MIN_THICKNESS * DEFAULT_CONSTANT_HARDNESS_FOR_SSA / \
(2.0 * pow(DEFAULT_TYPICAL_STRAIN_RATE,2./3.)); # Pa s m
# COMPARE: 30.0 * 1e6 * secpera = 9.45e14 is Ritz et al (2001) value of
#          30 MPa yr for \bar\nu

class pross(ssa.SSATestCase):
  def __init__(self,comm,rank,size,config,Mx,My,boot_file,riggs_file):
    # Set the file names so that grid construction knows where to find
    # the boot_file
    self.boot_file = boot_file
    self.riggs_file = riggs_file

    ssa.SSATestCase.__init__(self,comm,rank,size,config,Mx,My)

    # the superclass built us a grid.
    grid = self.grid
    
    self.obsAzimuth = PISM.IceModelVec2S()
    self.obsAzimuth.create(grid, "azi_obs", True);
    self.obsAzimuth.set_attrs("", "observed ice velocity azimuth", "degrees_east", "");

    self.obsMagnitude = PISM.IceModelVec2S()
    self.obsMagnitude.create(grid, "mag_obs", True)
    self.obsMagnitude.set_attrs("", "observed ice velocity magnitude", "m s-1", "");
    self.obsMagnitude.set_glaciological_units("m year-1");
    self.obsMagnitude.write_in_glaciological_units = True;

    self.obsAccurate = PISM.IceModelVec2S()
    self.obsAccurate.create(grid, "accur", True)
    self.obsAccurate.set_attrs("", "flag for accurate observed velocity","", "");

    self.longitude = ssa.standardLongitudeVec(grid)
    self.latitude  = ssa.standardLatitudeVec(grid)
    vars = set([self.longitude, self.latitude,self.obsAzimuth,self.obsAccurate,self.obsMagnitude])
    for v in vars:
      v.regrid(boot_file,True)

  def initGrid(self,Mx,My):
    ssa.grid2d_init_from_file(self.grid,self.boot_file,Mx,My)
    self.grid.printInfo(1)

  def initPhysics(self):
    secpera = PISM.secpera
    self.ice = PISM.CustomGlenIce(self.grid.com,"",config)
    self.ice.setHardness(DEFAULT_CONSTANT_HARDNESS_FOR_SSA);
    self.basal = PISM.IceBasalResistancePlasticLaw(  config.get("plastic_regularization") / secpera, 
                                           config.get_flag("do_pseudo_plastic_till"),
                                           config.get("pseudo_plastic_q"),
                                           config.get("pseudo_plastic_uthreshold") / secpera)
    self.basal.printInfo(1,self.grid.com)
    self.enthalpyconverter = PISM.EnthalpyConverter(config)

  def initSSACoefficients(self):
    solver = self.solver
    solver.allocateBCs(velname='bar',maskname='bcflag')

    solver.readCoeffsFromFile(self.boot_file,omit=[solver.tauc,solver.enthalpy,solver.surface])
    solver.enthalpy.set(528668.35)  # Corresponds to 263.15 Kelvin at depth=0.
                                    # The CustomGlenIce flow law does not use it.
    # set the basal yield stress (does not matter; everything is floating)
    solver.tauc.set(0.0)

    ice_rho = config.get("ice_density")
    ocean_rho = config.get("sea_water_density")
    solver.thickness.copy_to(solver.surface)
    solver.surface.scale(1.0 - ice_rho / ocean_rho)

    solver.ssa.strength_extension.set_min_thickness(DEFAULT_MIN_THICKNESS);
    solver.ssa.strength_extension.set_notional_strength(DEFAULT_nuH);

  def report(self):
    uerr = 0.0; verr=0.0; relvecerr=0.0; accN=0.0 
    accArea=0.0; maxcComputed=0.0; vecErrAcc = 0.0;
    
    grid = self.grid
    area = grid.dx*grid.dy
  
    mask = self.solver.ice_mask; mask.begin_access();
    H = self.solver.thickness; H.begin_access();
    azi=self.obsAzimuth; azi.begin_access();
    mag = self.obsMagnitude; mag.begin_access();
    acc = self.obsAccurate; acc.begin_access();
    vel_ssa = self.solver.solution(); vel_ssa.begin_access()
    
    for i in xrange(grid.xs,grid.xs+grid.xm):
      for j in xrange(grid.ys,grid.ys+grid.ym):
        if mask.is_floating(i,j) and H[i,j] > 1.0:
          ccomputed = vel_ssa[i,j].magnitude()
          maxcComputed = max(maxcComputed,ccomputed)
          if( abs(acc[i,j]-1.0) < 0.1):
            accN += 1.0
            accArea += area
            uobs = mag[i,j] * math.sin((math.pi/180.0) * azi[i,j]);
            vobs = mag[i,j] * math.cos((math.pi/180.0) * azi[i,j]);
            Dv = abs(vobs-vel_ssa[i,j].v)
            Du = abs(uobs-vel_ssa[i,j].u)
            verr += Dv; uerr += Du
            relvecerr += (Dv*Dv+Du*Du) / (vobs*vobs+uobs*uobs)
            vecErrAcc += (Dv*Dv+Du*Du) * area
    mag.end_access()
    azi.end_access()
    acc.end_access()
    vel_ssa.end_access()
    mask.end_access()
    H.end_access()

    gmaxcComputed = PISM.globalMax(maxcComputed,grid.com)
    gaccN         = PISM.globalSum(accN, grid.com)
    gaccArea      = PISM.globalSum(accArea, grid.com)
    gverr         = PISM.globalSum(verr, grid.com)
    guerr         = PISM.globalSum(uerr, grid.com)
    grelvecerr    = PISM.globalSum(relvecerr,grid.com)
    gvecErrAcc    = PISM.globalSum(vecErrAcc, grid.com)

    secpera = PISM.secpera
    r=ssa.Reporter(self.grid.com,verbosity=2)
    r.println("maximum computed speed in ice shelf is %10.3f (m/a)", gmaxcComputed * secpera);
    r.println("ERRORS relative to observations of Ross Ice Shelf:");
    r.println("  [number of grid points in 'accurate observed area' = %d]", int(gaccN))
    r.println("  [area of 'accurate observed area' = %9.4f (km^2)]", gaccArea / 1e6)
    r.println("  following are average errors computed over 'accurate observed area':");
    r.println("  average error in x-comp of vel       = %9.3f (m/a)", (gverr * secpera) / gaccN);
    r.println("  average error in y-comp of vel       = %9.3f (m/a)",(guerr * secpera) / gaccN);
    r.println("  average relative error in vector vel = %9.5f", grelvecerr / gaccN)
    gvecErrAcc = secpera * math.sqrt(gvecErrAcc) / math.sqrt(gaccArea);
    r.println("  rms average error in vector vel      = %9.3f (m/a)\n", gvecErrAcc);
    

  def write(self,filename):
    ssa.SSATestCase.write(self,filename)

    vars = set([self.longitude, self.latitude,self.obsAzimuth,self.obsAccurate,self.obsMagnitude])
    for v in vars:
      v.write(filename)

    cbar = PISM.IceModelVec2S();
    cbar.create(self.grid,"cbar",False)
    cbar.set_attrs( "diagnostic", 
                    "magnitude of vertically-integrated horizontal velocity of ice",
                    "m s-1", "");
    cbar.set_glaciological_units("m year-1")
    cbar.write_in_glaciological_units = True;
    self.solver.solution().magnitude(cbar)
    cbar.write(filename)


# The main code for a run follows:

com = PETSc.COMM_WORLD
rank = PETSc.Comm.getRank(com)
size = PETSc.Comm.getSize(com)


PISM.verbosityLevelFromOptions()
PISM.verbPrintf(2,com,"PROSS %s (EISMINT-Ross diagnostic velocity computation mode)\n", PISM.PISM_Revision)
PISM.stop_on_version_option()
usage = \
"""  pross -boot_file IN.nc -Mx number -My number [-o file.nc] [-riggs file.nc]
where:
-boot_file  IN.nc is input file in NetCDF format: contains PISM-written model state
  -Mx         number of grid points in the x direction
  -My         number of grid points in the y direction
  -riggs      read RIGGS data from a file
notes:
  * -boot_file is required
"""

PISM.show_usage_check_req_opts(com,"pross",["-boot_file"],usage)

config = PISM.NCConfigVariable(); overrides = PISM.NCConfigVariable();
PISM.init_config(com, rank, config, overrides)
config.set_flag("use_ssa_velocity", True)
config.set_flag("use_ssa_when_grounded", False)
config.set_flag("use_constant_nuh_for_ssa", False)
config.set("epsilon_ssa", 0.0)

#config.scalar_from_option("ssa_rtol", "ssa_relative_convergence")

# YUCK
optDB = PETSc.Options()
Mx = optDB.getInt("-Mx",default=0)
if Mx == 0: Mx = None
My = optDB.getInt("-My",default=0)
if My == 0: My = None
boot_file = optDB.getString("-boot_file")
riggs_file = None

output_file = optDB.getString("-o",default="ross_computed.nc")

test_case = pross(com,rank,size,config,Mx,My,boot_file,riggs_file)
test_case.solve()
test_case.report()
test_case.write(output_file)
