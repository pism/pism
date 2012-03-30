#! /usr/bin/env python
#
# Copyright (C) 2011, 2012 Ed Bueler and Constantine Khroulev and David Maxwell
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


import sys, petsc4py
petsc4py.init(sys.argv)
from petsc4py import PETSc

import PISM, math


DEFAULT_MIN_THICKNESS = 5.0 # meters
DEFAULT_CONSTANT_HARDNESS_FOR_SSA = 1.9e8  # Pa s^{1/3}; see p. 49 of MacAyeal et al 1996
DEFAULT_TYPICAL_STRAIN_RATE = (100.0 / PISM.secpera) / (100.0 * 1.0e3)  # typical strain rate is 100 m/yr per 
DEFAULT_nuH = DEFAULT_MIN_THICKNESS * DEFAULT_CONSTANT_HARDNESS_FOR_SSA / \
(2.0 * pow(DEFAULT_TYPICAL_STRAIN_RATE,2./3.)); # Pa s m
# COMPARE: 30.0 * 1e6 * secpera = 9.45e14 is Ritz et al (2001) value of
#          30 MPa yr for \bar\nu

class pross(PISM.ssa.SSARun):
  def __init__(self,Mx,My,boot_file,riggs_file):
    PISM.ssa.SSARun.__init__(self)
    self.Mx = Mx; self.My = My;
    
    # Set the file names so that grid construction knows where to find
    # the boot_file
    self.boot_file = boot_file
    self.riggs_file = riggs_file

  def _initGrid(self):
    self.grid = PISM.Context().newgrid()
    grid = self.grid
    
    pio = PISM.PIO(grid.com, grid.rank, "netcdf3")
    pio.open(self.boot_file, PISM.NC_NOWRITE)
    grid_info = PISM.grid_info();
    pio.inq_grid_info( "land_ice_thickness", grid_info)
    pio.close()

    grid.Mx = grid_info.x_len;
    grid.My = grid_info.y_len;
    grid.Mz = 2;
    grid.x0 = (grid_info.x_max+grid_info.x_min)/2.;
    grid.y0 = (grid_info.y_max+grid_info.y_min)/2.;
    grid.Lx = (grid_info.x_max-grid_info.x_min)/2.;
    grid.Ly = (grid_info.x_max-grid_info.x_min)/2.;

    if not self.Mx is None:
      grid.Mx = self.Mx
    if not self.My is None:
      grid.My = self.My

    grid.compute_nprocs();
    grid.compute_ownership_ranges();
    grid.compute_horizontal_spacing();
    grid.compute_vertical_levels();
    grid.createDA()


    self.grid.printInfo(1)

  def _initPhysics(self):
    secpera = PISM.secpera
    enthalpyconverter = PISM.EnthalpyConverter(config)

    config.set_string("ssa_flow_law", "isothermal_glen")
    config.set("ice_softness", pow(DEFAULT_CONSTANT_HARDNESS_FOR_SSA, -config.get("Glen_exponent")))

    basal = PISM.IceBasalResistancePlasticLaw(  config.get("plastic_regularization") / secpera, 
                                           config.get_flag("do_pseudo_plastic_till"),
                                           config.get("pseudo_plastic_q"),
                                           config.get("pseudo_plastic_uthreshold") / secpera)
    basal.printInfo(1,self.grid.com)
    self.modeldata.setPhysics(basal,enthalpyconverter)

  def _initSSACoefficients(self):
    self._allocStdSSACoefficients()
    self._allocateBCs(velname='bar',maskname='bcflag')
    vecs = self.modeldata.vecs
    grid = self.modeldata.grid

    obsAzimuth = PISM.IceModelVec2S()
    obsAzimuth.create(grid, "azi_obs", True);
    obsAzimuth.set_attrs("", "observed ice velocity azimuth", "degrees_east", "");
    vecs.add(obsAzimuth,"obsAzimuth")

    obsMagnitude = PISM.IceModelVec2S()
    obsMagnitude.create(grid, "mag_obs", True)
    obsMagnitude.set_attrs("", "observed ice velocity magnitude", "m s-1", "");
    obsMagnitude.set_glaciological_units("m year-1");
    obsMagnitude.write_in_glaciological_units = True;
    vecs.add(obsMagnitude,"obsMagnitude")

    obsAccurate = PISM.IceModelVec2S()
    obsAccurate.create(grid, "accur", True)
    obsAccurate.set_attrs("", "flag for accurate observed velocity","", "");
    vecs.add(obsAccurate,"obsAccurate")

    longitude = PISM.util.standardLongitudeVec(grid)
    latitude  = PISM.util.standardLatitudeVec(grid)
    vecs.add(longitude,"longitude")
    vecs.add(latitude,"latitude")
    
    for v in [longitude, latitude,obsAzimuth,obsAccurate,obsMagnitude]:
      v.regrid(self.boot_file,True)

    thickness = vecs.thickness; bed = vecs.bed; enthalpy = vecs.enthalpy
    mask = vecs.ice_mask; surface = vecs.surface; tauc = vecs.tauc

    for v in [bed,mask,thickness,vecs.vel_bc,vecs.bc_mask]:
      v.regrid(self.boot_file,True)

    enthalpy.set(528668.35)  # Corresponds to 263.15 Kelvin at depth=0.
                             # The CustomGlenIce flow law does not use it.
    # set the basal yield stress (does not matter; everything is floating)
    tauc.set(0.0)

    surface.copy_from(thickness)
    ice_rho = config.get("ice_density")
    ocean_rho = config.get("sea_water_density")
    surface.scale(1.0 - ice_rho / ocean_rho)

  def _initSSA(self):
    self.ssa.strength_extension.set_min_thickness(DEFAULT_MIN_THICKNESS);
    self.ssa.strength_extension.set_notional_strength(DEFAULT_nuH);

  def report(self):
    uerr = 0.0; verr=0.0; relvecerr=0.0; accN=0.0 
    accArea=0.0; maxcComputed=0.0; vecErrAcc = 0.0;
    
    grid = self.grid
    area = grid.dx*grid.dy
  
    vecs = self.modeldata.vecs
    
    mask = vecs.ice_mask;
    H = vecs.thickness;
    azi= vecs.obsAzimuth;
    mag = vecs.obsMagnitude;
    acc = vecs.obsAccurate;
    vel_ssa = vecs.vel_ssa

    m = PISM.MaskQuery(mask)
    
    with PISM.util.Access([mask,H,azi,mag,acc,vel_ssa]):
      for (i,j) in grid.points():
        if m.ocean(i,j) and H[i,j] > 1.0:
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

    gmaxcComputed = PISM.globalMax(maxcComputed,grid.com)
    gaccN         = PISM.globalSum(accN, grid.com)
    gaccArea      = PISM.globalSum(accArea, grid.com)
    gverr         = PISM.globalSum(verr, grid.com)
    guerr         = PISM.globalSum(uerr, grid.com)
    grelvecerr    = PISM.globalSum(relvecerr,grid.com)
    gvecErrAcc    = PISM.globalSum(vecErrAcc, grid.com)

    secpera = PISM.secpera
    r=PISM.VerbPrintf(self.grid.com,verbosity=2)
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
    
    if not self.riggs_file is None:
      self.report_riggs()
    
  def report_riggs(self):
    grid = self.grid
    r=PISM.VerbPrintf(grid.com,verbosity=2)
    r.println("comparing to RIGGS data in %s ...\n",self.riggs_file); 

    riggsvars = [ "riggslat", "riggslon", "riggsmag", "riggsu", "riggsv"]
    (latdata,longdata,magdata,udata,vdata) = \
         (PISM.Timeseries(grid,varname,"count") for varname in riggsvars)
    riggsvars = [latdata, longdata, magdata, udata, vdata]

    for v in [magdata, udata, vdata]:
      v.set_units("m year-1", "")

    for v in riggsvars:
      v.read(self.riggs_file, False)

    length = latdata.length();

    vecs = self.modeldata.vecs
    vel_ssa = vecs.vel_ssa;
    clat = vecs.latitude; clon = vecs.longitude; mask = vecs.ice_mask;

    secpera=PISM.secpera
    with PISM.util.Access([clat,clon,mask,vel_ssa]):
      goodptcount = 0.0; ChiSqr = 0.0;    
      for k in xrange(length):
        (lat,lon,mag,u,v) = [v[k] for v in riggsvars]
        r.printlnv(4," RIGGS[%3d]: lat = %7.3f, lon = %7.3f, mag = %7.2f, u = %7.2f, v = %7.2f",
                          k,lat,lon,mag,u,v)
        origdlat = (-5.42445 - (-12.3325)) / 110.0;
        lowlat = -12.3325 - origdlat * 46.0;
        dlat = (-5.42445 - lowlat) / (float) (grid.My - 1);        
        lowlon = -5.26168;
        dlon = (3.72207 - lowlon) / (float) (grid.Mx - 1);
        cj = int( math.floor((lat - lowlat) / dlat) )
        ci = int( math.floor((lon - lowlon) / dlon) )
      
        if ((ci >= grid.xs) and (ci < grid.xs+grid.xm) and (cj >= grid.ys) and (cj < grid.ys+grid.ym)):
          vel = vel_ssa[ci,cj]
          cu = secpera * vel.u
          cv = secpera * vel.v
          cmag = math.sqrt(cu*cu + cv*cv)
          PISM.verbPrintf(4,PETSc.COMM_SELF,
                          " PISM%d[%3d]: lat = %7.3f, lon = %7.3f, mag = %7.2f, u = %7.2f, v = %7.2f\n",
                          grid.rank,k,clat[ci,cj],clon[ci,cj],cmag,cu,cv)
          if mask[ci,cj] == PISM.MASK_FLOATING:
            goodptcount += 1.0;
            ChiSqr += (u-cu)*(u-cu)+(v-cv)*(v-cv)
    # end with

    ChiSqr = ChiSqr / (30.0*30.0) # see page 48 of MacAyeal et al
    g_goodptcount = PISM.globalSum(goodptcount,grid.com)
    g_ChiSqr      = PISM.globalSum(ChiSqr, grid.com)
    r.printlnv( 4, """number of RIGGS data points = %d
number of RIGGS points in computed ice shelf region = %8.2f""", length , g_goodptcount);
    r.println("Chi^2 statistic for computed results compared to RIGGS is %10.3f",
                      g_ChiSqr * (156.0 / g_goodptcount))


# The main code for a run follows:
if __name__ == '__main__':
  context = PISM.Context()
  com = context.com

  PISM.set_abort_on_sigint(True)
  PISM.verbosityLevelFromOptions()
  PISM.verbPrintf(2,com,"PROSS %s (EISMINT-Ross diagnostic velocity computation mode)\n",
                  PISM.PISM_Revision) 
  PISM.stop_on_version_option()
  usage = \
"""This is a python implementation of the PISM 'pross' example.  It's usage
is identical:  

pross.py -boot_file IN.nc -Mx number -My number [-o file.nc] [-riggs file.nc]
or (at python prompt)
  run pross -boot_file IN.nc -Mx number -My number [-o file.nc] [-riggs file.nc]
where:
  -boot_file  IN.nc is input file in NetCDF format: 
              contains PISM-written model state
  -Mx         number of grid points in the x direction
  -My         number of grid points in the y direction
  -riggs      read RIGGS data from a file
notes:
  * -boot_file is required
"""

  PISM.show_usage_check_req_opts(com,"python pross.py",["-boot_file"],usage)

  config = context.config
  config.set_flag("use_ssa_velocity", True)
  config.set_flag("use_ssa_when_grounded", False)
  config.set_flag("use_constant_nuh_for_ssa", False)
  config.set("epsilon_ssafd", 0.0)

  for o in PISM.OptionsGroup(com,"","PROSS"):
    Mx = PISM.optionsInt("-Mx","Number of grid points in the X-direction",default=None)
    My = PISM.optionsInt("-My","Number of grid points in the X-direction",default=None)
    boot_file = PISM.optionsString("-boot_file","file to bootstrap from")
    riggs_file = PISM.optionsString("-riggs","file with riggs measurements")
    output_file = PISM.optionsString("-o","output file",default="ross_computed.nc")

  test_case = pross(Mx,My,boot_file,riggs_file)
  test_case.setup()
  test_case.solve()
  test_case.report()
  test_case.write(output_file)

