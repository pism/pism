// Copyright (C) 2006-2010 Ed Bueler and Constantine Khroulev
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 2 of the License, or (at your option) any later
// version.
//
// PISM is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License
// along with PISM; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

#include <cstring>
#include <cstdio>
#include <cmath>
#include "../base/iceModel.hh"
#include "../base/PISMIO.hh"
#include "iceROSSModel.hh"
#include "../base/Timeseries.hh"

IceROSSModel::IceROSSModel(IceGrid &g, NCConfigVariable &conf, NCConfigVariable &conf_overrides)
  : IceModel(g, conf, conf_overrides) {  // do nothing; note derived classes must have constructors

  config.set_flag("use_ssa_velocity", true);
  computeSIAVelocities = PETSC_FALSE;

  // further settings for velocity computation 
  config.set_flag("use_constant_nuh_for_ssa", false);
  // compute the effective viscosity in usual shear-thinning way (except will
  // extend shelf using constantNuHForSSA below, also as usual)
  shelvesDragToo = PETSC_FALSE;            // exactly zero drag under shelves

  config.set("epsilon_ssa", 0.0);  // don't use this lower bound on effective viscosity

  iceFactory.setType(ICE_CUSTOM);  // ICE_CUSTOM has easy setting of ice hardness
  
  ssaStrengthExtend.set_min_thickness(5.0); // m
  const PetscReal
    DEFAULT_CONSTANT_HARDNESS_FOR_SSA = 1.9e8,  // Pa s^{1/3}; see p. 49 of MacAyeal et al 1996
    DEFAULT_TYPICAL_STRAIN_RATE = (100.0 / secpera) / (100.0 * 1.0e3),  // typical strain rate is 100 m/yr per 
    DEFAULT_nuH = ssaStrengthExtend.min_thickness_for_extension() * DEFAULT_CONSTANT_HARDNESS_FOR_SSA
                       / (2.0 * pow(DEFAULT_TYPICAL_STRAIN_RATE,2./3.)); // Pa s m
          // COMPARE: 30.0 * 1e6 * secpera = 9.45e14 is Ritz et al (2001) value of
          //          30 MPa yr for \bar\nu
  ssaStrengthExtend.set_notional_strength(DEFAULT_nuH);

}

PetscErrorCode IceROSSModel::setFromOptions() {
  PetscErrorCode ierr;

  config.set_flag("do_cold_ice_methods", true);

  ierr = IceModel::setFromOptions(); CHKERRQ(ierr);
  
  return 0;
}

PetscErrorCode IceROSSModel::createVecs() {
  PetscErrorCode ierr;

  ierr = obsAzimuth.create(grid, "azi_obs", true); CHKERRQ(ierr);
  ierr = obsAzimuth.set_attrs("", "observed ice velocity azimuth",
                              "degrees_east", ""); CHKERRQ(ierr);

  ierr = obsMagnitude.create(grid, "mag_obs", true); CHKERRQ(ierr);
  ierr = obsMagnitude.set_attrs("", "observed ice velocity magnitude",
                                "m s-1", ""); CHKERRQ(ierr);
  ierr = obsMagnitude.set_glaciological_units("m year-1"); CHKERRQ(ierr);

  ierr = obsAccurate.create(grid, "accur", true); CHKERRQ(ierr);
  ierr = obsAccurate.set_attrs("", "flag for accurate observed velocity",
                               "", ""); CHKERRQ(ierr);

  ierr = IceModel::createVecs(); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode IceROSSModel::set_vars_from_options() {
  PetscErrorCode ierr;

  ierr = verbPrintf(2,grid.com, 
    "initializing EISMINT-Ross ice shelf velocity computation ... \n");
  CHKERRQ(ierr);

  // This reads the -boot_from option and does the bootstrapping:
  ierr = IceModel::set_vars_from_options(); CHKERRQ(ierr);

  // fill in temperatures at depth according to special rule: temp in column
  // equals temp at surface
  ierr = fillinTemps();  CHKERRQ(ierr);

  return 0;
}


PetscErrorCode IceROSSModel::init_physics() {
  PetscErrorCode ierr;

  // This initializes the IceFlowLawFactory and calls IceFlowLawFactory.create()
  ierr = IceModel::init_physics(); CHKERRQ(ierr);

  CustomGlenIce *cgi = dynamic_cast<CustomGlenIce*>(ice);
  if (cgi) {
    ierr = cgi->setHardness(1.9e8);CHKERRQ(ierr); // Pa s^{1/3}; (MacAyeal et al 1996) value
    // This sets a default hardness parameter if CustomGlenIce is chosen (the default we set in the constructor).
    // The user will still be able to override it without a warning
  }
  ierr = ice->setFromOptions();CHKERRQ(ierr);
  if (!cgi) {
    ierr = verbPrintf(2,grid.com,
                      "WARNING: Not using CustomGlenIce so cannot set hardness defaults\n"
                      "         (Perhaps you chose something else with -ice_type xxx)\n"
                      "         Details on your chosen ice follows\n"); CHKERRQ(ierr);
    ierr = ice->printInfo(2);CHKERRQ(ierr);
  }

  return 0;
}


PetscErrorCode IceROSSModel::misc_setup() {
  PetscErrorCode  ierr;

  ierr = IceModel::misc_setup(); CHKERRQ(ierr);

  double ocean_rho = config.get("sea_water_density");

  // update surface elev (this has to happen after
  // updateSurfaceElevationAnsMask is called in IceModel::misc_setup() above):
  ierr = verbPrintf(2,grid.com, 
     "EIS-Ross: applying floatation criterion everywhere to get smooth surface ...\n");
     CHKERRQ(ierr);
  ierr = vH.copy_to(vh); CHKERRQ(ierr);
  ierr = vh.scale(1.0 - ice->rho / ocean_rho ); CHKERRQ(ierr);

  // in preparation for SSA b.c. read; zero out uvbar; SIA velocities will not
  //    be computed so this will stay
  ierr = uvbar.set(0.0); CHKERRQ(ierr);

  // read SSA b.c. from file
  bool   ssaBCset;
  string ssaBCfile;
  ierr = PetscOptionsBegin(grid.com, "", "IceROSSModel options", ""); CHKERRQ(ierr);
  {
    ierr = PISMOptionsString("-ssaBC", "Specifies the file containing SSA boundary conditions",
			     ssaBCfile, ssaBCset); CHKERRQ(ierr);
  }
  ierr = PetscOptionsEnd(); CHKERRQ(ierr);

  if (ssaBCset) {
     ierr = verbPrintf(2, grid.com,
             "EIS-Ross: reading SSA boundary condition file %s and setting bdry conds\n",
		       ssaBCfile.c_str()); CHKERRQ(ierr);
     ierr = readShelfStreamBCFromFile(ssaBCfile.c_str()); CHKERRQ(ierr);
     config.set_flag("use_ssabc", true);
     config.set_string("ssabc_filename", ssaBCfile);
  } else {
     config.set_flag("use_ssabc", false);
  }

  ierr = verbPrintf(2,grid.com, "EIS-Ross: computing velocity ...\n"); CHKERRQ(ierr);
  return 0;
}

// This method is called by the pismd driver.
PetscErrorCode IceROSSModel::finishROSS() {
  PetscErrorCode  ierr;

  if (config.get_flag("use_ssabc") == false)
    return 0;

  double ssaEpsilon = config.get("epsilon_ssa");
  string ssaBCfile = config.get_string("ssabc_filename");

  ierr = verbPrintf(2,grid.com, "\nEIS-Ross: reading observed velocities from %s...\n",
                    ssaBCfile.c_str()); CHKERRQ(ierr);
  ierr = readObservedVels(ssaBCfile.c_str()); CHKERRQ(ierr);

  ierr = computeErrorsInAccurateRegion(); CHKERRQ(ierr);
  ierr = readRIGGSandCompare(); CHKERRQ(ierr);

  // if we want to show observed velocities, they should be written 
  // to an nc file by additional code; old revisions had code to show in diagnostic
  // viewers by overwriting ubar,vbar
  
  IceModelVec2S  myvNu[2];
  ierr = myvNu[0].create(grid, "myvNu", true); CHKERRQ(ierr);
  ierr = myvNu[1].create(grid, "myvNu", true); CHKERRQ(ierr);
  ierr = computeEffectiveViscosity(myvNu, ssaEpsilon); CHKERRQ(ierr);
  ierr = update_nu_viewers(myvNu,myvNu,false); CHKERRQ(ierr);
  
  PetscInt    pause_time = 0;
  ierr = PetscOptionsGetInt(PETSC_NULL, "-pause", &pause_time, PETSC_NULL); CHKERRQ(ierr);
  if (pause_time > 0) {
    ierr = verbPrintf(2,grid.com,"\npausing for %d secs ...\n",pause_time); CHKERRQ(ierr);
    ierr = PetscSleep(pause_time); CHKERRQ(ierr);
  }
  
  return 0;
}


PetscErrorCode IceROSSModel::fillinTemps() {
  PetscErrorCode      ierr;

  // fill in all temps with artm:
  if (surface != PETSC_NULL) {
    ierr = surface->ice_surface_temperature(grid.year, 0.0, artm); CHKERRQ(ierr);
  } else {
    SETERRQ(1,"PISM ERROR: surface == PETSC_NULL");
  }

  ierr = artm.begin_access();  CHKERRQ(ierr);
  ierr = T3.begin_access(); CHKERRQ(ierr);
  ierr = Tb3.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      ierr =  T3.setColumn(i, j, artm(i,j)); CHKERRQ(ierr);
      ierr = Tb3.setColumn(i, j, artm(i,j)); CHKERRQ(ierr);
    }
  }
  ierr = T3.end_access(); CHKERRQ(ierr);
  ierr = Tb3.end_access(); CHKERRQ(ierr);
  ierr = artm.end_access();  CHKERRQ(ierr);

  return 0;
}


PetscErrorCode IceROSSModel::readObservedVels(const char *filename) {
  PetscErrorCode  ierr;
  PISMIO nc(&grid);

  ierr = nc.open_for_reading(filename); CHKERRQ(ierr);

  // will create "local interpolation context" from dimensions, limits, and
  //   lengths extracted from bootstrap file and from information about the
  //   part of the grid owned by this processor;
  grid_info g;
  ierr = nc.get_grid_info_2d(g); CHKERRQ(ierr);
  ierr = nc.close(); CHKERRQ(ierr);

  LocalInterpCtx lic(g, NULL, NULL, grid); // 2D only

  ierr =  obsAccurate.regrid(filename, lic, true); CHKERRQ(ierr);
  ierr = obsMagnitude.regrid(filename, lic, true); CHKERRQ(ierr);
  ierr =   obsAzimuth.regrid(filename, lic, true); CHKERRQ(ierr);
  
  return 0;
}


PetscErrorCode IceROSSModel::computeErrorsInAccurateRegion() {
  // average over grid, where observed velocity is accurate *according to
  // 111by147grid.dat*, the difference between computed and observed u,v
  PetscErrorCode  ierr;
  PetscScalar  uerr=0.0, verr=0.0, relvecerr=0.0, accN=0.0, 
               accArea=0.0, maxcComputed=0.0, vecErrAcc = 0.0;
  PetscScalar  **azi, **mag, **acc, **H;
  
  const PetscScalar area = grid.dx * grid.dy;
  ierr = vMask.begin_access(); CHKERRQ(ierr);
  ierr = vH.get_array(H); CHKERRQ(ierr);
  ierr = obsAzimuth.get_array(azi); CHKERRQ(ierr);    
  ierr = obsMagnitude.get_array(mag); CHKERRQ(ierr);    
  ierr = obsAccurate.get_array(acc); CHKERRQ(ierr);    
  ierr = vel_bar.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (vMask.is_floating(i,j) && (H[i][j] > 1.0)) {
        const PetscScalar ccomputed = sqrt(PetscSqr(vel_bar(i,j).v) + PetscSqr(vel_bar(i,j).u));
        maxcComputed = PetscMax(maxcComputed,ccomputed);
        if (PetscAbs(acc[i][j] - 1.0) < 0.1) {
          accN += 1.0;
          accArea += area;
          const PetscScalar uobs = mag[i][j] * sin((pi/180.0) * azi[i][j]);
          const PetscScalar vobs = mag[i][j] * cos((pi/180.0) * azi[i][j]);
          // compare from readme.txt
          // uxbar(i0,j0) = magvel(i0,j0)* sin(3.1415926/180.*azvel(i0,j0))
          // uybar(i0,j0) = magvel(i0,j0)* cos(3.1415926/180.*azvel(i0,j0))
          // debug:
          //verbPrintf(1,grid.com,"i,j=%d,%d:  observed = (%5.4f,%5.4f),   computed =  (%5.4f,%5.4f)\n",
          //           i,j,uobs*secpera,vobs*secpera,vel_bar(i,j).u*secpera,vel_bar(i,j).v*secpera);
          const PetscScalar Dv = PetscAbs(vobs - vel_bar(i,j).v);
          const PetscScalar Du = PetscAbs(uobs - vel_bar(i,j).u);
          verr += Dv;
          uerr += Du;
          relvecerr += (PetscSqr(Dv) + PetscSqr(Du)) / (PetscSqr(vobs) + PetscSqr(uobs));
          vecErrAcc += (PetscSqr(Dv) + PetscSqr(Du)) * area;
        }
      }
    }
  }
  ierr = obsMagnitude.end_access(); CHKERRQ(ierr);
  ierr = obsAzimuth.end_access(); CHKERRQ(ierr);
  ierr = obsAccurate.end_access(); CHKERRQ(ierr);
  ierr = vel_bar.end_access(); CHKERRQ(ierr);
  ierr = vMask.end_access(); CHKERRQ(ierr);
  ierr = vH.end_access(); CHKERRQ(ierr);

  PetscScalar  guerr, gverr, grelvecerr, gaccN, 
               gaccArea, gmaxcComputed, gvecErrAcc;
  ierr = PetscGlobalMax(&maxcComputed, &gmaxcComputed, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalSum(&accN, &gaccN, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalSum(&accArea, &gaccArea, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalSum(&verr, &gverr, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalSum(&uerr, &guerr, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalSum(&relvecerr, &grelvecerr, grid.com); CHKERRQ(ierr);
  ierr = PetscGlobalSum(&vecErrAcc, &gvecErrAcc, grid.com); CHKERRQ(ierr);

  ierr = verbPrintf(2,grid.com,"maximum computed speed in ice shelf is %10.3f (m/a)\n",
             gmaxcComputed * secpera); CHKERRQ(ierr);
  ierr = verbPrintf(2,grid.com,"ERRORS relative to observations of Ross Ice Shelf:\n");
             CHKERRQ(ierr);
  ierr = verbPrintf(2,grid.com,
             "  [number of grid points in 'accurate observed area' = %d]\n",
             (int) gaccN); CHKERRQ(ierr);
  ierr = verbPrintf(2,grid.com,
             "  [area of 'accurate observed area' = %9.4f (km^2)]\n",
             gaccArea / 1e6); CHKERRQ(ierr);
  ierr = verbPrintf(2,grid.com,
             "  following are average errors computed over 'accurate observed area':\n");
             CHKERRQ(ierr);
  ierr = verbPrintf(2,grid.com,
             "  average error in x-comp of vel       = %9.3f (m/a)\n",
             (gverr * secpera) / gaccN); CHKERRQ(ierr);
  ierr = verbPrintf(2,grid.com,
             "  average error in y-comp of vel       = %9.3f (m/a)\n",
             (guerr * secpera) / gaccN); CHKERRQ(ierr);
  ierr = verbPrintf(2,grid.com,
             "  average relative error in vector vel = %9.5f\n",
             grelvecerr / gaccN); CHKERRQ(ierr);
  gvecErrAcc = secpera * sqrt(gvecErrAcc) / sqrt(gaccArea);
  ierr = verbPrintf(2,grid.com,
             "  rms average error in vector vel      = %9.3f (m/a)\n",
             gvecErrAcc); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceROSSModel::readRIGGSandCompare() {
  PetscErrorCode  ierr;
  PetscTruth      riggsSet;
  char            riggsfile[PETSC_MAX_PATH_LEN];

  ierr = PetscOptionsGetString(PETSC_NULL, "-riggs", riggsfile,
                               PETSC_MAX_PATH_LEN, &riggsSet); CHKERRQ(ierr);
  if (riggsSet == PETSC_FALSE) {
    return 0;
  } else {
      ierr = verbPrintf(2,grid.com,"comparing to RIGGS data in %s ...\n",
             riggsfile); CHKERRQ(ierr);

      Timeseries latdata(&grid, "riggslat", "count"),
	londata(&grid, "riggslon", "count"),
	magdata(&grid, "riggsmag", "count"),
	udata(&grid, "riggsu", "count"),
	vdata(&grid, "riggsv", "count");
      PetscInt    len;
      PetscScalar **clat, **clon;

      ierr = vLongitude.get_array(clon); CHKERRQ(ierr);
      ierr =  vLatitude.get_array(clat); CHKERRQ(ierr);
      ierr = vel_bar.begin_access(); CHKERRQ(ierr);
      ierr = vMask.begin_access();  CHKERRQ(ierr);

      ierr = magdata.set_units("m year-1", ""); CHKERRQ(ierr);
      ierr =   udata.set_units("m year-1", ""); CHKERRQ(ierr);
      ierr =   vdata.set_units("m year-1", ""); CHKERRQ(ierr);

      ierr = latdata.read(riggsfile); CHKERRQ(ierr);
      ierr = londata.read(riggsfile); CHKERRQ(ierr);
      ierr = magdata.read(riggsfile); CHKERRQ(ierr);
      ierr =   udata.read(riggsfile); CHKERRQ(ierr);
      ierr =   vdata.read(riggsfile); CHKERRQ(ierr);
      
      // same length for all vars here
      len = latdata.length();
      PetscScalar  goodptcount = 0.0, ChiSqr = 0.0;
      for (PetscInt k = 0; k<len; k++) {
        PetscScalar lat, lon, mag, u, v;
        lat = latdata[k];
	lon = londata[k];
	mag = magdata[k];
	u   = udata[k];
	v   = vdata[k];
        ierr = verbPrintf(4,grid.com,
                 " RIGGS[%3d]: lat = %7.3f, lon = %7.3f, mag = %7.2f, u = %7.2f, v = %7.2f\n",
                 k,lat,lon,mag,u,v); CHKERRQ(ierr); 
        const PetscScalar origdlat = (-5.42445 - (-12.3325)) / 110.0;
        const PetscScalar lowlat = -12.3325 - origdlat * 46.0;
        const PetscScalar dlat = (-5.42445 - lowlat) / (float) (grid.My - 1);        
        const PetscScalar lowlon = -5.26168;
        const PetscScalar dlon = (3.72207 - lowlon) / (float) (grid.Mx - 1);
        const int         cj = (int) floor((lat - lowlat) / dlat);
        const int         ci = (int) floor((lon - lowlon) / dlon);
        if ((ci >= grid.xs) && (ci < grid.xs+grid.xm) && (cj >= grid.ys) && (cj < grid.ys+grid.ym)) {
          const PetscScalar cu = secpera * vel_bar(ci,cj).u;
          const PetscScalar cv = secpera * vel_bar(ci,cj).v;
          const PetscScalar cmag = sqrt(PetscSqr(cu)+PetscSqr(cv));
          ierr = verbPrintf(4,PETSC_COMM_SELF,
                 " PISM%d[%3d]: lat = %7.3f, lon = %7.3f, mag = %7.2f, u = %7.2f, v = %7.2f\n",
                 grid.rank,k,clat[ci][cj],clon[ci][cj],cmag,cu,cv); CHKERRQ(ierr); 
          if (vMask.value(ci,cj) == MASK_FLOATING) {
            goodptcount += 1.0;
            ChiSqr += PetscSqr(u-cu)+PetscSqr(v-cv);
          }
        }
      }
      ChiSqr = ChiSqr / PetscSqr(30.0); // see page 48 of MacAyeal et al
      PetscScalar g_goodptcount, g_ChiSqr;
      ierr = PetscGlobalSum(&goodptcount, &g_goodptcount, grid.com); CHKERRQ(ierr);
      ierr = PetscGlobalSum(&ChiSqr, &g_ChiSqr, grid.com); CHKERRQ(ierr);
      ierr = verbPrintf(4,grid.com,"number of RIGGS data points = %d\n"
             "number of RIGGS points in computed ice shelf region = %8.2f\n",
             len, g_goodptcount); CHKERRQ(ierr);
      ierr = verbPrintf(2,grid.com,"Chi^2 statistic for computed results compared to RIGGS is %10.3f\n",
             g_ChiSqr * (156.0 / g_goodptcount)); CHKERRQ(ierr);

      ierr = vLongitude.end_access(); CHKERRQ(ierr);
      ierr =  vLatitude.end_access(); CHKERRQ(ierr);
      ierr = vMask.end_access(); CHKERRQ(ierr);

      ierr = vel_bar.end_access(); CHKERRQ(ierr);
  }
  return 0;
}

