// Copyright (C) 2006-2009 Ed Bueler and Constantine Khroulev
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
#include "../base/grid.hh"
#include "../base/materials.hh"
#include "../base/iceModel.hh"
#include "../base/nc_util.hh"
#include "../coupler/forcing.hh"
#include "iceROSSModel.hh"


IceROSSModel::IceROSSModel(IceGrid &g)
  : IceModel(g) {  // do nothing; note derived classes must have constructors

  useSSAVelocity= PETSC_TRUE;
  computeSIAVelocities = PETSC_FALSE;
  doMassConserve = PETSC_FALSE;  // diagnostic calculation

  // further settings for velocity computation 
  useConstantNuHForSSA = PETSC_FALSE; // compute the effective viscosity in usual
           // shear-thinning way (except will extend shelf using constantNuHForSSA below,
           // also as usual)
  shelvesDragToo = PETSC_FALSE;            // exactly zero drag under shelves
  ssaEpsilon = 0.0;  // don't use this lower bound on effective viscosity

  // Ross validation uses isothermal ice by default
  iceFactory.setType(ICE_CUSTOM);
  // Most of these are defaults anyway, but this is how to go about setting defaults
  shelfExtension.setThickness(5);
  shelfExtension.setTemperature(263.15);
  // Typical strain rate is 100 m/yr per 100km in an ice shelf or fast ice stream
  shelfExtension.setStrainRate((100.0 / secpera) / (100.0 * 1.0e3));
}

PetscErrorCode IceROSSModel::createVecs() {
  PetscErrorCode ierr;

  ierr = obsAzimuth.create(grid, "azi_obs", true); CHKERRQ(ierr);
  ierr = obsAzimuth.set_attrs(NULL, "observed ice velocity azimuth",
                              "degrees_east", NULL); CHKERRQ(ierr);

  ierr = obsMagnitude.create(grid, "mag_obs", true); CHKERRQ(ierr);
  ierr = obsMagnitude.set_attrs(NULL, "observed ice velocity magnitude",
                                "m s-1", NULL); CHKERRQ(ierr);
  ierr = obsMagnitude.set_glaciological_units("m year-1"); CHKERRQ(ierr);

  ierr = obsAccurate.create(grid, "accur", true); CHKERRQ(ierr);
  ierr = obsAccurate.set_attrs(NULL, "flag for accurate observed velocity",
                               NULL, NULL); CHKERRQ(ierr);

  ierr = IceModel::createVecs(); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode IceROSSModel::destroyVecs() {
  PetscErrorCode ierr;

  ierr = obsAzimuth.destroy(); CHKERRQ(ierr);
  ierr = obsMagnitude.destroy(); CHKERRQ(ierr);
  ierr = obsAccurate.destroy(); CHKERRQ(ierr);

  ierr = IceModel::destroyVecs(); CHKERRQ(ierr);
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

  // This initializes the IceFactory and calls IceFactory.create()
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

  // update surface elev (this has to happen after
  // updateSurfaceElevationAnsMask is called in IceModel::misc_setup() above):
  ierr = verbPrintf(2,grid.com, 
     "EIS-Ross: applying floatation criterion everywhere to get smooth surface ...\n");
     CHKERRQ(ierr);
  ierr = vH.copy_to(vh); CHKERRQ(ierr);
  ierr = vh.scale(1.0 - ice->rho / ocean.rho ); CHKERRQ(ierr);

  // in preparation for SSA b.c. read; zero out vuvbar; SIA velocities will not
  //    be computed so this will stay
  ierr = vuvbar[0].set(0.0); CHKERRQ(ierr);
  ierr = vuvbar[1].set(0.0); CHKERRQ(ierr);

  // read SSA b.c. from file
  PetscTruth  ssaBCset;
  char        ssaBCfile[PETSC_MAX_PATH_LEN];
  ierr = PetscOptionsGetString(PETSC_NULL, "-ssaBC", ssaBCfile,
                               PETSC_MAX_PATH_LEN, &ssaBCset); CHKERRQ(ierr);
  if (ssaBCset == PETSC_TRUE) {
     ierr = verbPrintf(2, grid.com,
             "EIS-Ross: reading SSA boundary condition file %s and setting bdry conds\n",
             ssaBCfile); CHKERRQ(ierr);
     ierr = readShelfStreamBCFromFile(ssaBCfile); CHKERRQ(ierr);
  }

  ierr = verbPrintf(2,grid.com, "EIS-Ross: computing velocity ...\n"); CHKERRQ(ierr);
  return 0;
}

// This method is called by the pismd driver.
PetscErrorCode IceROSSModel::finishROSS() {
  PetscErrorCode  ierr;

  PetscTruth  ssaBCset;
  char        ssaBCfile[PETSC_MAX_PATH_LEN];
  ierr = PetscOptionsGetString(PETSC_NULL, "-ssaBC", ssaBCfile,
                               PETSC_MAX_PATH_LEN, &ssaBCset); CHKERRQ(ierr);
  if (ssaBCset == PETSC_FALSE)
    return 0;
    
  ierr = verbPrintf(2,grid.com, "\nEIS-Ross: reading observed velocities from %s...\n",
                    ssaBCfile); CHKERRQ(ierr);
  ierr = readObservedVels(ssaBCfile); CHKERRQ(ierr);

  ierr = computeErrorsInAccurateRegion(); CHKERRQ(ierr);
  ierr = readRIGGSandCompare(); CHKERRQ(ierr);

  // if we want to show observed velocities, they should be written 
  // to an nc file by additional code; old revisions had code to show in diagnostic
  // viewers by overwriting ubar,vbar
  
  IceModelVec2  myvNu[2];
  ierr = myvNu[0].create(grid, "myvNu", true); CHKERRQ(ierr);
  ierr = myvNu[1].create(grid, "myvNu", true); CHKERRQ(ierr);
  ierr = computeEffectiveViscosity(myvNu, ssaEpsilon); CHKERRQ(ierr);
  ierr = updateNuViewers(myvNu,myvNu,false); CHKERRQ(ierr);
  ierr = myvNu[0].destroy(); CHKERRQ(ierr);
  ierr = myvNu[1].destroy(); CHKERRQ(ierr);
  
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
  PetscScalar         **Ts;

  // fill in all temps with Ts
  IceModelVec2    *pccTs;
  if (atmosPCC != PETSC_NULL) {
    // call sets pccTs to point to IceModelVec2 with current surface temps
    ierr = atmosPCC->updateSurfTempAndProvide(
              grid.year, 0.0, (void*)(&info_atmoscoupler), pccTs); CHKERRQ(ierr);
  } else {  SETERRQ(1,"PISM ERROR: atmosPCC == PETSC_NULL  in  IceROSSModel::fillinTemps()");  }

  ierr = pccTs->get_array(Ts);  CHKERRQ(ierr);
  ierr = T3.begin_access(); CHKERRQ(ierr);
  ierr = Tb3.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      ierr = T3.setColumn(i,j,Ts[i][j]); CHKERRQ(ierr);
      ierr = Tb3.setColumn(i,j,Ts[i][j]); CHKERRQ(ierr);
    }
  }
  ierr = T3.end_access(); CHKERRQ(ierr);
  ierr = Tb3.end_access(); CHKERRQ(ierr);
  ierr = pccTs->end_access();  CHKERRQ(ierr);

  return 0;
}


PetscErrorCode IceROSSModel::readObservedVels(const char *filename) {
  PetscErrorCode  ierr;
  NCTool nc(&grid);

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
  PetscScalar  **azi, **mag, **acc, **ubar, **vbar, **H, **mask;
  
  const PetscScalar area = grid.dx * grid.dy;
  ierr = vMask.get_array(mask); CHKERRQ(ierr);
  ierr = vH.get_array(H); CHKERRQ(ierr);
  ierr = obsAzimuth.get_array(azi); CHKERRQ(ierr);    
  ierr = obsMagnitude.get_array(mag); CHKERRQ(ierr);    
  ierr = obsAccurate.get_array(acc); CHKERRQ(ierr);    
  ierr = vubar.get_array(ubar); CHKERRQ(ierr);
  ierr = vvbar.get_array(vbar); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
     if ((PismModMask(mask[i][j]) == MASK_FLOATING) && (H[i][j] > 1.0)) {
        const PetscScalar ccomputed = sqrt(PetscSqr(vbar[i][j]) + PetscSqr(ubar[i][j]));
        maxcComputed = PetscMax(maxcComputed,ccomputed);
        if (PetscAbs(acc[i][j] - 1.0) < 0.1) {
          accN += 1.0;
          accArea += area;
          const PetscScalar vobs = mag[i][j] * sin((pi/180.0) * azi[i][j]);
          const PetscScalar uobs = mag[i][j] * cos((pi/180.0) * azi[i][j]);
          // compare from readme.txt
          // uxbar(i0,j0) = magvel(i0,j0)* sin(3.1415926/180.*azvel(i0,j0))
          // uybar(i0,j0) = magvel(i0,j0)* cos(3.1415926/180.*azvel(i0,j0))
          // the difference is the fundamental transpose in IceGrid::createDA()
          // debug:
          //verbPrintf(1,grid.com,"i,j=%d,%d:  observed = (%5.4f,%5.4f),   computed =  (%5.4f,%5.4f)\n",
          //           i,j,uobs*secpera,vobs*secpera,ubar[i][j]*secpera,vbar[i][j]*secpera);
          const PetscScalar Dv = PetscAbs(vobs - vbar[i][j]);
          const PetscScalar Du = PetscAbs(uobs - ubar[i][j]);
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
  ierr = vubar.end_access(); CHKERRQ(ierr);
  ierr = vvbar.end_access(); CHKERRQ(ierr);
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

      Data1D      latdata, londata, magdata, udata, vdata;
      PetscInt    len;
      PetscScalar **ubar, **vbar, **clat, **clon, **mask;

      ierr = vLongitude.get_array(clon); CHKERRQ(ierr);
      ierr =  vLatitude.get_array(clat); CHKERRQ(ierr);
      ierr = vubar.get_array(ubar); CHKERRQ(ierr);
      ierr = vvbar.get_array(vbar); CHKERRQ(ierr);
      ierr = vMask.get_array(mask); CHKERRQ(ierr);
      
      ierr = latdata.readData(grid.com,grid.rank, riggsfile, "count", "riggslat"); CHKERRQ(ierr);
      ierr = londata.readData(grid.com,grid.rank, riggsfile, "count", "riggslon"); CHKERRQ(ierr);
      ierr = magdata.readData(grid.com,grid.rank, riggsfile, "count", "riggsmag"); CHKERRQ(ierr);
      ierr = udata.readData(grid.com,grid.rank, riggsfile, "count", "riggsu"); CHKERRQ(ierr);
      ierr = vdata.readData(grid.com,grid.rank, riggsfile, "count", "riggsv"); CHKERRQ(ierr);
      ierr = latdata.getIndexMax(&len); CHKERRQ(ierr);  // same length for all vars here
      PetscScalar  goodptcount = 0.0, ChiSqr = 0.0;
      for (PetscInt k = 0; k<len; k++) {
        PetscScalar lat, lon, mag, u, v;
        ierr = latdata.getIndexedDataValue(k, &lat); CHKERRQ(ierr);
        ierr = londata.getIndexedDataValue(k, &lon); CHKERRQ(ierr);
        ierr = magdata.getIndexedDataValue(k, &mag); CHKERRQ(ierr);
        ierr = udata.getIndexedDataValue(k, &u); CHKERRQ(ierr);
        ierr = vdata.getIndexedDataValue(k, &v); CHKERRQ(ierr);
        ierr = verbPrintf(4,grid.com,
                 " RIGGS[%3d]: lat = %7.3f, lon = %7.3f, mag = %7.2f, u = %7.2f, v = %7.2f\n",
                 k,lat,lon,mag,u,v); CHKERRQ(ierr); 
        const PetscScalar origdlat = (-5.42445 - (-12.3325)) / 110.0;
        const PetscScalar lowlat = -12.3325 - origdlat * 46.0;
        const PetscScalar dlat = (-5.42445 - lowlat) / (float) (grid.Mx - 1);        
        const PetscScalar lowlon = -5.26168;
        const PetscScalar dlon = (3.72207 - lowlon) / (float) (grid.My - 1);
        const int         ci = (int) floor((lat - lowlat) / dlat);
        const int         cj = (int) floor((lon - lowlon) / dlon);
        if ((ci >= grid.xs) && (ci < grid.xs+grid.xm) && (cj >= grid.ys) && (cj < grid.ys+grid.ym)) {
          const PetscScalar cu = secpera * vbar[ci][cj];  // note switched meaning
          const PetscScalar cv = secpera * ubar[ci][cj];
          const PetscScalar cmag = sqrt(PetscSqr(cu)+PetscSqr(cv));
          ierr = verbPrintf(4,PETSC_COMM_SELF,
                 " PISM%d[%3d]: lat = %7.3f, lon = %7.3f, mag = %7.2f, u = %7.2f, v = %7.2f\n",
                 grid.rank,k,clat[ci][cj],clon[ci][cj],cmag,cu,cv); CHKERRQ(ierr); 
          if (PismIntMask(mask[ci][cj]) == MASK_FLOATING) {
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

      ierr = vubar.end_access(); CHKERRQ(ierr);
      ierr = vvbar.end_access(); CHKERRQ(ierr);
  }
  return 0;
}

