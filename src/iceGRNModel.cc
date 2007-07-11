// Copyright (C) 2004-2007 Jed Brown and Ed Bueler
//
// This file is part of Pism.
//
// Pism is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 2 of the License, or (at your option) any later
// version.
//
// Pism is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License
// along with Pism; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

#include <petscvec.h>
#include <petscda.h>
#include <cstring>
#include <netcdf.h>
#include "nc_util.hh"
#include "grid.hh"
#include "materials.hh"
#include "iceModel.hh"
#include "iceGRNModel.hh"

const PetscScalar DEFAULT_ABLATION_IN_OCEAN0 = 20.0;   // m/a
const PetscScalar EISMINT_G_geothermal   = 0.050;      // J/m^2 s; geo. heat flux
const int CONST_PIECE_FWD_INTERP = 0;
const int CONST_PIECE_BCK_INTERP = 1;
const int LINEAR_INTERP = 2;


IceGRNModel::IceGRNModel(IceGrid &g, IceType &i) : IceModel(g, i) {
  // only call parent's constructor; do all classes need constructors?
}


PetscErrorCode IceGRNModel::setFromOptions() {
  PetscErrorCode ierr;
  PetscTruth ssl2Set, ssl3Set, ccl3Set, gwl3Set, gkSet, flowlawSet;
  PetscInt lawNum;

  testnum = 1;  // SSL2 is the default
  ierr = PetscOptionsHasName(PETSC_NULL, "-ssl2", &ssl2Set); CHKERRQ(ierr);
  if (ssl2Set == PETSC_TRUE)   testnum = 1;
  ierr = PetscOptionsHasName(PETSC_NULL, "-ssl3", &ssl3Set); CHKERRQ(ierr);
  if (ssl3Set == PETSC_TRUE)   testnum = 2;
  ierr = PetscOptionsHasName(PETSC_NULL, "-ccl3", &ccl3Set); CHKERRQ(ierr);
  if (ccl3Set == PETSC_TRUE)   testnum = 3;
  ierr = PetscOptionsHasName(PETSC_NULL, "-gwl3", &gwl3Set); CHKERRQ(ierr);
  if (gwl3Set == PETSC_TRUE)   testnum = 4;

  ierr = PetscOptionsHasName(PETSC_NULL, "-gk", &gkSet); CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(PETSC_NULL, "-law", &lawNum, &flowlawSet); CHKERRQ(ierr);

  if (ssl3Set == PETSC_TRUE) {
    enhancementFactor = 1;
    if (gkSet == PETSC_FALSE && ((flowlawSet == PETSC_TRUE && lawNum != 4) || flowlawSet == PETSC_FALSE)) {
      verbPrintf(1, grid.com, "WARNING: SSL3 specified, but not -gk\n");
    }
  } else {
    enhancementFactor = 3;
  }
  ierr = IceModel::setFromOptions(); CHKERRQ(ierr);  // note: user value for -e will override e=3
  return 0;
}


PetscErrorCode IceGRNModel::initFromOptions() {
  PetscErrorCode ierr;
  int stat;

  ierr = IceModel::initFromOptions(); CHKERRQ(ierr);
  ierr = saveOrigVecs(); CHKERRQ(ierr);

  char inFile[PETSC_MAX_PATH_LEN], coreFile[PETSC_MAX_PATH_LEN];
  int usrYear;
  PetscTruth inFileSet, bootFileSet, forceSet, nopddSet, yearSet, flowLawSet;
  
  verbPrintf(3, grid.com,"geothermal flux vGhf is being set to: %f\n",
             EISMINT_G_geothermal);
  ierr = VecSet(vGhf, EISMINT_G_geothermal); CHKERRQ(ierr);
  

  ierr = PetscOptionsGetString(PETSC_NULL, "-if", inFile,
                               PETSC_MAX_PATH_LEN, &inFileSet); CHKERRQ(ierr);

  ierr = PetscOptionsGetString(PETSC_NULL, "-bif", inFile,
                               PETSC_MAX_PATH_LEN, &bootFileSet); CHKERRQ(ierr);

  ierr = PetscOptionsGetString(PETSC_NULL, "-forcing", coreFile,
                               PETSC_MAX_PATH_LEN, &forceSet); CHKERRQ(ierr);

  ierr = PetscOptionsGetInt(PETSC_NULL, "-ys", &usrYear,
                            &yearSet); CHKERRQ(ierr);

  ierr = PetscOptionsHasName(PETSC_NULL, "-no_pdd", &nopddSet); CHKERRQ(ierr);

  ierr = PetscOptionsHasName(PETSC_NULL, "-law", &flowLawSet); CHKERRQ(ierr);

  if (yearSet == PETSC_TRUE) {
    grid.p->year = usrYear;
  }
 
  if (nopddSet == PETSC_TRUE) {
    doPDD = PETSC_FALSE;
  } else { 
    // in this case, turn the PDD on for this derived class, so no "-pdd" is needed to turn it on
    doPDD = PETSC_TRUE;
    if (pddStuffCreated == PETSC_FALSE) {
      ierr = initPDDFromOptions(); CHKERRQ(ierr);
    }
    PetscTruth pddSummerWarmingSet, pddStdDevSet;
    ierr = PetscOptionsHasName(PETSC_NULL, "-pdd_summer_warming", &pddSummerWarmingSet); CHKERRQ(ierr);
    if (pddSummerWarmingSet == PETSC_TRUE) { // note IceGRNModel::getSummerWarming() is below
      verbPrintf(1, grid.com, 
         "WARNING: -pdd_summer_warming option ignored.  Using EISMINT-GREENLAND\n"
         "         summer temperature formula\n");
    }
    ierr = PetscOptionsHasName(PETSC_NULL, "-pdd_std_dev", &pddStdDevSet); CHKERRQ(ierr);
    if (pddStdDevSet == PETSC_FALSE) {
      pddStdDev = 5.0;  // EISMINT-GREENLAND default; note we allow user to override
    }
  }
  
  if (testnum == 2 && flowLawSet == PETSC_TRUE) {
    verbPrintf(1, grid.com, "WARNING: using flow law specified, not gk for SSL3\n");
  }

  if (inFileSet == PETSC_TRUE) {
    if (bootFileSet) {
      verbPrintf(1, grid.com, "WARNING: -bif and -if given; using -if\n");
    }
  } else if (bootFileSet == PETSC_TRUE) {
    // after we set the new temperatures, we need
    // to set the 3D temps again
    ierr = updateTs(); CHKERRQ(ierr);
    ierr = putTempAtDepth(); CHKERRQ(ierr);
    ierr = cleanExtraLand(); CHKERRQ(ierr);
  } else {
    SETERRQ(1, "ERROR: IceGRNModel needs an input file\n");
  }

  if (forceSet == PETSC_TRUE) {
    int ncid, v_dT, v_dSea, v_timeT, v_timeSea;

    iceCoreIdxT = 0;
    iceCoreIdxSea = 0;
    // read in core data from file
    if (grid.rank == 0) {
      stat = nc_open(coreFile, 0, &ncid); CHKERRQ(nc_check(stat));
      stat = nc_inq_varid(ncid, "delta_T", &v_dT); CHKERRQ(nc_check(stat));
      stat = nc_inq_varid(ncid, "delta_Sea", &v_dSea); CHKERRQ(nc_check(stat));
      stat = nc_inq_varid(ncid, "delta_T_t", &v_timeT); CHKERRQ(nc_check(stat));
      stat = nc_inq_varid(ncid, "delta_Sea_t", &v_timeSea); CHKERRQ(nc_check(stat));
      ierr = getInterpolationCode(ncid, v_dT, &gripDeltaTInterp); CHKERRQ(ierr);
      ierr = getInterpolationCode(ncid, v_dSea, &gripDeltaSeaInterp); CHKERRQ(ierr);
    }
    MPI_Bcast(&gripDeltaTInterp, 1, MPI_INT, 0, grid.com);
    MPI_Bcast(&gripDeltaSeaInterp, 1, MPI_INT, 0, grid.com);
    ierr = ncVarBcastVec(ncid, v_dT, &vIceCoreDeltaT); CHKERRQ(ierr);  // creates this Vec
    ierr = ncVarBcastVec(ncid, v_dSea, &vIceCoreDeltaSea); CHKERRQ(ierr);  // creates this Vec
    ierr = ncVarBcastVec(ncid, v_timeT, &vIceCoreTimeT); CHKERRQ(ierr);  // creates this Vec
    ierr = ncVarBcastVec(ncid, v_timeSea, &vIceCoreTimeSea); CHKERRQ(ierr);  // creates this Vec
   
    ierr = initGRIPVecIdx(); CHKERRQ(ierr);

  } else if (testnum == 3) {
    SETERRQ(1, "ERROR: EISMINT-GREENLAND experiment CCL3 needs ice core data\n");
  }

  if (testnum == 4 && forceSet == PETSC_TRUE) {
    SETERRQ(1, "ERROR: EISMINT-GREENLAND experiment GWL3 cannot use climate forcing option.\n");
  }

  if (!isInitialized()) {
    SETERRQ(1, "Model has not been initialized.\n");
  }
  return 0;
}

PetscErrorCode IceGRNModel::initGRIPVecIdx() {
  PetscErrorCode ierr;
  PetscScalar *yearsT, *yearsSea;

  ierr = VecGetArray(vIceCoreTimeT, &yearsT); CHKERRQ(ierr);
  ierr = VecGetArray(vIceCoreTimeSea, &yearsSea); CHKERRQ(ierr);

  ierr = VecGetLocalSize(vIceCoreTimeT, &iceCoreLenT); CHKERRQ(ierr);
  ierr = VecGetLocalSize(vIceCoreTimeSea, &iceCoreLenSea); CHKERRQ(ierr);

  int r, l=0;
  r = iceCoreLenT;
  // do a binary search to find where our year fits in.
  while (r > l + 1) {
    PetscInt j = (r + l)/2;
    if(grid.p->year < yearsT[j]) {
      r = j;
    } else {
      l = j;
    }
  }    
  iceCoreIdxT = l;
  // maybe we are already past our place.
  if (l >= iceCoreLenT) {
    climateForcing = PETSC_FALSE;
    verbPrintf(1, grid.com, "No climate forcing data for temperature. Using last DeltaTs and DeltaSea.\n");
  } else {
    climateForcing = PETSC_TRUE;
  }

  l=0;
  r = iceCoreLenSea;
  // do a binary search to find where our year fits in.
  while (r > l + 1) {
    PetscInt j = (r + l)/2;
    if(grid.p->year < yearsSea[j]) {
      r = j;
    } else {
      l = j;
    }
  }    
  iceCoreIdxSea = l;

  // maybe we are already past our place.
  if (l >= iceCoreLenSea) {
    climateForcing = PETSC_FALSE;
    verbPrintf(1, grid.com, "No climate forcing data for sea level. Using last DeltaTs and DeltaSea.\n");
  } else {
    climateForcing = PETSC_TRUE;
  }

  ierr = VecRestoreArray(vIceCoreTimeT, &yearsT); CHKERRQ(ierr);
  ierr = VecRestoreArray(vIceCoreTimeSea, &yearsSea); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode IceGRNModel::createVecs() {
  PetscErrorCode ierr;

  ierr = IceModel::createVecs(); CHKERRQ(ierr);

  ierr = VecDuplicate(vh, &vOrigBed); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceGRNModel::destroyVecs() {
  PetscErrorCode ierr;

  ierr = IceModel::destroyVecs(); CHKERRQ(ierr);

  ierr = VecDestroy(vIceCoreDeltaT); CHKERRQ(ierr);
  ierr = VecDestroy(vIceCoreDeltaSea); CHKERRQ(ierr);
  ierr = VecDestroy(vIceCoreTimeT); CHKERRQ(ierr);
  ierr = VecDestroy(vIceCoreTimeSea); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceGRNModel::copyOrigBed() {
  PetscErrorCode ierr;
  ierr = VecCopy(vOrigBed, vbed); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceGRNModel::getInterpolationCode(int ncid, int vid, int *code) {
  int stat;
  char attr[NC_MAX_NAME+1];
  size_t len;

  stat = nc_get_att_text(ncid, vid, "interpolation", attr); CHKERRQ(nc_check(stat));
  stat = nc_inq_attlen(ncid, vid, "interpolation", &len); CHKERRQ(nc_check(stat));

  attr[len] = '\0';

  if (strcmp(attr, "constant_piecewise_forward") == 0) {
    *code = CONST_PIECE_FWD_INTERP;
  } else if (strcmp(attr, "constant_piecewise_backward") == 0) {
    *code = CONST_PIECE_BCK_INTERP;
  } else if (strcmp(attr, "linear") == 0) {
    *code = LINEAR_INTERP;
  } else {
    verbPrintf(1, grid.com, "Interpolation '%s' is unknown, defaulting to linear.\n", attr);
    *code = LINEAR_INTERP;
  }
  return 0;
}


PetscErrorCode IceGRNModel::saveOrigVecs() {
  PetscErrorCode ierr;

  // we need to save the original bed information and copy
  // it back before we write it to the nc file.
  ierr = VecCopy(vbed, vOrigBed); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode IceGRNModel::additionalAtStartTimestep() {
    PetscErrorCode ierr;
    PetscScalar *deltaT, *iceCoreTimeT, *iceCoreTimeSea, *deltaSea;
    PetscScalar **Ts, **bed, **origBed;
    PetscScalar TsChange, seaChange;

    // at the beginning of each time step we need to recompute
    // surface temperatures from surface elevation and latitude
    ierr = updateTs(); CHKERRQ(ierr);

    // apply global warming temperature forcing
    if (testnum == 4) {
      PetscScalar t_increase;
      PetscScalar age = grid.p->year - startYear;
      if (age <= 80) {
        t_increase = (age) * .035;
      } else if (age <= 500) {
        t_increase = 2.8 + (age - 80) * .0017;
      } else {
        t_increase = 3.514;
      }
      ierr = DAVecGetArray(grid.da2, vTs, &Ts); CHKERRQ(ierr);

      //verbPrintf(2, grid.com, "age: %f, t_increase: %f\n", age, t_increase);

      for (PetscInt i = grid.xs; i<grid.xs+grid.xm; ++i) {
        for (PetscInt j = grid.ys; j<grid.ys+grid.ym; ++j) {
          Ts[i][j] += t_increase;
        }
      }

      ierr = DAVecRestoreArray(grid.da2, vTs, &Ts); CHKERRQ(ierr);

      ierr = DALocalToLocalBegin(grid.da2, vTs, INSERT_VALUES, vTs); CHKERRQ(ierr);
      ierr = DALocalToLocalEnd(grid.da2, vTs, INSERT_VALUES, vTs); CHKERRQ(ierr);

    }
    
    if (testnum == 3 && climateForcing == PETSC_TRUE) {

      // apply the climate forcing
      ierr = VecGetArray(vIceCoreDeltaT, &deltaT); CHKERRQ(ierr);
      ierr = VecGetArray(vIceCoreDeltaSea, &deltaSea); CHKERRQ(ierr);
      ierr = VecGetArray(vIceCoreTimeT, &iceCoreTimeT); CHKERRQ(ierr);
      ierr = VecGetArray(vIceCoreTimeSea, &iceCoreTimeSea); CHKERRQ(ierr);
      ierr = DAVecGetArray(grid.da2, vTs, &Ts); CHKERRQ(ierr);
      ierr = DAVecGetArray(grid.da2, vbed, &bed); CHKERRQ(ierr);
      ierr = DAVecGetArray(grid.da2, vOrigBed, &origBed); CHKERRQ(ierr);

      // if there was a large time step, it is possible
      // that we skip over multiple entries
      while (iceCoreIdxT < iceCoreLenT && grid.p->year > iceCoreTimeT[iceCoreIdxT]) {
         iceCoreIdxT++;
      }
      while (iceCoreIdxSea < iceCoreLenSea && grid.p->year > iceCoreTimeSea[iceCoreIdxSea]) {
         iceCoreIdxSea++;
      }
      if (iceCoreIdxT >= iceCoreLenT || iceCoreIdxSea >= iceCoreLenSea) {
        verbPrintf(1, grid.com, "ATTENTION: no more climate data. Turning off Climate Forcing.\n");
        climateForcing = PETSC_FALSE;
      } else {

        // if we have exact data, use it
        if (grid.p->year == iceCoreTimeT[iceCoreIdxT]) {
          TsChange = deltaT[iceCoreIdxT];
        } else { // otherwise we need to interpolate
          PetscScalar y0, y1;
          switch (gripDeltaTInterp) {
            case CONST_PIECE_BCK_INTERP:
              // use the data point we are infront of
              if (iceCoreIdxT == 0) {
                TsChange = 0;
                verbPrintf(2, grid.com, "model year precedes beginning of climate Data; using TsChange=0\n");
              } else {
                TsChange = deltaT[iceCoreIdxT-1];
              }
              break;
            case CONST_PIECE_FWD_INTERP:
              TsChange = deltaT[iceCoreIdxT];
              break;
            case LINEAR_INTERP:
              if (iceCoreIdxT == 0) {
                TsChange = 0;
                verbPrintf(2, grid.com, "model year precedes beginning of climate Data; using TsChange=0\n");
              } else {
                y0 = deltaT[iceCoreIdxT-1];
                y1 = deltaT[iceCoreIdxT];
                TsChange = y0 + ((y1-y0) / (iceCoreTimeT[iceCoreIdxT]-iceCoreTimeT[iceCoreIdxT-1]))
                                    * (grid.p->year-iceCoreTimeT[iceCoreIdxT-1]);
              }
              break;
            default:
              SETERRQ(1, "Unknown Interpolation method");
          }
        }
        if (grid.p->year == iceCoreTimeSea[iceCoreIdxSea]) {
          seaChange = deltaSea[iceCoreIdxSea];
        } else { // otherwise, we need to interpolate.
          PetscScalar y0, y1;
          switch (gripDeltaSeaInterp) {
            case CONST_PIECE_BCK_INTERP:
              // use the data point we are infront of
              if (iceCoreIdxSea == 0) {
                seaChange = 0;
                verbPrintf(2, grid.com, "model year precedes beginning of climate Data; using seaChange=0\n");
              } else {
                seaChange = deltaSea[iceCoreIdxSea-1];
              }
              break;
            case CONST_PIECE_FWD_INTERP:
              seaChange = deltaSea[iceCoreIdxSea];
              break;
            case LINEAR_INTERP:
              if (iceCoreIdxSea == 0) {
                seaChange = 0;
                verbPrintf(2, grid.com, "model year precedes beginning of climate Data; using seaChange=0\n");
              } else {
                y0 = deltaSea[iceCoreIdxSea-1];
                y1 = deltaSea[iceCoreIdxSea];
                seaChange = y0 + ((y1-y0) / (iceCoreTimeSea[iceCoreIdxSea]-iceCoreTimeSea[iceCoreIdxSea-1]))
                                    * (grid.p->year-iceCoreTimeSea[iceCoreIdxSea-1]);
              }
              break;
            default:
              SETERRQ(1, "Unknown Interpolation method");
          }
        }

        //verbPrintf(2, grid.com, "For year: %f, TsChange: %f\n", grid.p->year, TsChange);       
        //verbPrintf(2, grid.com, "For year: %f, seaChange: %f\n", grid.p->year, seaChange);       
 
        if (seaChange != 0 || TsChange != 0) {
          for (PetscInt i = grid.xs; i<grid.xs+grid.xm; ++i) {
            for (PetscInt j = grid.ys; j<grid.ys+grid.ym; ++j) {
              Ts[i][j] += TsChange;
              bed[i][j] = origBed[i][j] - seaChange;
            }
          }
        }
      }
      ierr = VecRestoreArray(vIceCoreDeltaSea, &deltaSea); CHKERRQ(ierr);
      ierr = VecRestoreArray(vIceCoreDeltaT, &deltaT); CHKERRQ(ierr);
      ierr = VecRestoreArray(vIceCoreTimeT, &iceCoreTimeT); CHKERRQ(ierr);
      ierr = VecRestoreArray(vIceCoreTimeSea, &iceCoreTimeSea); CHKERRQ(ierr);
      ierr = DAVecRestoreArray(grid.da2, vTs, &Ts); CHKERRQ(ierr);
      ierr = DAVecRestoreArray(grid.da2, vbed, &bed); CHKERRQ(ierr);
      ierr = DAVecRestoreArray(grid.da2, vOrigBed, &origBed); CHKERRQ(ierr);

      ierr = DALocalToLocalBegin(grid.da2, vbed, INSERT_VALUES, vbed); CHKERRQ(ierr);
      ierr = DALocalToLocalEnd(grid.da2, vbed, INSERT_VALUES, vbed); CHKERRQ(ierr);
      ierr = DALocalToLocalBegin(grid.da2, vTs, INSERT_VALUES, vTs); CHKERRQ(ierr);
      ierr = DALocalToLocalEnd(grid.da2, vTs, INSERT_VALUES, vTs); CHKERRQ(ierr);

      if (seaChange != 0) {
        updateSurfaceElevationAndMask();
      }
 
    }
    return 0;
}


PetscErrorCode IceGRNModel::calculateMeanAnnual(PetscScalar h, PetscScalar lat, PetscScalar *val) {
  //EISMINT surface temperature model
  PetscScalar Z = PetscMax(h, 20 * (lat - 65));
  *val = 49.13 - 0.007992 * Z - 0.7576 * (lat);
  return 0;
}


PetscScalar IceGRNModel::getSummerWarming(
       const PetscScalar elevation, const PetscScalar latitude, const PetscScalar Ta) const {
  // this is virtual in IceModel
  // EISMINT summer surface temperature model (expressed as warming above mean annual)
  // Ta, Ts in degrees C; Ta is mean annual, Ts is summer peak
  const PetscScalar Ts = 30.38 - 0.006277 * elevation - 0.3262 * latitude;
  return Ts - Ta;
}


PetscErrorCode IceGRNModel::updateTs() {
  PetscErrorCode ierr;
  PetscScalar val;
  PetscScalar **Ts, **lat, **h;
  
  ierr = DAVecGetArray(grid.da2, vTs, &Ts); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vLatitude, &lat); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vh, &h); CHKERRQ(ierr);
  for (PetscInt i = grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j<grid.ys+grid.ym; ++j) {
      calculateMeanAnnual(h[i][j], lat[i][j], &val);
      Ts[i][j] = val + ice.meltingTemp;
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vTs, &Ts); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vLatitude, &lat); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vh, &h); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceGRNModel::ellePiecewiseFunc(PetscScalar lon, PetscScalar *lat) {
  // first function
  float l1_x1 = -68.18, l1_y1 = 80.1;
  float l1_x2 = -62, l1_y2 = 82.24;
  float m, b;  // piecewise boundaries

  m = (l1_y1 - l1_y2) / (l1_x1 - l1_x2);
  b = (l1_y2) - m * (l1_x2);
  *lat = m * lon + b;
  return 0;
}


PetscErrorCode IceGRNModel::cleanExtraLand(){
  PetscErrorCode ierr;
  PetscScalar lat_line;
  // remove mask SE of the following point
  float ice_lon = 30, ice_lat = 67;
  PetscScalar **lat, **lon, **mask;

  ierr = DAVecGetArray(grid.da2, vLatitude, &lat); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vLongitude, &lon); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vMask, &mask); CHKERRQ(ierr);
  for (PetscInt i = grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j<grid.ys+grid.ym; ++j) {
      ellePiecewiseFunc(lon[i][j], &lat_line);
      if (lat[i][j]>lat_line) {
          mask[i][j] = MASK_FLOATING_OCEAN0;
      } else if (lat[i][j] < ice_lat && lon[i][j] > -ice_lon) {
        mask[i][j] = MASK_FLOATING_OCEAN0;
      }
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vLatitude, &lat); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vLongitude, &lon); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vMask, &mask); CHKERRQ(ierr);

  ierr = DALocalToLocalBegin(grid.da2, vMask, INSERT_VALUES, vMask); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2, vMask, INSERT_VALUES, vMask); CHKERRQ(ierr);

  return 0;
}

