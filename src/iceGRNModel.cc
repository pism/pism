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

#include <cstring>
#include <netcdf.h>
#include "nc_util.hh"
#include "grid.hh"
#include "materials.hh"
#include "iceModel.hh"
#include "iceGRNModel.hh"

const int STD_DEV = 5;
const PetscScalar POS_DEGREE_DAY_FACTOR_SNOW = .003; // m d^-1 C^-1
const PetscScalar POS_DEGREE_DAY_FACTOR_ICE = .008;  // m d^-1 C^-1
const PetscScalar REFREEZE_FACTOR = .6;
const PetscScalar DEFAULT_ABLATION_IN_OCEAN0 = 20.0;   // m/a
const PetscScalar   G_geothermal   = 0.050;      // J/m^2 s; geo. heat flux
const int CONST_PIECE_FWD_INTERP = 0;
const int CONST_PIECE_BCK_INTERP = 1;
const int LINEAR_INTERP = 2;

PetscErrorCode nc_check(int stat) {
  if (stat) {
    SETERRQ1(1, "NC_ERR: %s\n", nc_strerror(stat));
  }
  return 0;
}

int IceGRNModel::getTestNum() {
  return testnum;
}

PetscErrorCode IceGRNModel::copySnowAccum() {
  PetscErrorCode ierr;
  ierr = VecCopy(vSnowAccum, vAccum); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode IceGRNModel::copyOrigBed() {
  PetscErrorCode ierr;
  ierr = VecCopy(vOrigBed, vbed); CHKERRQ(ierr);
  return 0;
}

IceGRNModel::IceGRNModel(IceGrid &g, IceType &i) : IceModel(g, i) {
  // only call parent's constructor
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

PetscErrorCode IceGRNModel::setFromOptions() {
  PetscErrorCode ierr;
  PetscTruth usrTestSet;
  int usr_test;

  ierr = PetscOptionsGetInt(PETSC_NULL, "-test", &usr_test,
                            &usrTestSet); CHKERRQ(ierr);
  if (usrTestSet) {
    testnum = usr_test;
  } else {
    testnum=1;
  }

  enhancementFactor = 3;

  ierr = IceModel::setFromOptions(); CHKERRQ(ierr);  // note -e will override e=3 if set
  return 0;
}


PetscErrorCode IceGRNModel::initFromOptions() {
  PetscErrorCode ierr;
  int stat;

  ierr = IceModel::initFromOptions(); CHKERRQ(ierr);
  ierr = saveOrigVecs(); CHKERRQ(ierr);

  char inFile[PETSC_MAX_PATH_LEN], coreFile[PETSC_MAX_PATH_LEN];
  int usr_year, usr_seed;
  PetscTruth inFileSet, bootFileSet, usrYearSet, usrSeedSet,
             forceSet;
  
  verbPrintf(3, grid.com,
             "Running Greenland Eismint test, so ghf is being set to: %f\n",
             G_geothermal);
  ierr = VecSet(vGhf, G_geothermal); CHKERRQ(ierr);
  

  ierr = PetscOptionsGetString(PETSC_NULL, "-if", inFile,
                               PETSC_MAX_PATH_LEN, &inFileSet); CHKERRQ(ierr);

  ierr = PetscOptionsGetString(PETSC_NULL, "-bif", inFile,
                               PETSC_MAX_PATH_LEN, &bootFileSet); CHKERRQ(ierr);

  ierr = PetscOptionsGetString(PETSC_NULL, "-forcing", coreFile,
                               PETSC_MAX_PATH_LEN, &forceSet); CHKERRQ(ierr);

  ierr = PetscOptionsGetInt(PETSC_NULL, "-ys", &usr_year,
                            &usrYearSet); CHKERRQ(ierr);

  ierr = PetscOptionsGetInt(PETSC_NULL, "-repeatable", &usr_seed,
                            &usrSeedSet); CHKERRQ(ierr);

  if (inFileSet == PETSC_TRUE) {
    if (bootFileSet) {
      verbPrintf(1, grid.com, "WARNING: -bif and -if given. using -if.\n");
    }
  } else if (bootFileSet == PETSC_TRUE) {
    // after we set the new temperatures, we need
    // to set the 3D temps again
    ierr = fillTs(); CHKERRQ(ierr);
    ierr = putTempAtDepth(); CHKERRQ(ierr);
    ierr = cleanExtraLand(); CHKERRQ(ierr);
  } else {
    SETERRQ(1, "Error: need input file.\n");
  }

  // If the user is using the ice core data, it is necessary
  // to set the current year to be in the past so that we can
  // use the ice core data to model the ice sheet until the
  // present
  if (usrYearSet) {
    grid.p->year = usr_year;
  }

  if (forceSet == PETSC_TRUE) {
    int ncid, v_dT, v_dSea, v_time;
    PetscScalar *years;

    iceCoreIdx = 0;
    // read in core data from file
    if (grid.rank == 0) {
      stat = nc_open(coreFile, 0, &ncid); CHKERRQ(nc_check(stat));
      stat = nc_inq_varid(ncid, "delta_T", &v_dT); CHKERRQ(nc_check(stat));
      stat = nc_inq_varid(ncid, "delta_Sea", &v_dSea); CHKERRQ(nc_check(stat));
      stat = nc_inq_varid(ncid, "t", &v_time); CHKERRQ(nc_check(stat));
      ierr = getInterpolationCode(ncid, v_dT, &gripDeltaTInterp); CHKERRQ(ierr);
      ierr = getInterpolationCode(ncid, v_dSea, &gripDeltaSeaInterp); CHKERRQ(ierr);
    }
    MPI_Bcast(&gripDeltaTInterp, 1, MPI_INT, 0, grid.com);
    MPI_Bcast(&gripDeltaSeaInterp, 1, MPI_INT, 0, grid.com);
    ierr = ncVarBcastVec(ncid, v_dT, &vIceCoreDeltaT); CHKERRQ(ierr);
    ierr = ncVarBcastVec(ncid, v_dSea, &vIceCoreDeltaSea); CHKERRQ(ierr);
    ierr = ncVarBcastVec(ncid, v_time, &vIceCoreTime); CHKERRQ(ierr);
   
    ierr = VecGetArray(vIceCoreTime, &years); CHKERRQ(ierr);
    int r, l=0;
    ierr = VecGetLocalSize(vIceCoreTime, &iceCoreLen); CHKERRQ(ierr);
    r = iceCoreLen;
    // do a binary search to find where our year fits in.
    while (r > l + 1) {
      PetscInt j = (r + l)/2;
      if(grid.p->year < years[j]) {
        r = j;
      } else {
        l = j;
      }
    }    
    iceCoreIdx = l;
    // maybe we are already past our place.
    if (l >= iceCoreLen) {
      climateForcing = PETSC_FALSE;
      verbPrintf(1, grid.com, "No climate forcing data. Using last DeltaTs and DeltaSea.\n");
    } else {
      climateForcing = PETSC_TRUE;
    }

  } else if (testnum == 2) {
    SETERRQ(1, "Error: need ice core data for this test.\n");
  }

  // initialize the random number generator
  #if (WITH_GSL)
  rand_gen = gsl_rng_alloc(gsl_rng_taus);
  gsl_rng_set(rand_gen, usrSeedSet ? usr_seed : (int)grid.p->year);
  #else
  verbPrintf(1, grid.com, "WARNING: No randomness in weather because GSL is not included.\n");
  #endif
                        
  if (!isInitialized()) {
    SETERRQ(1, "Model has not been initialized.\n");
  }       

  return 0;
}

PetscErrorCode IceGRNModel::createVecs() {
  PetscErrorCode ierr;

  ierr = IceModel::createVecs(); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode IceGRNModel::destroyVecs() {
  PetscErrorCode ierr;

  ierr = IceModel::destroyVecs(); CHKERRQ(ierr);

  ierr = VecDestroy(vIceCoreDeltaT); CHKERRQ(ierr);
  ierr = VecDestroy(vIceCoreDeltaSea); CHKERRQ(ierr);
  ierr = VecDestroy(vIceCoreTime); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode IceGRNModel::saveOrigVecs() {
  PetscErrorCode ierr;

  // the data that represents the snow fall is actually
  // in vAccum, so we need to make an equivalent vector
  // and copy over the data.
  ierr = VecDuplicate(vAccum, &vSnowAccum); CHKERRQ(ierr);
  ierr = VecCopy(vAccum, vSnowAccum); CHKERRQ(ierr);

  // we need to save the original bed information and copy
  // it back before we write it to the nc file.
  ierr = VecDuplicate(vbed, &vOrigBed); CHKERRQ(ierr);
  ierr = VecCopy(vbed, vOrigBed); CHKERRQ(ierr);
  

  ierr = calculateNetAccum(); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode IceGRNModel::additionalAtStartTimestep() {
    PetscErrorCode ierr;
    PetscTruth integralYear = PETSC_FALSE;
    PetscScalar *deltaT, *iceCoreTime, *deltaSea;
    PetscScalar **Ts, **bed, **origBed;
    PetscScalar TsChange, seaChange;

    // at the beginning of each time step
    // we need to redo the model of the surface
    // temperatures.
    ierr = fillTs(); CHKERRQ(ierr);
    if (testnum==2 && climateForcing == PETSC_TRUE) {

      // apply the climate forcing
      ierr = VecGetArray(vIceCoreDeltaT, &deltaT); CHKERRQ(ierr);
      ierr = VecGetArray(vIceCoreDeltaSea, &deltaSea); CHKERRQ(ierr);
      ierr = VecGetArray(vIceCoreTime, &iceCoreTime); CHKERRQ(ierr);
      ierr = DAVecGetArray(grid.da2, vTs, &Ts); CHKERRQ(ierr);
      ierr = DAVecGetArray(grid.da2, vbed, &bed); CHKERRQ(ierr);
      ierr = DAVecGetArray(grid.da2, vOrigBed, &origBed); CHKERRQ(ierr);

      // if there was a large time step, it is possible
      // that we skip over multiple entries
      while (iceCoreIdx < iceCoreLen && grid.p->year > iceCoreTime[iceCoreIdx]) {
         iceCoreIdx++;
      }
      if (iceCoreIdx >= iceCoreLen) {
        verbPrintf(1, grid.com, "No more climate data. Turning off Climate Forcing.\n");
        climateForcing = PETSC_FALSE;
      } else {

        // if we have exact data, use it
        if (grid.p->year == iceCoreTime[iceCoreIdx]) {
          TsChange = deltaT[iceCoreIdx];
          seaChange = deltaSea[iceCoreIdx];
        } else { // otherwise, we need to interpolate.
          PetscScalar y0, y1;
          switch (gripDeltaTInterp) {
            case CONST_PIECE_BCK_INTERP:
              // use the data point we are infront of
              if (iceCoreIdx == 0) {
                TsChange = 0;
                seaChange = 0;
                verbPrintf(2, grid.com, "Year before Climate Data, using TsChange 0\n");
              } else {
                TsChange = deltaT[iceCoreIdx-1];
                seaChange = deltaSea[iceCoreIdx-1];
              }
              break;
            case CONST_PIECE_FWD_INTERP:
              TsChange = deltaT[iceCoreIdx];
              seaChange = deltaSea[iceCoreIdx];
              break;
            case LINEAR_INTERP:
              if (iceCoreIdx == 0) {
                TsChange = 0;
                seaChange = 0;
                verbPrintf(2, grid.com, "Year before Climate Data, using TsChange 0\n");
              } else {
                y0 = deltaT[iceCoreIdx-1];
                y1 = deltaT[iceCoreIdx];
                TsChange = y0+(y1-y0)/(iceCoreTime[iceCoreIdx]-iceCoreTime[iceCoreIdx-1])*(grid.p->year-iceCoreTime[iceCoreIdx-1]);

                y0 = deltaSea[iceCoreIdx-1];
                y1 = deltaSea[iceCoreIdx];
                seaChange = y0+(y1-y0)/(iceCoreTime[iceCoreIdx]-iceCoreTime[iceCoreIdx-1])*(grid.p->year-iceCoreTime[iceCoreIdx-1]);
              }
              break;
            default:
              SETERRQ(1, "Unknown Interpolation method");
          }
        }

        //verbPrintf(2, grid.com, "For year: %f, TsChange: %f\n", grid.p->year, TsChange);       
        //verbPrintf(2, grid.com, "For year: %f, seaChange: %f\n", grid.p->year, seaChange);       
 
        for (PetscInt i = grid.xs; i<grid.xs+grid.xm; ++i) {
          for (PetscInt j = grid.ys; j<grid.ys+grid.ym; ++j) {
            Ts[i][j] += TsChange;
            bed[i][j] = origBed[i][j] - seaChange;
          }
        }
      }
      ierr = VecRestoreArray(vIceCoreDeltaT, &deltaT); CHKERRQ(ierr);
      ierr = VecRestoreArray(vIceCoreDeltaSea, &deltaSea); CHKERRQ(ierr);
      ierr = VecRestoreArray(vIceCoreTime, &iceCoreTime); CHKERRQ(ierr);
      ierr = DAVecRestoreArray(grid.da2, vTs, &Ts); CHKERRQ(ierr);
      ierr = DAVecRestoreArray(grid.da2, vbed, &bed); CHKERRQ(ierr);
      ierr = DAVecRestoreArray(grid.da2, vOrigBed, &origBed); CHKERRQ(ierr);

      ierr = DALocalToLocalEnd(grid.da2, vTs, INSERT_VALUES, vTs); CHKERRQ(ierr);
      ierr = DALocalToLocalEnd(grid.da2, vbed, INSERT_VALUES, vbed); CHKERRQ(ierr);
      
    }
    if ((fmod(grid.p->year, 1.0) + 5e-6) < 1e-5 ) {
      integralYear = PETSC_TRUE;
    }

    // the net accumulation depends on Ta, Ts, lat,
    // d18O, and snow accumulation
    if (integralYear) {
      calculateNetAccum();
    }

    if (grid.p->year<0) {
      maxdt_temporary = fmod(fabs(grid.p->year), 1.0) * secpera;
    } else {
      maxdt_temporary = (1.0 - fmod(fabs(grid.p->year), 1.0)) * secpera;
    }
    
    return 0;
}


PetscErrorCode IceGRNModel::calculateMeanAnnual(PetscScalar h, PetscScalar lat, PetscScalar *val) {
  PetscScalar Z;

  Z = PetscMax(h, 20 * (lat - 65));
  *val = 49.13 - 0.007992 * Z - 0.7576 * (lat);
  return 0;
}

PetscErrorCode IceGRNModel::calculateSummerTemp(PetscScalar h, PetscScalar lat, PetscScalar *val) {
  *val = 30.38 - .006277 * h - .3262 * lat;
  return 0;
}

PetscErrorCode IceGRNModel::fillTs() {
  PetscErrorCode ierr;
  PetscScalar val;
  PetscScalar **Ts, **lat, **h, **H, **b;
  
  ierr = DAVecGetArray(grid.da2, vTs, &Ts); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vLatitude, &lat); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vh, &h); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vbed, &b); CHKERRQ(ierr);
  for (PetscInt i = grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j<grid.ys+grid.ym; ++j) {
      //EISMINT surface temperature model
      calculateMeanAnnual(h[i][j], lat[i][j], &val);
      Ts[i][j] = val + ice.meltingTemp;
    }
  }
  
  ierr = DAVecRestoreArray(grid.da2, vTs, &Ts); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vLatitude, &lat); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vh, &h); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vbed, &b); CHKERRQ(ierr);

  ierr = DALocalToLocalBegin(grid.da2, vTs, INSERT_VALUES, vTs); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2, vTs, INSERT_VALUES, vTs); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode IceGRNModel::calculateTempFromCosine(PetscScalar Ts, PetscScalar Ta,
                                                    int t, PetscScalar *temp) {
  float summer_offset = 243; // approximately August 1st
  // NOTE: Ta means "Mean annual Temperature" and
  //       Ts means "Summer Temperature"
  // These values are computed according to the
  // functions given in the EISMINT Greenland
  // experiment. also, they can be computed using
  // the follow functions:
  // Ta: calculateMeanAnnual
  // Ts: calculateSummerTemp
  *temp = (Ts-Ta) * cos((2 * PETSC_PI * (t - summer_offset)) / 365) + Ta;
  return 0;
}

PetscErrorCode IceGRNModel::calculateNetAccum() {
  PetscErrorCode ierr;
  PetscScalar **accum, **surface_t, **h, **lat, **snow_accum;
  PetscScalar val, pos_sum;
  PetscScalar mean_annual, summer;
  float rand_values[365];
  PetscScalar snow_melt, ice_melt, snow_year;

  // calculate a random number for each day
  // in the year. The same numbers will be
  // used for all points on the grid
  #if (WITH_GSL)
  for (int x=0; x<365; x++) {
    rand_values[x] = gsl_ran_gaussian(rand_gen, STD_DEV);
  }
  #endif

  ierr = DAVecGetArray(grid.da2, vLatitude, &lat); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vh, &h); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vTs, &surface_t); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vAccum, &accum); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vSnowAccum, &snow_accum); CHKERRQ(ierr);
  
  for (PetscInt i = grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j<grid.ys+grid.ym; ++j) {
      if (snow_accum[i][j]*secpera == (float)- DEFAULT_ABLATION_IN_OCEAN0) {
        continue;
      }
      pos_sum = 0;
      ierr = calculateMeanAnnual(h[i][j], lat[i][j], &mean_annual); CHKERRQ(ierr);
      ierr = calculateSummerTemp(h[i][j], lat[i][j], &summer); CHKERRQ(ierr);
      for (int x=0; x<365; x++){
        ierr = calculateTempFromCosine(summer, mean_annual, x, &val); CHKERRQ(ierr);
        val += rand_values[x];
        if (val > 0) {
          pos_sum += val;
        }
      }
      snow_year = snow_accum[i][j]*secpera;
      snow_melt = pos_sum * POS_DEGREE_DAY_FACTOR_SNOW;
      if (snow_melt!=0){
      }
      if (snow_melt <= snow_year) {
        accum[i][j] = ((snow_year-snow_melt) + (snow_melt * REFREEZE_FACTOR))/secpera;
      } else {
        accum[i][j] = (snow_accum[i][j] * REFREEZE_FACTOR);
        ice_melt = (snow_melt-snow_year)/POS_DEGREE_DAY_FACTOR_SNOW*POS_DEGREE_DAY_FACTOR_ICE;
        accum[i][j] -= ice_melt/secpera;
      }
    }
  }

  ierr = DAVecRestoreArray(grid.da2, vLatitude, &lat); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vh, &h); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vTs, &surface_t); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vAccum, &accum); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vSnowAccum, &snow_accum); CHKERRQ(ierr);

  ierr = DALocalToLocalBegin(grid.da2, vAccum, INSERT_VALUES, vAccum); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2, vAccum, INSERT_VALUES, vAccum); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode IceGRNModel::ellePiecewiseFunc(PetscScalar lon, PetscScalar *lat) {
  // first function
  float l1_x1 = -68.18, l1_y1 = 80.1;
  float l1_x2 = -62, l1_y2 = 82.24;

  // piecewise boundaries
  float m, b;

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

