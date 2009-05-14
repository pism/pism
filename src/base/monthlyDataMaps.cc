// Copyright (C) 2009 Ed Bueler
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


#include <petscda.h>
#include "pism_const.hh"
#include "grid.hh"
#include "iceModelVec.hh"
#include "monthlyDataMaps.hh"


#define MDMDEBUG 1
PetscErrorCode printIfDebug(bool debug, const char *message) {
  PetscErrorCode ierr;
  if (debug) {  ierr = verbPrintf(1,PETSC_COMM_WORLD,message); CHKERRQ(ierr);  }
  return 0;
}


MonthlyDataMaps::MonthlyDataMaps()  {
  grid = NULL;
}


MonthlyDataMaps::~MonthlyDataMaps() {
  for (PetscInt j = 0; j < 12; ++j) {
    vdata[j].destroy();  // destroys if allocated (created)
  }
}


PetscErrorCode MonthlyDataMaps::initFromFile(
            const char* prefix, const char* standard_name, 
            const char* filename, IceGrid* g) {
  PetscErrorCode ierr;
  printIfDebug(MDMDEBUG,"entering MonthlyDataMaps::initFromFile()\n");

  grid = g;

  // check file is really there; if so, read the dimensions of computational
  // grid so that we can set up a LocalInterpCtx for actual reading of climate data
  NCTool nc(grid);
  grid_info gi;
  ierr = nc.open_for_reading(filename); CHKERRQ(ierr);
  ierr = nc.get_grid_info_2d(gi); CHKERRQ(ierr);
  ierr = nc.close(); CHKERRQ(ierr);
  printIfDebug(MDMDEBUG,"  MonthlyDataMaps::initFromFile(): found file and read grid\n");
  
  LocalInterpCtx lic(gi, NULL, NULL, *grid); // 2D only; destructed by end of scope
  printIfDebug(MDMDEBUG,"  MonthlyDataMaps::initFromFile(): created lic\n");

  // for each month, create an IceModelVec2 and assign attributes
  for (PetscInt j = 0; j < 12; ++j) {
    char monthlyName[20], mstring[100];
    snprintf(monthlyName, 20, "%s%d", prefix, j+1);  // = "foo1" if prefix=foo and j=0
    ierr = vdata[j].create(*grid, monthlyName, false); CHKERRQ(ierr); // global; no ghosts
    snprintf(mstring, 100, 
             "some data for month %d of {1,..,12}", j+1); // FIXME
    ierr = vdata[j].set_attrs(
               "", // FIXME: need way to read pism_intent from "foo1" and set for all
               mstring,
               "",  // FIXME: need way to read units from "foo1" and set for all
               standard_name);  CHKERRQ(ierr);
    ierr = verbPrintf(2, grid->com, 
       "  reading month %d data '%s' ...", j+1, monthlyName); CHKERRQ(ierr); 
    ierr = vdata[j].regrid(filename, lic, true); CHKERRQ(ierr); // it *is* critical
  }
  ierr = verbPrintf(2, grid->com, "\n"); CHKERRQ(ierr); 
  
  printIfDebug(MDMDEBUG,"leaving MonthlyDataMaps::initFromFile()\n");
  return 0;
}


/*!
Assumes file is prepared in the sense that it exists and that global attributes are
already written.  See IceModel::dumpToFile() for how main PISM output file is
prepared.  Calls here do handle opening and closing the file.  We write in FLOAT 
not DOUBLE because these are expected to be imprecise at that level and not be
essential for restart accuracy.
 */
PetscErrorCode MonthlyDataMaps::write(const char *filename) {
  PetscErrorCode ierr;
  
  for (PetscInt j = 0; j < 12; ++j) {
    if (vdata[j].was_created()) {
      ierr = vdata[j].write(filename, NC_FLOAT); CHKERRQ(ierr);
    } else {
      ierr = verbPrintf(2, grid->com, 
        "MonthlyDataMaps ERROR:  vdata[%d].was_created() returned FALSE; ending ...\n", 
        j); CHKERRQ(ierr); 
      PetscEnd();
    }
  }
  return 0;
}


PetscErrorCode MonthlyDataMaps::getIndicesFromMonth(
              PetscScalar month, 
              PetscInt &currIndex, PetscInt &nextIndex, PetscScalar &lambda) {
  PetscErrorCode ierr;

  if ((month < 0.0) || (month >= 12.0)) {
    ierr = PetscPrintf(grid->com, 
       "MonthlyDataMaps ERROR:  invalid input month = %.5f in getIndicesFromMonth();\n"
       "      Should be in range 0 <= month < 12.0.  ENDING ...\n",month); CHKERRQ(ierr);
    PetscEnd();
  }

  lambda = month - floor(month);            // 0 <= lambda < 1
  currIndex = (int) floor(month);           // in {0,1,2,...,10,11}
  nextIndex = currIndex + 1;
  if (nextIndex == 12)   nextIndex = 0;     // in {0,1,2,...,10,11}
  return 0;
}


PetscErrorCode MonthlyDataMaps::getIndicesFromDay(
              PetscScalar myday, 
              PetscInt &currIndex, PetscInt &nextIndex, PetscScalar &lambda) {
  PetscErrorCode ierr;
  PetscScalar dummy;
  const PetscScalar
     sperd = 8.64e4,               // exact; from UDUNITS
     daypera = secpera / sperd;    // = 365.24...; FIXME: replace by config.get() call?
  PetscScalar // NOTE: modf returns 0.56 if given 23.56 but returns -0.56 if given -23.56
     month = 12.0 * modf(myday / daypera, &dummy); // -12 < month < 12
  if (month < 0)  month = month + 12.0;     // 0 <= month < 12
  ierr = getIndicesFromMonth(month, currIndex, nextIndex, lambda); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode MonthlyDataMaps::getIndicesFromTime(
              PetscScalar mytime, 
              PetscInt &currIndex, PetscInt &nextIndex, PetscScalar &lambda) {
  PetscErrorCode ierr;
  PetscScalar dummy;
  PetscScalar // NOTE: modf returns 0.56 if given 23.56 but returns -0.56 if given -23.56
     month = 12.0 * modf(mytime / secpera, &dummy); // -12 < month < 12
  if (month < 0)  month = month + 12.0;     // 0 <= month < 12
  ierr = getIndicesFromMonth(month, currIndex, nextIndex, lambda); CHKERRQ(ierr);
  return 0;
}


PetscScalar MonthlyDataMaps::interpolateMonthlyData(
              PetscInt i, PetscInt j,
              PetscScalar **currData, PetscScalar **nextData, PetscScalar lambda) {
  return  (1.0 - lambda) * currData[i][j] + lambda * nextData[i][j];
}

