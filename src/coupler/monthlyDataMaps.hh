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

#ifndef __monthlyDataMaps_hh
#define __monthlyDataMaps_hh

#include <petsc.h>
#include "../base/grid.hh"
#include "../base/LocalInterpCtx.hh"
#include "../base/iceModelVec.hh"


//! A class for reading/writing monthly data maps to/from a NetCDF file.
/*!
The data is monthly but is assumed (for now) to be periodic with period of one year.
For example, every January is the same.

Provides maps as an array of IceModelVec2.  Abstracts the linear interpolation
in time of these maps.  No spatial interpolation is intended, except at time of
reading the maps where the regrid() method of IceModelVec is used.

An example use in PISMAtmosphereCoupler is to read monthly snow temperatures
which are used in a PDD.  These monthly snow temperatures are for the layers
above many firn processes, so they are not directly meaningful as the upper
surface boundary condition for the conservation of energy equation in the ice.
They are typically the output of an atmospheric model.  They correspond to the
temperature which controls melting and refreezing processes within the snow.

Use this way:
\code
    IceGrid grid;
    PetscInt curr, next;
    PetscScalar lam, **currmap, **nextmap;
    PetscScalar day = 90;  // April Fool's day, approximately

    MonthlyDataMaps mdm;

    mdm.initFromFile("foo", "STD_NAME_IF_APPROPRIATE", "file.nc", &grid);
    mdm.getIndicesFromDay(day,curr,next,lam);
    mdm.vdata[curr].get_array(currmap);
    mdm.vdata[next].get_array(nextmap);
    for (PetscInt i = grid->xs; i<grid->xs+grid->xm; ++i) {
      for (PetscInt j = grid->ys; j<grid->ys+grid->ym; ++j) {
        result = mdm.interpolateMonthlyData(i,j,currmap,nextmap,lam);
      }
    }
    mdm.vdata[curr].end_access();
    mdm.vdata[next].end_access();
\endcode
 */
class MonthlyDataMaps {

public:
  MonthlyDataMaps();  //!< Sets pointers to NULL.  Does not allocate.
  virtual ~MonthlyDataMaps();  //!< Deallocates if allocated.

  /*! FIXME:  There should be a "month" dimension in the NetCDF?
      Needs a prefix, e.g. "foo".  Then data fields must have names 
      \c foo1, ...,\c foo12.  Units and other metadata will be read from file
      for foo1 and used for all months.  Used empty string "" for 
      mystandard_name if there isn't one.  */ 
  virtual PetscErrorCode initFromFile(
            const char* prefix, const char* standard_name, 
            const char* filename, IceGrid* g);

  /*! Writes data according to prefix, to a prepared file.  E.g. if prefix = "foo" 
      then written fields will have names \c foo1, ...,\c foo12.  Metadata will be
      the same for all fields.  */
  virtual PetscErrorCode write(const char *filename);

  /*! Input month must satisfy 0.0 <= month < 12.0, but can be any real value
      in that range.  The returned value lambda is used by interpolateMonthlyData(). */
  virtual PetscErrorCode getIndicesFromMonth(
              PetscScalar month, 
              PetscInt &currIndex, PetscInt &nextIndex, PetscScalar &lambda);

  /*! The day is treated mod 365.24.  In fact, the first action by this 
      procedure is to compute \code modf(day / 365.24, &dummy) \endcode
      to get fraction of year.  Calls getIndicesFromMonth().  The returned 
      value lambda is used by interpolateMonthlyData(). */
  virtual PetscErrorCode getIndicesFromDay(
              PetscScalar myday, 
              PetscInt &currIndex, PetscInt &nextIndex, PetscScalar &lambda);

  /*! The time in seconds is treated mod secpera (=3.155e7).
      Calls getIndicesFromMonth().  The returned value lambda is used by
      interpolateMonthlyData(). */
  virtual PetscErrorCode getIndicesFromTime(
              PetscScalar mytime, // seconds
              PetscInt &currIndex, PetscInt &nextIndex, PetscScalar &lambda);

  /*! Linear interpolation in time at grid point i,j.  Call getIndices..() 
      methods first to get PetscScalar** for current and next months data.  */
  virtual PetscScalar interpolateMonthlyData(// FIXME CAUTION: does not add TsOffset!!
              PetscInt i, PetscInt j,
              PetscScalar **currData, PetscScalar **nextData, PetscScalar lambda);

public:
  IceModelVec2 vdata[12];

protected:
  IceGrid*     grid;
};

#endif
