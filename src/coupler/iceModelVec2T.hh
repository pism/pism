// Copyright (C) 2009 Constantine Khroulev
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

#ifndef __IceModelVec2T_hh
#define __IceModelVec2T_hh

#include <petsc.h>
#include <string>
#include <vector>
#include "NCVariable.hh"
#include "iceModelVec.hh"
#include "nc_util.hh"

// use namespace std BUT remove trivial namespace browser from doxygen-erated HTML source browser
/// @cond NAMESPACE_BROWSER
using namespace std;
/// @endcond

//! A class for storing and accessing 2D time-series (for climate forcing)
/*!
  
 */
class IceModelVec2T : public IceModelVec {
public:
  IceModelVec2T();
  IceModelVec2T(const IceModelVec2T &other);
  virtual PetscErrorCode create(IceGrid &mygrid, const char my_short_name[], int N);
  virtual PetscErrorCode destroy();
  virtual PetscErrorCode get_arrays(PetscScalar** &a2, PetscScalar*** &a3);
  virtual PetscErrorCode init(string filename, string dim_name);
  virtual PetscErrorCode update(double t_years, double dt_years);
  virtual PetscErrorCode set_record(int n);
  virtual PetscErrorCode get_record(int n);
  virtual PetscErrorCode max_timestep(double t_years, double &dt_years);
  virtual PetscErrorCode interp(double t_years);			 // done, needs testing
  virtual PetscErrorCode interp(int i, int j, int N,
				PetscScalar *times, PetscScalar *results);
  virtual PetscErrorCode write(string filename, double t_years, nc_type nctype);
  virtual PetscErrorCode begin_access();
  virtual PetscErrorCode end_access();
private:
  NCTimeseries dimension;
  vector<double> times,		//!< all the times available in filename
    T;				//!< times stored in this IceModelVec2T
  string filename;
  LocalInterpCtx *lic;		//!< We store the context because we might need
				//!< to regrid many times
  DA da3;
  Vec v3;			//!< a 3D Vec used to store records
  void ***array3;
  int n_records;
  
  double j_min, j_max;

  virtual PetscErrorCode update(int start); // done, needs optimization
};


#endif // __IceModelVec2T_hh
