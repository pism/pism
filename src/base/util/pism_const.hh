// Copyright (C) 2007--2012 Jed Brown, Ed Bueler and Constantine Khroulev
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

#ifndef __pism_const_hh
#define __pism_const_hh

#include <gsl/gsl_math.h>
#include <petsc.h>
#include <string>
#include <vector>
#include <set>

#include "udunits.h"

// use namespace std BUT remove trivial namespace browser from doxygen-erated HTML source browser
/// @cond NAMESPACE_BROWSER
using namespace std;
/// @endcond

extern const char *PISM_Revision;
extern const char *PISM_DefaultConfigFile;

const PetscScalar secpera    = 3.15569259747e7; // The constant used in UDUNITS
						// (src/udunits/pismudunits.dat)
const PetscScalar pi         = M_PI;		// defined in gsl/gsl_math.h

enum PismMask {
  MASK_UNKNOWN          = -1,
  MASK_ICE_FREE_BEDROCK = 0,
  MASK_GROUNDED   = 2,
  MASK_FLOATING         = 3,
  MASK_ICE_FREE_OCEAN   = 4
};

enum PismIcebergMask {
  ICEBERGMASK_NO_ICEBERG = -3,
  ICEBERGMASK_NOT_SET = 0,
  ICEBERGMASK_ICEBERG_CAND = 2,
  ICEBERGMASK_STOP_OCEAN = 3,
  ICEBERGMASK_STOP_ATTACHED = 4
};

const PetscInt TEMPORARY_STRING_LENGTH = 32768; // 32KiB ought to be enough.

bool is_increasing(const vector<double> &a);

PetscErrorCode setVerbosityLevel(PetscInt level);
PetscInt       getVerbosityLevel();
PetscErrorCode verbPrintf(const int thresh, MPI_Comm comm,const char format[],...);

void endPrintRank();

//void PISMEnd()  __attribute__((noreturn));
//void PISMEndQuiet()  __attribute__((noreturn));
void PISMEnd();
void PISMEndQuiet();

string pism_timestamp();
string pism_username_prefix(MPI_Comm com);
string pism_args_string();
string pism_filename_add_suffix(string filename, string separator, string suffix);

bool ends_with(string str, string suffix);

inline bool set_contains(set<string> S, string name) {
  return (S.find(name) != S.end());
}

//! \brief Convert a quantity from unit1 to unit2.
/*!
 * Example: convert(1, "m/year", "m/s").
 *
 * Please avoid using in computationally-intensive code.
 */
inline double convert(double value, const char spec1[], const char spec2[]) {
  utUnit unit1, unit2;
  double slope, intercept;
  int errcode;

  errcode = utScan(spec1, &unit1);
  if (errcode != 0) {
#if (PISM_DEBUG==1)
    PetscPrintf(MPI_COMM_SELF, "utScan failed trying to parse %s\n", spec1);
    PISMEnd();
#endif
    return GSL_NAN;
  }

  errcode = utScan(spec2, &unit2);
  if (errcode != 0) {
#if (PISM_DEBUG==1)
    PetscPrintf(MPI_COMM_SELF, "utScan failed trying to parse %s\n", spec2);
    PISMEnd();
#endif
    return GSL_NAN;
  }

  errcode = utConvert(&unit1, &unit2, &slope, &intercept);
  if (errcode != 0) {
#if (PISM_DEBUG==1)
    PetscPrintf(MPI_COMM_SELF, "utConvert failed trying to convert %s to %s\n", spec1, spec2);
    PISMEnd();
#endif
    return GSL_NAN;
  }

  return value * slope + intercept;
}

inline PetscErrorCode PISMGlobalMin(PetscReal *local, PetscReal *result, MPI_Comm comm)
{
  return MPI_Allreduce(local,result,1,MPIU_REAL,MPI_MIN,comm);
}

inline PetscErrorCode PISMGlobalMax(PetscReal *local, PetscReal *result, MPI_Comm comm)
{
  return MPI_Allreduce(local,result,1,MPIU_REAL,MPI_MAX,comm);
}

inline PetscErrorCode PISMGlobalSum(PetscReal *local, PetscReal *result, MPI_Comm comm)
{
  return MPI_Allreduce(local,result,1,MPIU_REAL,MPI_SUM,comm);
}

#endif
