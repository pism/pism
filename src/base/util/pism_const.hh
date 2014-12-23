// Copyright (C) 2007--2014 Jed Brown, Ed Bueler and Constantine Khroulev
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
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

#include <petsc.h>
#include <string>
#include <vector>
#include <set>
#include <map>

#ifdef PISM_USE_TR1
#include <tr1/memory>
#else
#include <memory>
#endif

#include "error_handling.hh"

namespace pism {

extern const char *PISM_Revision;
extern const char *PISM_DefaultConfigFile;

const int TEMPORARY_STRING_LENGTH = 32768; // 32KiB ought to be enough.

bool is_increasing(const std::vector<double> &a);

PetscErrorCode setVerbosityLevel(int level);
int       getVerbosityLevel();
void verbPrintf(const int thresh, MPI_Comm comm,const char format[],...);

std::string pism_timestamp();
std::string pism_username_prefix(MPI_Comm com);
std::string pism_args_string();
std::string pism_filename_add_suffix(std::string filename, std::string separator, std::string suffix);

PetscErrorCode GetTime(PetscLogDouble *result);

bool ends_with(std::string str, std::string suffix);

inline bool set_contains(std::set<std::string> S, std::string name) {
  return (S.find(name) != S.end());
}

inline void GlobalReduce(MPI_Comm comm, double *local, double *result, int count, MPI_Op op) {
  int err = MPI_Allreduce(local, result, count, MPIU_REAL, op, comm);
  PISM_CHK(err, 0, "MPI_Allreduce");
}

inline void GlobalMin(MPI_Comm comm, double *local, double *result, int count) {
  GlobalReduce(comm, local, result, count, MPI_MIN);
}

inline void GlobalMax(MPI_Comm comm, double *local, double *result, int count) {
  GlobalReduce(comm, local, result, count, MPI_MAX);
}

inline void GlobalSum(MPI_Comm comm, double *local, double *result, int count) {
  GlobalReduce(comm, local, result, count, MPI_SUM);
}

inline double GlobalMin(MPI_Comm comm, double local) {
  double result;
  GlobalMin(comm, &local, &result, 1);
  return result;
}

inline double GlobalMax(MPI_Comm comm, double local) {
  double result;
  GlobalMax(comm, &local, &result, 1);
  return result;
}

inline double GlobalSum(MPI_Comm comm, double local) {
  double result;
  GlobalSum(comm, &local, &result, 1);
  return result;
}

class Profiling {
public:
  Profiling();
  void begin(const char *name) const;
  void end(const char *name) const;
  void stage_begin(const char *name) const;
  void stage_end(const char *name) const;
private:
  PetscClassId m_classid;
  mutable std::map<std::string, PetscLogEvent> m_events;
  mutable std::map<std::string, PetscLogStage> m_stages;
};

} // end of namespace pism

#endif
