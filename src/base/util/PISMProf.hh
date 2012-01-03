// Copyright (C) 2010, 2011 Constantine Khroulev
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

#ifndef __PISMProf_hh
#define __PISMProf_hh

#ifndef PISM_PROFILE
#define PISM_PROFILE
#endif

#include <string>
#include <vector>
#include <petsc.h>
#include "NetCDF3Wrapper.hh"

/// @cond NAMESPACE_BROWSER
using namespace std;
/// @endcond

//! \brief A class storing and writing PISM profiling event data.
class PISMEvent {
public:
  PISMEvent();
  string name,			//!< NetCDF variable name
    description,		//!< NetCDF variable long_name attribute
    units;                      //!< NetCDF variable units
  int parent;			//!< index of the parent event
  PetscLogDouble start_time;	//!< event start time
  double total_time;		//!< total time spent in an event; includes
				//!< time spent doing nothing
  PetscLogEvent petsc_event;
};

//! PISM profiler class.
/*!
  Usage example:

  \code
  PISMProf *prof;
  prof = new PISMProf(grid, config);
  int event;

  event = prof->create("event_varname", "event_description");

  prof->begin(event);
  // do stuff
  prof->end(event);

  ierr = prof->save_report("prof.nc"); CHKERRQ(ierr); 

  delete prof;
  \endcode
 */
class PISMProf {
public:
  PISMProf(MPI_Comm c, PetscMPIInt r, PetscMPIInt s);
  ~PISMProf() {}
  int create(string name, string description);
  int get(string name);
  void begin(int index);
  void end(int index);
  PetscErrorCode barrier();
  PetscErrorCode save_report(string filename);
  void set_grid_size(int n);
  int Nx, Ny;
protected:
  vector<PISMEvent> events;
  int current_event;
  PetscMPIInt rank, size;
  MPI_Comm com;

  PetscErrorCode save_report(int index, const NetCDF3Wrapper &nc, int varid);
  PetscErrorCode find_variables(NetCDF3Wrapper &nc, string name, int &varid);
  PetscErrorCode define_variable(const NetCDF3Wrapper &nc, string name, int &varid);
  PetscErrorCode create_dimensions(const NetCDF3Wrapper &nc);
};

#endif // __PISMProf_hh
