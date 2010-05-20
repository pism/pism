// Copyright (C) 2010 Constantine Khroulev
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
#include "grid.hh"
#include "PISMIO.hh"

//! A class storing and writing PISM profiling event data.
class PISMEvent {
public:
  PISMEvent();
  string name,			//!< NetCDF variable name
    description;		//!< NetCDF variable long_name attribute
  int parent;			//!< index of the parent event
  PetscLogDouble start_time,	//!< event start time
    start_cpu_time;		//!< event start CPU time
  double total_time,		//!< total time spent in an event; includes
				//!< time spent doing nothing
    total_cpu_time;		//!< total CPU time spent in an event
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
  PISMProf(IceGrid &grid, const NCConfigVariable &config);
  ~PISMProf();
  int create(string name, string description);
  void begin(int index);
  void end(int index);
  PetscErrorCode barrier();
  PetscErrorCode save_report(string filename);
protected:
  vector<PISMEvent> events;
  int current_event;
  IceGrid *g;
  PetscErrorCode create_profiling_grid(IceGrid &grid, const NCConfigVariable &config);
  PetscErrorCode save_report(int index, int ncid, int varid, int varid_cpu);
  PetscErrorCode find_variables(PISMIO &nc, string name, int &varid, int &varid_cpu);
  PetscErrorCode define_variable(int ncid, string name, int &varid);
  PetscErrorCode put_att_text(int ncid, int varid, string name, string text);
};

#endif // __PISMProf_hh
