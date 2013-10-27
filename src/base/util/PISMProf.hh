// Copyright (C) 2010, 2011, 2012, 2013 Constantine Khroulev
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

#ifndef __PISMProf_hh
#define __PISMProf_hh

#ifndef PISM_PROFILE
#define PISM_PROFILE
#endif

#include <string>
#include <vector>
#include <petsc.h>
#include "PISMNCFile.hh"

//! \brief A class storing and writing PISM profiling event data.
class PISMEvent {
public:
  PISMEvent();
  std::string name,			//!< NetCDF variable name
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
  int create(std::string name, std::string description);
  int get(std::string name);
  void begin(int index);
  void end(int index);
  PetscErrorCode barrier();
  PetscErrorCode save_report(std::string filename);
  void set_grid_size(int n);
  int Nx, Ny;
protected:
  std::vector<PISMEvent> events;
  int current_event;
  PetscMPIInt rank, size;
  MPI_Comm com;

  PetscErrorCode save_report(int index, const PISMNCFile &nc, std::string name);
  PetscErrorCode find_variables(PISMNCFile &nc, std::string name, bool &exists);
  PetscErrorCode define_variable(const PISMNCFile &nc, std::string name);
  PetscErrorCode create_dimensions(const PISMNCFile &nc);
};

#endif // __PISMProf_hh
