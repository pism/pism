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

#include "PISMProf.hh"
#include "nc_util.hh"

/// PISMEvent
PISMEvent::PISMEvent() {
  name = "unknown";
  start_time = -1;
  total_time = 0;
  parent = -1;
}

/// PISMProf
PISMProf::PISMProf(MPI_Comm c, PetscMPIInt r, PetscMPIInt s,
                   int my_Nx, int my_Ny) {
  com = c;
  rank = r;
  size = s;
  Nx = my_Nx;
  Ny = my_Ny;
  current_event = -1;
}

//! Create a profiling event.
/*!
 * Checks if an event with this name already exists.
 */
int PISMProf::create(string name, string description) {
  PISMEvent tmp;
  int index = get(name);
  
  if (index != -1)
    return index;

  tmp.name = name;
  tmp.description = description;

  events.push_back(tmp);

  return events.size() - 1;
}

//! \brief Get an integer (index) corresponding to an event.
/*!
 * Returns -1 if an event was not found.
 */
int PISMProf::get(string name) {
  
  for (int i = 0; i < events.size(); ++i)
    if (events[i].name == name)
      return i;
  
  return -1;
}

// This ensures that begin() and end() code is optimized away if
// PISM_PROFILE is not defined.
#ifndef PISM_PROFILE
void PISMProf::begin(int) {}
#else
//! Begin a profiling event.
void PISMProf::begin(int index) {
  PISMEvent &event = events[index];

  // fprintf(stderr, "Rank %d: Start '%s'\n", rank, events[index].description.c_str());

  event.parent = current_event;

  PetscGetTime(&event.start_time);

  current_event = index;
}
#endif

#ifndef PISM_PROFILE
void PISMProf::end(int) {}
#else
//! End a profiling event.
void PISMProf::end(int /*index*/) {
  PetscLogDouble time;
  PISMEvent &event = events[current_event];

  PetscGetTime(&time);

  event.total_time += (time - event.start_time);

  // fprintf(stderr, "Rank %d: End   '%s' (%3.9f s lapsed)\n", rank,
  // 	  events[current_event].description.c_str(), event.total_time);

  current_event = event.parent;
}
#endif

//! Save a profiling report to a file.
PetscErrorCode PISMProf::save_report(string filename) {
  PetscErrorCode ierr;
  int varid;
  NCTool nc(com, rank);

  ierr = nc.move_if_exists(filename.c_str()); CHKERRQ(ierr);

  ierr = nc.open_for_writing(filename.c_str()); CHKERRQ(ierr);

  ierr = create_dimensions(nc); CHKERRQ(ierr); 

  for (unsigned int j = 0; j < events.size(); ++j) {
    string &name = events[j].name;
    if (name == "unknown")
      continue;

    ierr = find_variables(nc, name, varid); CHKERRQ(ierr);

    ierr = save_report(j, nc, varid); CHKERRQ(ierr);

    int parent_index = events[j].parent;
    string parent;
    if (parent_index == -1)
      parent = "root";
    else
      parent = events[parent_index].name;
    string descr = events[j].description;

    ierr = put_att_text(nc, varid, "units", "seconds"); CHKERRQ(ierr); 
    ierr = put_att_text(nc, varid, "parent", parent); CHKERRQ(ierr); 
    ierr = put_att_text(nc, varid, "long_name", descr); CHKERRQ(ierr);
  }

  ierr = nc.close();
  return 0;
}

//! Save the report corresponding to an event events[index].
PetscErrorCode PISMProf::save_report(int index, const NCTool &nc, int varid) {
  PetscErrorCode ierr;
  const int data_tag = 1;
  MPI_Status mpi_stat;
  double data[1];
  size_t start[2];

  data[0] = events[index].total_time;

  ierr = nc.data_mode(); CHKERRQ(ierr);

  if (rank == 0) { // on processor 0: receive messages from every other
    for (int proc = 0; proc < size; proc++) {
      if (proc != 0) {
	MPI_Recv(data, 1, MPI_DOUBLE, proc, data_tag, com, &mpi_stat);
      }

      start[0] = proc % (Ny);
      start[1] = proc / (Ny);

      ierr = nc_put_var1_double(nc.get_ncid(), varid, start, &data[0]);
      CHKERRQ(check_err(ierr,__LINE__,__FILE__));

    }
  } else {  // all other processors: send data to processor 0
    MPI_Send(data, 1, MPI_DOUBLE, 0, data_tag, com);
  }

  return 0;
}

//! Find NetCDF variables to save total and CPU times into.
PetscErrorCode PISMProf::find_variables(NCTool &nc, string name, int &varid) {
  PetscErrorCode ierr;
  bool exists;

  ierr = nc.find_variable(name, &varid, exists); CHKERRQ(ierr); 
  if (!exists) {
    ierr = define_variable(nc, name, varid); CHKERRQ(ierr); 
  }

  return 0;
}

//! Define a NetCDF variable to store a profiling report in.
PetscErrorCode PISMProf::define_variable(const NCTool &nc, string name, int &varid) {
  PetscErrorCode ierr;
  int dimids[2];

  if (rank != 0) return 0;

  ierr = nc.define_mode(); CHKERRQ(ierr);

  ierr = nc_inq_dimid(nc.get_ncid(), "y", &dimids[0]); CHKERRQ(check_err(ierr,__LINE__,__FILE__));
  ierr = nc_inq_dimid(nc.get_ncid(), "x", &dimids[1]); CHKERRQ(check_err(ierr,__LINE__,__FILE__));
  
  ierr = nc_def_var(nc.get_ncid(), name.c_str(), NC_DOUBLE, 2, dimids, &varid);
  CHKERRQ(check_err(ierr,__LINE__,__FILE__));

  return 0;
}

//! Put a text attribute.
PetscErrorCode PISMProf::put_att_text(const NCTool &nc, int varid,
				      string name, string text) {
  PetscErrorCode ierr;

  if (rank != 0) return 0;

  ierr = nc.define_mode(); CHKERRQ(ierr); 

  ierr = nc_put_att_text(nc.get_ncid(), varid,
                         name.c_str(), text.size(), text.c_str());
  CHKERRQ(check_err(ierr,__LINE__,__FILE__));

  return 0;
}

//! Barrier synchroniation. (Calls MPI_Barrier.)
PetscErrorCode PISMProf::barrier() {

#ifdef PISM_PROFILE
  PetscErrorCode ierr;
  ierr = MPI_Barrier(com); CHKERRQ(ierr);
#endif

  return 0;
}

//! Creates x and y dimensions.
PetscErrorCode PISMProf::create_dimensions(const NCTool &nc) {

  int stat, dimid;
  if (rank == 0) {

    stat = nc.define_mode(); CHKERRQ(stat); 

    stat = nc_def_dim(nc.get_ncid(), "x", Nx, &dimid); CHKERRQ(check_err(stat,__LINE__,__FILE__));
    stat = nc_def_dim(nc.get_ncid(), "y", Ny, &dimid); CHKERRQ(check_err(stat,__LINE__,__FILE__));
  }

  return 0;
}
