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
#include "PISMIO.hh"

/// PISMEvent
PISMEvent::PISMEvent() {
  name = "unknown";
  start_time = -1;
  start_cpu_time = -1;
  total_time = total_cpu_time = 0;
  parent = -1;
}

/// PISMProf

PISMProf::PISMProf(IceGrid &grid, const NCConfigVariable &config) {
  create_profiling_grid(grid, config);
  current_event = -1;
}

PISMProf::~PISMProf() {
  delete g;
}

//! Create a profiling event.
int PISMProf::create(string name, string description) {
  PISMEvent tmp;
  tmp.name = name;
  tmp.description = description;

  events.push_back(tmp);

  return events.size() - 1;
}

// This ensures that start_event and end_event code is optimized away if
// PISM_PROFILE is not defined.
#ifndef PISM_PROFILE
void PISMProf::begin(int) {}
#else
//! Begin a profiling event.
void PISMProf::begin(int index) {
  PISMEvent &event = events[index];

  // fprintf(stderr, "Rank %d: Start of event %d\n", g->rank, index);

  event.parent = current_event;

  PetscGetTime(&event.start_time);
  PetscGetCPUTime(&event.start_cpu_time);

  current_event = index;
}
#endif

#ifndef PISM_PROFILE
void PISMProf::end(int) {}
#else
//! End a profiling event.
void PISMProf::end(int index) {
  PetscLogDouble time, cpu_time;
  PISMEvent &event = events[current_event];

  PetscGetTime(&time);
  PetscGetCPUTime(&cpu_time);

  event.total_time += (time - event.start_time);
  event.total_cpu_time += (cpu_time - event.start_cpu_time);

  // fprintf(stderr, "Rank %d: End of event %d (%3.9f s lapsed)\n", g->rank, current_event, event.total_time);

  current_event = event.parent;
}
#endif

//! Create a grid that contains 2D information necessary to create x and y
//! dimensions (and variables) in the profiling output file.
PetscErrorCode PISMProf::create_profiling_grid(IceGrid &grid,
					       const NCConfigVariable &config) {
  g = new IceGrid(grid.com, grid.rank, grid.size, config);

  g->xs = grid.xs;
  g->ys = grid.ys;
  g->xm = grid.xm;
  g->ym = grid.ym;
  g->zlevels = new PetscScalar;
  g->zlevels[0] = 0;
  g->zblevels = new PetscScalar;
  g->zblevels[0] = 0;

  g->x0 = grid.x0;
  g->y0 = grid.y0;
  g->Lx = grid.Lx;
  g->Ly = grid.Ly;

  // Set Mx and My to Nx and Ny:
  g->Mx = grid.Nx;
  g->My = grid.Ny;

  g->Nx = grid.Nx;
  g->Ny = grid.Ny;

  g->dx = g->dy  = 0;
  g->Lz = g->Lbz = 0;
  g->Mz = g->Mbz = 1;

  return 0;
}

//! Save a profiling report to a file.
PetscErrorCode PISMProf::save_report(string filename) {
  PetscErrorCode ierr;
  int ncid, varid, varid_cpu;
  PISMIO nc(g);

  ierr = nc.open_for_writing(filename.c_str(),
			     false, // append
			     true); // check_dimensions
  CHKERRQ(ierr);

  ncid = nc.get_ncid();

  for (unsigned int j = 0; j < events.size(); ++j) {
    string &name = events[j].name;
    if (name == "unknown")
      continue;

    ierr = find_variables(nc, name, varid, varid_cpu); CHKERRQ(ierr);

    ierr = save_report(j, ncid, varid, varid_cpu); CHKERRQ(ierr);

    int parent_index = events[j].parent;
    string parent;
    if (parent_index == -1)
      parent = "root";
    else
      parent = events[parent_index].name;
    string descr = events[j].description,
      descr_cpu = descr + " (CPU time)",
      parent_cpu = parent + "_cpu";

    ierr = put_att_text(ncid, varid, "parent", parent); CHKERRQ(ierr); 
    ierr = put_att_text(ncid, varid, "long_name", descr); CHKERRQ(ierr);

    ierr = put_att_text(ncid, varid_cpu, "parent", parent_cpu); CHKERRQ(ierr); 
    ierr = put_att_text(ncid, varid_cpu, "long_name", descr_cpu); CHKERRQ(ierr);
  }

  ierr = nc.close();
  return 0;
}

//! Save the report corresponding to an event events[index].
PetscErrorCode PISMProf::save_report(int index, int ncid, int varid, int varid_cpu) {
  PetscErrorCode ierr;
  const int data_tag = 1;
  MPI_Status mpi_stat;
  double data[2];
  size_t start[2];

  data[0] = events[index].total_time;
  data[1] = events[index].total_cpu_time;

  if (g->rank == 0) { // on processor 0: receive messages from every other
    for (int proc = 0; proc < g->size; proc++) {
      if (proc != 0) {
	MPI_Recv(data, 2, MPI_DOUBLE, proc, data_tag, g->com, &mpi_stat);
      }

      // fprintf(stderr, "Event %d: proc %d: data: %f, %f\n", index, proc, data[0], data[1]);

      start[0] = proc % (g->Ny);
      start[1] = proc / (g->Ny);

      ierr = nc_put_var1_double(ncid, varid,     start, &data[0]);
      CHKERRQ(check_err(ierr,__LINE__,__FILE__));

      ierr = nc_put_var1_double(ncid, varid_cpu, start, &data[1]);
      CHKERRQ(check_err(ierr,__LINE__,__FILE__));

    }
  } else {  // all other processors: send data to processor 0
    MPI_Send(data, 2, MPI_DOUBLE, 0, data_tag, g->com);
  }

  return 0;
}

//! Find NetCDF variables to save total and CPU times into.
PetscErrorCode PISMProf::find_variables(PISMIO &nc, string name,
					int &varid, int &varid_cpu) {
  PetscErrorCode ierr;
  bool exists;
  int ncid;

  ncid = nc.get_ncid();

  ierr = nc.find_variable(name, &varid, exists); CHKERRQ(ierr); 
  if (!exists) {
    ierr = define_variable(ncid, name, varid); CHKERRQ(ierr); 
  }

  string name_cpu = name + "_cpu";
  ierr = nc.find_variable(name_cpu, &varid_cpu, exists); CHKERRQ(ierr);
  if (!exists) {
    ierr = define_variable(ncid, name_cpu, varid_cpu); CHKERRQ(ierr); 
  }

  return 0;
}

//! Define a NetCDF variable to store a profiling report in.
PetscErrorCode PISMProf::define_variable(int ncid, string name, int &varid) {
  PetscErrorCode ierr;
  int dimids[2];

  if (g->rank == 0) {
    ierr = nc_redef(ncid); CHKERRQ(check_err(ierr,__LINE__,__FILE__));

    ierr = nc_inq_dimid(ncid, "y", &dimids[0]); CHKERRQ(check_err(ierr,__LINE__,__FILE__));
    ierr = nc_inq_dimid(ncid, "x", &dimids[1]); CHKERRQ(check_err(ierr,__LINE__,__FILE__));
  
    ierr = nc_def_var(ncid, name.c_str(), NC_DOUBLE, 2, dimids, &varid);
    CHKERRQ(check_err(ierr,__LINE__,__FILE__));

    ierr = nc_enddef(ncid); CHKERRQ(check_err(ierr,__LINE__,__FILE__));
  }

  return 0;
}

//! Put a text attribute.
PetscErrorCode PISMProf::put_att_text(int ncid, int varid,
				      string name, string text) {
  PetscErrorCode ierr;

  if (g->rank != 0) return 0;

  ierr = nc_redef(ncid); CHKERRQ(check_err(ierr,__LINE__,__FILE__));
  {
    ierr = nc_put_att_text(ncid, varid,
			   name.c_str(), text.size(), text.c_str());
    CHKERRQ(check_err(ierr,__LINE__,__FILE__));
  }
  ierr = nc_enddef(ncid); CHKERRQ(check_err(ierr,__LINE__,__FILE__));

  return 0;
}

//! Barrier synchroniation. (Calls MPI_Barrier.)
PetscErrorCode PISMProf::barrier() {

#ifdef PISM_PROFILE
  PetscErrorCode ierr;
  ierr = MPI_Barrier(g->com); CHKERRQ(ierr);
#endif

  return 0;
}
