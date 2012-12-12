// Copyright (C) 2012 PISM Authors
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

#include "PISMNC4_Par.hh"

// netcdf_par.h has to be included *after* mpi.h
extern "C" {
#include <netcdf_par.h>
}

int PISMNC4_Par::open(string fname, int mode) {
  MPI_Info info = MPI_INFO_NULL;
  int stat;

  m_filename = fname;

  stat = nc_open_par(m_filename.c_str(),
                     mode | NC_MPIIO,
                     com, info, &ncid);

  define_mode = false;

  return stat;
}

int PISMNC4_Par::create(string fname) {
  MPI_Info info = MPI_INFO_NULL;
  int stat;

  m_filename = fname;

  stat = nc_create_par(m_filename.c_str(),
                       NC_NETCDF4 | NC_MPIIO,
                       com, info, &ncid);
  define_mode = true;

  return stat;
}

int PISMNC4_Par::set_access_mode(int varid, bool mapped) const {
  int stat;

  if (mapped) {
    // Use independent parallel access mode because it works. It would be
    // better to use collective mode, but I/O performance is ruined by
    // "mapping" anyway.

    stat = nc_var_par_access(ncid, varid, NC_INDEPENDENT); check(stat);
  } else {
    // Use collective parallel access mode because it is faster (and because it
    // works in this case).
    stat = nc_var_par_access(ncid, varid, NC_COLLECTIVE); check(stat);
  }

  return stat;
}


