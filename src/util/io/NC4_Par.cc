// Copyright (C) 2012, 2013, 2014, 2015, 2016, 2017, 2019 PISM Authors
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

#include "NC4_Par.hh"
#include "pism/util/error_handling.hh"

// netcdf_par.h has to be included *after* mpi.h and after netcdf.h
//
// note that we don't need to define MPI_INCLUDED because this code is built *only* if we
// have a parallel NetCDF library.
extern "C" {
#include <netcdf.h>
#include <netcdf_par.h>
}

namespace pism {
namespace io {

//! \brief Prints an error message; for debugging.
static void check(const ErrorLocation &where, int return_code) {
  if (return_code != NC_NOERR) {
    throw RuntimeError(where, nc_strerror(return_code));
  }
}

void NC4_Par::open_impl(const std::string &fname,
                        IO_Mode mode,
                        int FileID,
                        const std::map<std::string, int> &dimsa) {
  MPI_Info info = MPI_INFO_NULL;
  int stat;

  int open_mode = mode == PISM_READONLY ? NC_NOWRITE : NC_WRITE;
  open_mode = open_mode | NC_MPIIO;

  stat = nc_open_par(fname.c_str(), open_mode, m_com, info, &m_file_id);

  check(PISM_ERROR_LOCATION, stat);
}

void NC4_Par::create_impl(const std::string &fname, int FileID) {
  MPI_Info info = MPI_INFO_NULL;
  int stat;

  stat = nc_create_par(fname.c_str(),
                       NC_NETCDF4 | NC_MPIIO,
                       m_com, info, &m_file_id);

  check(PISM_ERROR_LOCATION, stat);
}

void NC4_Par::set_access_mode(int varid, bool transposed) const {
  int stat;

  if (transposed) {
    // Use independent parallel access mode because it works. It would be
    // better to use collective mode, but I/O performance is ruined by
    // the transpose anyway.
    //
    // See https://bugtracking.unidata.ucar.edu/browse/NCF-152 for the description of the bug we're
    // avoiding here.
    stat = nc_var_par_access(m_file_id, varid, NC_INDEPENDENT); check(PISM_ERROR_LOCATION, stat);
  } else {
    // Use collective parallel access mode because it is faster (and because it
    // works in this case).
    stat = nc_var_par_access(m_file_id, varid, NC_COLLECTIVE); check(PISM_ERROR_LOCATION, stat);
  }
}



} // end of namespace io
} // end of namespace pism
