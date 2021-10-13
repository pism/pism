// Copyright (C) 2012, 2013, 2014, 2015, 2019, 2020 PISM Authors
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

#include "NC4_Serial.hh"

// The following is a stupid kludge necessary to make NetCDF 4.x work in
// serial mode in an MPI program:
#ifndef MPI_INCLUDED
#define MPI_INCLUDED 1
#endif
#include <netcdf.h>

#include "pism/util/error_handling.hh"
#include "pism/util/pism_utilities.hh"

namespace pism {
namespace io {

//! \brief Prints an error message; for debugging.
static void check(const ErrorLocation &where, int return_code) {
  if (return_code != NC_NOERR) {
    throw RuntimeError(where, nc_strerror(return_code));
  }
}

NC4_Serial::NC4_Serial(MPI_Comm c)
  : NC_Serial(c), m_compression_level(0) {
  // empty
}

void NC4_Serial::create_impl(const std::string &fname) {
  int stat = NC_NOERR;

  if (m_rank == 0) {
    stat = nc_create(fname.c_str(), NC_CLOBBER | NC_NETCDF4, &m_file_id);
  }

  MPI_Barrier(m_com);
  MPI_Bcast(&m_file_id, 1, MPI_INT, 0, m_com);
  MPI_Bcast(&stat, 1, MPI_INT, 0, m_com);

  check(PISM_ERROR_LOCATION, stat);
}

void NC4_Serial::set_compression_level_impl(int level) const {
  m_compression_level = pism::clip(level, 0, 9);
}

void NC4_Serial::def_var_impl(const std::string &name,
                              IO_Type nctype,
                              const std::vector<std::string> &dims) const {

  // use the parent class to define the variable
  NC_Serial::def_var_impl(name, nctype, dims);

  int stat = NC_NOERR;

  // set compression level (except for scalars and time-series)
  if (m_compression_level > 0 and dims.size() > 1) {
    if (m_rank == 0) {
      int varid = get_varid(name);

      stat = nc_def_var_deflate(m_file_id, varid, 0, 1, m_compression_level);

      // The NetCDF version used by PISM may not support compression.
      if (stat == NC_EINVAL) {
        stat = NC_NOERR;
      }
    }
  }

  MPI_Barrier(m_com);
  MPI_Bcast(&stat, 1, MPI_INT, 0, m_com);

  check(PISM_ERROR_LOCATION, stat);
}

} // end of namespace io
} // end of namespace pism
