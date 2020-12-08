// Copyright (C) 2012, 2013, 2014, 2015, 2019 PISM Authors
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

#include "PISMNC4_Serial.hh"

// The following is a stupid kludge necessary to make NetCDF 4.x work in
// serial mode in an MPI program:
#ifndef MPI_INCLUDED
#define MPI_INCLUDED 1
#endif
#include <netcdf.h>

namespace pism {
namespace io {

int NC4_Serial::open_impl(const std::string &fname, IO_Mode mode, const std::map<std::string, int> &varsi, int FileID) {
  int open_mode = mode == PISM_READONLY ? NC_NOWRITE : NC_WRITE;

  int stat = nc_open(fname.c_str(), open_mode, &m_file_id);
  check(PISM_ERROR_LOCATION, stat);
}

int NC4_Serial::create_impl(const std::string &fname, int FileID) {
  int stat = nc_create(fname.c_str(), NC_NETCDF4, &m_file_id);
  check(PISM_ERROR_LOCATION, stat);
}

} // end of namespace io
} // end of namespace pism
