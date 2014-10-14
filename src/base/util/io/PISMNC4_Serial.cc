// Copyright (C) 2012, 2013, 2014 PISM Authors
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

int NC4_Serial::integer_open_mode(IO_Mode input) const {
  if (input == PISM_READONLY) {
    return NC_NOWRITE;
  } else {
    return NC_WRITE;
  }
}

int NC4_Serial::open_impl(const std::string &fname, IO_Mode mode) {
  int stat;

  m_filename = fname;

  int nc_mode = integer_open_mode(mode);
  stat = nc_open(m_filename.c_str(), nc_mode, &m_file_id);

  m_define_mode = false;

  return stat;
}

int NC4_Serial::create_impl(const std::string &fname) {
  int stat;

  m_filename = fname;

  stat = nc_create(m_filename.c_str(), NC_NETCDF4, &m_file_id);

  m_define_mode = true;

  return stat;
}

} // end of namespace pism
