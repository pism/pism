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

#include "PISMNC4_Serial.hh"

// The following is a stupid kludge necessary to make NetCDF 4.x work in
// serial mode in an MPI program:
#ifndef MPI_INCLUDED
#define MPI_INCLUDED 1
#endif
#include <netcdf.h>

int PISMNC4_Serial::open(string fname, int mode) {
  int stat;

  m_filename = fname;

  stat = nc_open(m_filename.c_str(), mode, &ncid);

  define_mode = false;

  return stat;
}

int PISMNC4_Serial::create(string fname) {
  int stat;

  m_filename = fname;

  stat = nc_create(m_filename.c_str(), NC_NETCDF4, &ncid);

  define_mode = true;

  return stat;
}
