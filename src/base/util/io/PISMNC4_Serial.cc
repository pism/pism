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

#include <netcdf.h>
#include "PISMNC4_Serial.hh"

int PISMNC4_Serial::open(string fname, int mode) {
  int stat;

  filename = fname;

  stat = nc_open(filename.c_str(), mode, &ncid);

  define_mode = false;

  return stat;
}

int PISMNC4_Serial::create(string fname) {
  int stat;

  filename = fname;

  stat = nc_create(filename.c_str(), NC_NETCDF4, &ncid);

  define_mode = true;

  return stat;
}
