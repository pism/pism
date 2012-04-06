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


//! Convert PISM's IO types into NetCDF types and back. Note that NC_* may be
//! macros, so you need to include the appropriate NetCDF header first.

static nc_type pism_type_to_nc_type(PISM_IO_Type input) {
  switch (input) {
  case PISM_BYTE:
    return NC_BYTE;
  case PISM_CHAR:
    return NC_CHAR;
  case PISM_SHORT:
    return NC_SHORT;
  case PISM_INT:
    return NC_INT;
  case PISM_FLOAT:
    return NC_FLOAT;
  case PISM_DOUBLE:
    return NC_DOUBLE;
  default:
    return NC_NAT;
  }
}

static PISM_IO_Type nc_type_to_pism_type(nc_type input) {
  switch (input) {
  case NC_BYTE:
    return PISM_BYTE;
  case NC_CHAR:
    return PISM_CHAR;
  case NC_SHORT:
    return PISM_SHORT;
  case NC_INT:
    return PISM_INT;
  case NC_FLOAT:
    return PISM_FLOAT;
  case NC_DOUBLE:
    return PISM_DOUBLE;
  default:
    return PISM_NAT;
  }
}
