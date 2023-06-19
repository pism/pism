// Copyright (C) 2012, 2014, 2015, 2016, 2023 PISM Authors
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


//! Convert PISM's IO types into NetCDF types and back. Note that NC_* may be
//! macros, so you need to include the appropriate NetCDF header first.
namespace pism {

static nc_type pism_type_to_nc_type(pism::io::Type input) {
  switch (input) {
  case io::PISM_BYTE:
    return NC_BYTE;
  case io::PISM_CHAR:
    return NC_CHAR;
  case io::PISM_SHORT:
    return NC_SHORT;
  case io::PISM_INT:
    return NC_INT;
  case io::PISM_FLOAT:
    return NC_FLOAT;
  case io::PISM_DOUBLE:
    return NC_DOUBLE;
  default:
    return NC_NAT;
  }
}

static pism::io::Type nc_type_to_pism_type(int input) {
  switch (input) {
  case NC_BYTE:
    return io::PISM_BYTE;
  case NC_CHAR:
  case NC_STRING:               // treat NC_CHAR and NC_STRING as equivalent
    return io::PISM_CHAR;
  case NC_SHORT:
    return io::PISM_SHORT;
  case NC_INT:
    return io::PISM_INT;
  case NC_FLOAT:
    return io::PISM_FLOAT;
  case NC_DOUBLE:
    return io::PISM_DOUBLE;
  default:
    return io::PISM_NAT;
  }
}

} // end of namespace pism
