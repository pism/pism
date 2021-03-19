// Copyright (C) 2012, 2014, 2015, 2016 PISM Authors
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
#include "pism/pism_config.hh"
#if (Pism_USE_CDIPIO==1)
extern "C"{
#include "cdi.h"
}
#endif

namespace pism {

static int pism_type_to_cdi_type(pism::IO_Type input) {
#if (Pism_USE_CDIPIO==1)
  switch (input) {
  case PISM_BYTE:
    return CDI_DATATYPE_INT8;
  case PISM_CHAR:
    return CDI_DATATYPE_TXT; //to be clarified
  case PISM_SHORT:
    return CDI_DATATYPE_INT16;
  case PISM_INT:
    return CDI_DATATYPE_FLT32;
  case PISM_FLOAT:
    return CDI_DATATYPE_FLT32;
  case PISM_DOUBLE:
    return CDI_DATATYPE_FLT64;
  default:
    return 0; //to be clarified
  }
#endif
}

static pism::IO_Type cdi_type_to_pism_type(int input) {
#if (Pism_USE_CDIPIO==1)
  switch (input) {
  case CDI_DATATYPE_INT8:
    return PISM_BYTE;
  case CDI_DATATYPE_TXT: // treat NC_CHAR and NC_STRING as equivalent
    return PISM_CHAR;
  case CDI_DATATYPE_INT16:
    return PISM_SHORT;
  case CDI_DATATYPE_INT32:
    return PISM_INT;
  case CDI_DATATYPE_FLT32:
    return PISM_FLOAT;
  case CDI_DATATYPE_FLT64:
    return PISM_DOUBLE;
  default:
    return PISM_NAT; //to be clarified
  }
#endif
}



} // end of namespace pism
