/* Copyright (C) 2014, 2015, 2018, 2019, 2020, 2021 PISM Authors
 *
 * This file is part of PISM.
 *
 * PISM is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option) any later
 * version.
 *
 * PISM is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PISM; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

#ifndef _IO_FLAGS_H_
#define _IO_FLAGS_H_

namespace pism {

// I/O Flags used by File and NCFile. They are used in both interfaces,
// but I want to be able to create Python wrappers for File without
// exposing NCFile, and NCFile should compile without File, so it does
// not belong in either File.hh or NCFile.hh.

// This is a subset of NetCDF data-types.
enum IO_Type {
  PISM_NAT    = 0,              /* NAT = 'Not A Type' (c.f. NaN) */
  PISM_BYTE   = 1,              /* signed 1 byte integer */
  PISM_CHAR   = 2,              /* ISO/ASCII character */
  PISM_SHORT  = 3,              /* signed 2 byte integer */
  PISM_INT    = 4,              /* signed 4 byte integer */
  PISM_FLOAT  = 5,              /* single precision floating point number */
  PISM_DOUBLE = 6               /* double precision floating point number */
};

enum IO_Backend {PISM_GUESS, PISM_NETCDF3, PISM_NETCDF4_SERIAL, PISM_NETCDF4_PARALLEL,
                 PISM_PNETCDF, PISM_PIO_PNETCDF, PISM_PIO_NETCDF, PISM_PIO_NETCDF4C, PISM_PIO_NETCDF4P,
                 PISM_CDI};

// This is a subset of NetCDF file modes. Use values that don't match
// NetCDF flags so that we can detect errors caused by passing these
// straight to NetCDF.
enum IO_Mode {
  //! open an existing file for reading only
  PISM_READONLY          = 7,
  //! open an existing file for reading and writing
  PISM_READWRITE         = 8,
  //! create a file for writing, overwrite if present
  PISM_READWRITE_CLOBBER = 9,
  //! create a file for writing, move foo.nc to foo.nc~ if present
  PISM_READWRITE_MOVE    = 10
};

// This is the special value corresponding to the "unlimited" dimension length.
// Gets cast to "int", so it should match the value used by NetCDF.
enum Dim_Length {
  PISM_UNLIMITED = 0
};

// "Fill" mode. Gets cast to "int", so it should match values used by NetCDF.
enum Fill_Mode {
  PISM_FILL   = 0,
  PISM_NOFILL = 0x100
};

enum RegriddingFlag {OPTIONAL, OPTIONAL_FILL_MISSING, CRITICAL, CRITICAL_FILL_MISSING};

} // end of namespace pism

#endif /* _IO_FLAGS_H_ */
