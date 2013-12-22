// Copyright (C) 2007--2011, 2013, 2014 Jed Brown, Ed Bueler and Constantine Khroulev
//
// This file is part of Pism.
//
// Pism is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
// version.
//
// Pism is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License
// along with Pism; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

#ifndef __lic_hh
#define __lic_hh

#include "IceGrid.hh"

//! The "local interpolation context" describes the processor's part of the source NetCDF file (for regridding).
/*!
The local interpolation context contains the details of how the processor's block
of the new computational domain fits into the domain of the netCDF file.  Note that each vertical column 
of the grid is owned by exactly one processor.

For any particular dimension, we have a new computational domain \f$[a,b]\f$ with
spacing \f$h\f$ so there are \f$n = (b - a) / h\f$ interior cells, indexed by \f$\{i_0, \dots, i_n\}\f$.
The local processor owns a range \f$\{i_m, \dots, i_{m'}\}\f$.  Suppose the netCDF file has
domain \f$[A,B]\f$, spacing \f$H\f$, and \f$N = (B - A) / H\f$ cells.  In order to interpolate 
onto these points, we need the indices \f$\{I_m, \dots, I_{m'}\}\f$ of the netCDF file so that

  \f[  [x(i_m), x(i_{m'})] \quad \text{is a subset of} \quad  [x(I_m), x(I_{m'})]  \f]

The arrays `start` and `count` have 4 integer entries, corresponding to the dimensions
\f$t, x, y, z(zb)\f$.
 */
class LocalInterpCtx {
public:
  unsigned int start[4], count[4]; // Indices in netCDF file.
  std::vector<int> x_left, x_right, y_left, y_right; // neighbors
  std::vector<double> x_alpha, y_alpha;
  double *a;                       //!< temporary buffer
  int a_len;                       //!< the size of the buffer
  std::vector<double> zlevels;          //!< input z levels
  bool report_range;
  MPI_Comm com;			//!< MPI Communicator (for printing, mostly)
  PetscMPIInt rank;		//!< MPI rank, to allocate a_raw on proc 0 only

public:
  LocalInterpCtx(grid_info g, const IceGrid &grid, PetscReal z_min, PetscReal z_max);
  ~LocalInterpCtx();
  PetscErrorCode printArray();
private:
 void print_grid_info(grid_info g, PISMUnitSystem s, int threshold);
};

#endif // __lic_hh
