// Copyright (C) 2007, 2008 Jed Brown, Ed Bueler and Constantine Khroulev
//
// This file is part of Pism.
//
// Pism is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 2 of the License, or (at your option) any later
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

#include <netcdf.h>
#include <petscmat.h>
#include "grid.hh"


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

In the code, \c xbdy\f$[2] = \{x(i_m), x(i_{m'})\}\f$.  We have obtained the netCDF bounds 
\f$x(I_0)\f$ and \f$x(I_N)\f$ in the \c bdy array and the
number of elements \f$N+1\f$ in the \c dim array.

This is a summary of the fields of a \c LocalInterpCtx, and a translation between scalars in the previous
paragraphs and the fields:
 - \f$I_m = \operatorname{floor}((x(i_m) - A) / H)\f$ is called \c start
 - \f$I_{m'} - I_m = \operatorname{ceil}((x(i_{m'}) - X(I_m)) / H\f$ is \c count - 1
 - \f$X(I_m)\f$ is called \c fstart
 - \f$H\f$ is called \c delta

\c delta and \c fstart are not used in the vertical dimension because the spacing is not 
generally equal.  Rather, \c zlevs and \c zblevs contain the needed information.

The arrays \c start and \c count have 5 integer entries, corresponding to the dimensions
\f$t, x, y, z, zb\f$.
 */
class LocalInterpCtx {
public:
  double fstart[3], delta[3];
  int start[5], count[5];    // Indices in netCDF file.
  double *a;
  int a_len;
  int nz, nzb;
  double *zlevs, *zblevs;
  bool regrid_2d_only, no_regrid_bedrock;

public:
  LocalInterpCtx(const size_t dim[], const double bdy[],
                 const double zlevsIN[], const double zblevsIN[], IceGrid &grid);
  ~LocalInterpCtx();
  int kBelowHeight(const double height, MPI_Comm com);
  int kbBelowHeight(const double elevation, MPI_Comm com);
  PetscErrorCode printGrid(MPI_Comm com);
  PetscErrorCode printArray(MPI_Comm com);
protected:
};

#endif // __lic_hh
