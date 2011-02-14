// Copyright (C) 2007--2011 Jed Brown, Ed Bueler and Constantine Khroulev
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

#include "grid.hh"

//! \brief Contains parameters of an input file grid.
class grid_info {
public:
  // dimension lengths
  int t_len, x_len, y_len, z_len, zb_len;
  double time,			//!< current time (seconds)
    x0,				//!< x-coordinate of the grid center
    y0,				//!< y-coordinate of the grid center
    Lx,				//!< half-width in the X-direction
    Ly,				//!< half-width in the Y-direction
    x_min,			//!< [x_min, x_max] is the X extent of the grid
    x_max,			//!< [x_min, x_max] is the X extent of the grid
    y_min,			//!< [y_min, y_max] is the Y extent of the grid
    y_max,			//!< [y_min, y_max] is the Y extent of the grid
    zb_min,			//!< minimal value of the zb dimension
    z_max;			//!< maximal value of the z dimension
  grid_info();
  PetscErrorCode print(MPI_Comm com, int threshold = 3);
};

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
  double *a;		     //!< temporary buffer
  int a_len;		     //!< the size of the buffer
  int nz,		     //!< number of z-levels
    nzb;		     //!< number of zb-levels 
  double *zlevs, *zblevs;
  bool no_regrid_ice, no_regrid_bedrock, report_range;
  MPI_Comm com;			//!< MPI Communicator (for printing, mostly)
  PetscMPIInt rank;		//!< MPI rank, to allocate a_raw on proc 0 only

public:
  LocalInterpCtx(grid_info g,
                 const double zlevsIN[], const double zblevsIN[], IceGrid &grid);
  ~LocalInterpCtx();
  int kBelowHeight(const double height);
  int kbBelowHeight(const double elevation);
  PetscErrorCode printGrid();
  PetscErrorCode printArray();
protected:
};

#endif // __lic_hh
