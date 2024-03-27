// Copyright (C) 2007--2011, 2013, 2014, 2015, 2017, 2018, 2023, 2024 Jed Brown, Ed Bueler and Constantine Khroulev
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

#ifndef PISM_LOCALINTERPCTX_H
#define PISM_LOCALINTERPCTX_H

#include <array>
#include <vector>
#include <memory>

namespace pism {

enum InterpolationType : int;
class Interpolation1D;

class Grid;

namespace grid {
class InputGridInfo;
}

//! The "local interpolation context" describes the processor's part of the source NetCDF
//! file (for regridding).
/*! The local interpolation context contains the details of how the processor's block of
  the new computational domain fits into the domain of the netCDF file. Note that each
  vertical column of the grid is owned by exactly one processor.

  For any particular dimension, we have a new computational domain \f$[a,b]\f$ with
  spacing \f$h\f$ so there are \f$n = (b - a) / h\f$ interior cells, indexed by \f$\{i_0, \dots, i_n\}\f$.
  The local processor owns a range \f$\{i_m, \dots, i_{m'}\}\f$.  Suppose the netCDF file has
  domain \f$[A,B]\f$, spacing \f$H\f$, and \f$N = (B - A) / H\f$ cells.  In order to interpolate 
  onto these points, we need the indices \f$\{I_m, \dots, I_{m'}\}\f$ of the netCDF file so that

  \f[  [x(i_m), x(i_{m'})] \quad \text{is a subset of} \quad  [x(I_m), x(I_{m'})]  \f]

  The arrays `start` and `count` have 4 integer entries, corresponding to the dimensions
  `t, x, y, z(zb)`.
*/
class LocalInterpCtx {
public:
  LocalInterpCtx(const grid::InputGridInfo &input_grid, const Grid &internal_grid,
                 InterpolationType type);
  LocalInterpCtx(const grid::InputGridInfo &input_grid, const Grid &internal_grid,
                 const std::vector<double> &z_internal, InterpolationType type);

  int buffer_size() const;

  // Indices in netCDF file.
  std::array<int, 4> start, count;
  // indexes and coefficients for 1D linear interpolation
  std::shared_ptr<Interpolation1D> x, y, z;
};

} // end of namespace pism

#endif // PISM_LOCALINTERPCTX_H
