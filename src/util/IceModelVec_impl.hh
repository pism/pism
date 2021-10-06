/* Copyright (C) 2020, 2021 PISM Authors
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
#ifndef ICEMODELVEC_IMPL_H
#define ICEMODELVEC_IMPL_H

#include <map>
#include <memory>
#include <string>
#include <vector>

#include <gsl/gsl_interp.h>

#include "pism/util/petscwrappers/Vec.hh"
#include "pism/util/interpolation.hh"

namespace pism {

struct IceModelVec::Impl {
  Impl() {
    access_counter = 0;

    da.reset();

    da_stencil_width = 1;
    dof = 1;
    begin_access_use_dof = false;

    ghosted = true;

    report_range = true;

    name = "uninitialized variable";

    zlevels = {0.0};

    state_counter = 0;
    interpolation_type = LINEAR;

    bsearch_accel = nullptr;
  }
  //! If true, report range when regridding.
  bool report_range;

  //! The array itself
  //!
  //! Note: do not access this directly (via `m_impl->v`). Use `vec()` instead.
  petsc::Vec v;

  //! Name of the field. In general this is *not* the name of the corresponding NetCDF
  //! variable.
  std::string name;

  //! Metadata (NetCDF variable attributes)
  std::vector<SpatialVariableMetadata> metadata;

  //! The computational grid
  IceGrid::ConstPtr grid;

  //! number of "degrees of freedom" per grid point
  unsigned int dof;

  //! stencil width supported by the DA
  unsigned int da_stencil_width;

  //! true if this IceModelVec is ghosted
  bool ghosted;

  //! distributed mesh manager (DM)
  std::shared_ptr<petsc::DM> da;

  //! If true, use DMDAVecGetArrayDOF() in begin_access()
  bool begin_access_use_dof;

  //! Map plane viewers. It is a map because a temporary IceModelVec can be used to view
  //! different quantities
  std::map<std::string,std::shared_ptr<petsc::Viewer> > map_viewers;

  // used in begin_access() and end_access()
  int access_counter;

  //! Internal IceModelVec "revision number"
  int state_counter;

  // 2D Interpolation type (used by regrid())
  InterpolationType interpolation_type;

  //! Vertical levels (for 3D fields)
  std::vector<double> zlevels;

  // binary search accelerator (used for interpolation in a column in 3D fields)
  gsl_interp_accel *bsearch_accel;
};

} // end of namespace pism

#endif /* ICEMODELVEC_IMPL_H */
