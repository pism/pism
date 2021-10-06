// Copyright (C) 2004-2017, 2021 PISM Authors
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

#ifndef PART_GRID_THRESHOLD_THICKNESS_H
#define PART_GRID_THRESHOLD_THICKNESS_H

#include "pism/util/stencils.hh"

namespace pism {

double part_grid_threshold_thickness(stencils::Star<int> Mask,
                                     stencils::Star<double> thickness,
                                     stencils::Star<double> surface_elevation,
                                     double bed_elevation);

} // end of namespace pism

#endif /* PART_GRID_THRESHOLD_THICKNESS_H */
