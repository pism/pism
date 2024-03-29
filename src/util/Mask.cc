// Copyright (C) 2011, 2014, 2015, 2016, 2017, 2018, 2022, 2023 Constantine Khroulev and David Maxwell
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

#include <cassert>

#include "pism/util/Mask.hh"
#include "pism/util/Grid.hh"
#include "pism/util/array/Scalar.hh"

namespace pism {

void GeometryCalculator::compute(const array::Scalar& sea_level,
                                 const array::Scalar& bed,
                                 const array::Scalar& thickness,
                                 array::Scalar& out_mask,
                                 array::Scalar& out_surface) const {
  compute_mask(sea_level, bed, thickness, out_mask);
  compute_surface(sea_level, bed, thickness, out_surface);
}

void GeometryCalculator::compute_mask(const array::Scalar &sea_level,
                                      const array::Scalar &bed,
                                      const array::Scalar &thickness,
                                      array::Scalar &result) const {
  array::AccessScope list{&sea_level, &bed, &thickness, &result};

  const Grid &grid = *bed.grid();

  const unsigned int stencil = result.stencil_width();
  assert(sea_level.stencil_width() >= stencil);
  assert(bed.stencil_width()       >= stencil);
  assert(thickness.stencil_width() >= stencil);

  for (auto p = grid.points(stencil); p; p.next()) {
    const int i = p.i(), j = p.j();

    result(i,j) = this->mask(sea_level(i, j), bed(i, j), thickness(i, j));
  }
}

void GeometryCalculator::compute_surface(const array::Scalar &sea_level,
                                         const array::Scalar &bed,
                                         const array::Scalar &thickness,
                                         array::Scalar &result) const {
  array::AccessScope list{&sea_level, &bed, &thickness, &result};

  const Grid &grid = *bed.grid();

  const unsigned int stencil = result.stencil_width();
  assert(sea_level.stencil_width() >= stencil);
  assert(bed.stencil_width()       >= stencil);
  assert(thickness.stencil_width() >= stencil);

  for (auto p = grid.points(stencil); p; p.next()) {
    const int i = p.i(), j = p.j();

    result(i,j) = this->surface(sea_level(i, j), bed(i, j), thickness(i, j));
  }
}

} // end of namespace pism
