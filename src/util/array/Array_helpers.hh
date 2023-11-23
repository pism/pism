// Copyright (C) 2011, 2013, 2014, 2016, 2017, 2020, 2021, 2022, 2023 PISM Authors
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

#ifndef PISM_ARRAY_HELPERS_H
#define PISM_ARRAY_HELPERS_H

#include "pism/util/Grid.hh"

namespace pism {

namespace array {

namespace details {

//! \brief Computes result = x + alpha * y, where x, y, and z are 2D
//! Arrays (scalar or vector).
/*!
 */
template <class V>
void add(const V &x, double alpha, const V &y, V &result, bool scatter = true) {

  array::AccessScope list{ &x, &y, &result };
  for (auto p = result.grid()->points(); p; p.next()) {
    const int i = p.i(), j = p.j();

    result(i, j) = x(i, j) + y(i, j) * alpha;
  }

  if (scatter) {
    result.update_ghosts();
  }

  result.inc_state_counter();
}

template <class V>
void copy(const V &input, V &result, bool scatter = true) {

  array::AccessScope list{ &input, &result };

  for (auto p = result.grid()->points(); p; p.next()) {
    const int i = p.i(), j = p.j();

    result(i, j) = input(i, j);
  }

  if (scatter) {
    result.update_ghosts();
  }

  result.inc_state_counter();
}

} // namespace details

} // namespace array

} // end of namespace pism

#endif /* PISM_ARRAY_HELPERS_H */
