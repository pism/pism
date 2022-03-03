/* Copyright (C) 2020, 2021, 2022 PISM Authors
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
#ifndef PISM_ARRAY2D_H
#define PISM_ARRAY2D_H

#include "Array.hh"
#include "Array_helpers.hh"

namespace pism {

namespace array {

//! A storage vector combining related fields in a struct
template<typename T>
class Array2D : public Array {
public:
  using value_type = T;

  Array2D(IceGrid::ConstPtr grid, const std::string &short_name,
          IceModelVecKind ghostedp, unsigned int stencil_width = 1)
    : Array(grid, short_name, ghostedp,
            sizeof(T) / sizeof(double), stencil_width, {0.0}) {
    set_begin_access_use_dof(false);
  }

  T** array() {
    return reinterpret_cast<T**>(m_array);
  }

  T const* const* array() const {
    return reinterpret_cast<T const* const*>(m_array);
  }

  inline T& operator()(int i, int j) {
#if (Pism_DEBUG==1)
    check_array_indices(i, j, 0);
#endif
    return static_cast<T**>(m_array)[j][i];
  }

  inline const T& operator()(int i, int j) const {
#if (Pism_DEBUG==1)
    check_array_indices(i, j, 0);
#endif
    return static_cast<T**>(m_array)[j][i];
  }

  void add(double alpha, const Array2D<T> &x) {
    vec::add(*this, alpha, x, *this);
  }

  void add(double alpha, const Array2D<T> &x, Array2D<T> &result) const {
    vec::add(*this, alpha, x, result);
  }

  void copy_from(const Array2D<T> &source) {
    return vec::copy(source, *this);
  }

protected:

  inline stencils::Star<T> star(int i, int j) const {
    const auto &self = *this;

    stencils::Star<T> result;

    result.ij = self(i,j);
    result.e =  self(i+1,j);
    result.w =  self(i-1,j);
    result.n =  self(i,j+1);
    result.s =  self(i,j-1);

    return result;
  }

  inline stencils::Box<T> box(int i, int j) const {
    const auto &x = *this;

    const int
      E = i + 1,
      W = i - 1,
      N = j + 1,
      S = j - 1;

    return {x(i, j), x(i, N), x(W, N), x(W, j), x(W, S),
            x(i, S), x(E, S), x(E, j), x(E, N)};
  }
};

} // end of namespace array
} // end of namespace pism

#endif /* PISM_ARRAY2D_H */
