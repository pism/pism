/* Copyright (C) 2020 PISM Authors
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
#ifndef PISM_ICEMODELVEC2_STRUCT_H
#define PISM_ICEMODELVEC2_STRUCT_H

#include "pism/util/iceModelVec.hh"

namespace pism {

//! A storage vector combining related fields in a struct
template<typename T>
class IceModelVec2Struct : public IceModelVec2 {
public:
  IceModelVec2Struct(IceGrid::ConstPtr grid, const std::string &short_name,
                     IceModelVecKind ghostedp, unsigned int stencil_width = 1)
    : IceModelVec2(grid, short_name, ghostedp, stencil_width,
                   sizeof(T) / sizeof(double)) {
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

  inline StarStencil<T> star(int i, int j) const {
    const auto &self = *this;

    StarStencil<T> result;

    result.ij = self(i,j);
    result.e =  self(i+1,j);
    result.w =  self(i-1,j);
    result.n =  self(i,j+1);
    result.s =  self(i,j-1);

    return result;
  }

  inline BoxStencil<T> box(int i, int j) const {
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

} // end of namespace pism

#endif /* PISM_ICEMODELVEC2_STRUCT_H */
