/* Copyright (C) 2022 PISM Authors
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

#ifndef PISM_ICEMODELVEC2INT_H
#define PISM_ICEMODELVEC2INT_H

#include <cmath>                // floor()

#include "IceModelVec2S.hh"

namespace pism {

//! \brief A simple class "hiding" the fact that the mask is stored as
//! floating-point scalars (instead of integers).
class IceModelVec2Int : public IceModelVec2S {
public:
  IceModelVec2Int(IceGrid::ConstPtr grid, const std::string &name,
                  IceModelVecKind ghostedp, int width = 1);

  typedef std::shared_ptr<IceModelVec2Int> Ptr;
  typedef std::shared_ptr<const IceModelVec2Int> ConstPtr;

  inline int as_int(int i, int j) const;
  inline stencils::Star<int> star_int(int i, int j) const;
  inline stencils::Box<int> box_int(int i, int j) const;
};

inline int IceModelVec2Int::as_int(int i, int j) const {
  const double &value = (*this)(i, j);
  return static_cast<int>(floor(value + 0.5));
}

inline stencils::Star<int> IceModelVec2Int::star_int(int i, int j) const {
  stencils::Star<int> result;

  result.ij = as_int(i,j);
  result.e =  as_int(i+1,j);
  result.w =  as_int(i-1,j);
  result.n =  as_int(i,j+1);
  result.s =  as_int(i,j-1);

  return result;
}

inline stencils::Box<int> IceModelVec2Int::box_int(int i, int j) const {
  const int
      E = i + 1,
      W = i - 1,
      N = j + 1,
      S = j - 1;

  return {as_int(i, j), as_int(i, N), as_int(W, N), as_int(W, j), as_int(W, S),
          as_int(i, S), as_int(E, S), as_int(E, j), as_int(E, N)};
}

} // end of namespace pism

#endif /* PISM_ICEMODELVEC2INT_H */
