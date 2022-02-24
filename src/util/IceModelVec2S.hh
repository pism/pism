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

#ifndef PISM_ICEMODELVEC2S_H
#define PISM_ICEMODELVEC2S_H

#include <cmath>                // floor()

#include "pism/util/IceModelVec2.hh"

namespace pism {

/** A class for storing and accessing scalar 2D fields.
    IceModelVec2S is just IceModelVec2 with "dof == 1" */
class IceModelVec2S : public IceModelVec2<double> {
public:
  typedef std::shared_ptr<IceModelVec2S> Ptr;
  typedef std::shared_ptr<const IceModelVec2S> ConstPtr;

  IceModelVec2S(IceGrid::ConstPtr grid, const std::string &name);

  inline int as_int(int i, int j) const;
  inline stencils::Star<int> star_int(int i, int j) const;
  inline stencils::Box<int> box_int(int i, int j) const;

protected:
  IceModelVec2S(IceGrid::ConstPtr grid, const std::string &name,
                int width);
};

inline int IceModelVec2S::as_int(int i, int j) const {
  const double &value = (*this)(i, j);
  return static_cast<int>(floor(value + 0.5));
}

template<int width>
class Array2SGhosted : public IceModelVec2S {
public:
  Array2SGhosted(IceGrid::ConstPtr grid, const std::string &name)
  : IceModelVec2S(grid, name, width) {
    // empty
  }

  // Allow implicit casting to a reference to an array with a smaller stencil width:
  template <int smaller_width>
  operator Array2SGhosted<smaller_width>&() {
    static_assert(smaller_width < width, "insufficient stencil width");
    return *this;
  }
};

inline stencils::Star<int> IceModelVec2S::star_int(int i, int j) const {
  stencils::Star<int> result;

  result.ij = as_int(i,j);
  result.e =  as_int(i+1,j);
  result.w =  as_int(i-1,j);
  result.n =  as_int(i,j+1);
  result.s =  as_int(i,j-1);

  return result;
}

inline stencils::Box<int> IceModelVec2S::box_int(int i, int j) const {
  const int
      E = i + 1,
      W = i - 1,
      N = j + 1,
      S = j - 1;

  return {as_int(i, j), as_int(i, N), as_int(W, N), as_int(W, j), as_int(W, S),
          as_int(i, S), as_int(E, S), as_int(E, j), as_int(E, N)};
}

std::shared_ptr<IceModelVec2S> duplicate(const IceModelVec2S &source);

// Finite-difference shortcuts. They may be slower than hard-coding FD approximations of x
// and y derivatives. Use with care.
double diff_x(const IceModelVec2S &array, int i, int j);
double diff_y(const IceModelVec2S &array, int i, int j);

// These take grid periodicity into account and use one-sided differences at domain edges.
double diff_x_p(const IceModelVec2S &array, int i, int j);
double diff_y_p(const IceModelVec2S &array, int i, int j);

double sum(const IceModelVec2S &input);
double min(const IceModelVec2S &input);
double max(const IceModelVec2S &input);
double absmax(const IceModelVec2S &input);

void apply_mask(const IceModelVec2S &M, double fill, IceModelVec2S &result);

void compute_magnitude(const IceModelVec2S &v_x,
                       const IceModelVec2S &v_y,
                       IceModelVec2S &result);

} // end of namespace pism

#endif /* PISM_ICEMODELVEC2S_H */
