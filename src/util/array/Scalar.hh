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

#ifndef PISM_ARRAY_SCALAR_H
#define PISM_ARRAY_SCALAR_H

#include <cmath>                // floor()

#include "pism/util/array/Array2D.hh"

namespace pism {
namespace array {

/** A class for storing and accessing scalar 2D fields. */
class Scalar : public Array2D<double> {
public:
  typedef std::shared_ptr<Scalar> Ptr;
  typedef std::shared_ptr<const Scalar> ConstPtr;

  Scalar(IceGrid::ConstPtr grid, const std::string &name);

  std::shared_ptr<Scalar> duplicate() const;

  inline int as_int(int i, int j) const;
  inline stencils::Star<int> star_int(int i, int j) const;
  inline stencils::Box<int> box_int(int i, int j) const;

protected:
  Scalar(IceGrid::ConstPtr grid, const std::string &name, int width);
};

inline int Scalar::as_int(int i, int j) const {
  const double &value = (*this)(i, j);
  return static_cast<int>(floor(value + 0.5));
}

/*!
 * Scalar 2D array supporting width=1 stencil computations
 */
class Scalar1 : public Scalar {
public:
  typedef std::shared_ptr<Scalar1> Ptr;
  typedef std::shared_ptr<const Scalar1> ConstPtr;

  Scalar1(IceGrid::ConstPtr grid, const std::string &name);
  using Array2D<double>::star;
  using Array2D<double>::box;
protected:
  Scalar1(IceGrid::ConstPtr grid, const std::string &name, int width);
};

/*!
 * Scalar 2D array supporting width=1 stencil computations
 */
class Scalar2 : public Scalar1 {
public:
  typedef std::shared_ptr<Scalar2> Ptr;
  typedef std::shared_ptr<const Scalar2> ConstPtr;

  Scalar2(IceGrid::ConstPtr grid, const std::string &name);
};

inline stencils::Star<int> Scalar::star_int(int i, int j) const {
  stencils::Star<int> result;

  result.ij = as_int(i,j);
  result.e =  as_int(i+1,j);
  result.w =  as_int(i-1,j);
  result.n =  as_int(i,j+1);
  result.s =  as_int(i,j-1);

  return result;
}

inline stencils::Box<int> Scalar::box_int(int i, int j) const {
  const int
      E = i + 1,
      W = i - 1,
      N = j + 1,
      S = j - 1;

  return {as_int(i, j), as_int(i, N), as_int(W, N), as_int(W, j), as_int(W, S),
          as_int(i, S), as_int(E, S), as_int(E, j), as_int(E, N)};
}
} // end of namespace array

// Finite-difference shortcuts. They may be slower than hard-coding FD approximations of x
// and y derivatives. Use with care.
double diff_x(const array::Scalar &array, int i, int j);
double diff_y(const array::Scalar &array, int i, int j);

// These take grid periodicity into account and use one-sided differences at domain edges.
double diff_x_p(const array::Scalar &array, int i, int j);
double diff_y_p(const array::Scalar &array, int i, int j);

double sum(const array::Scalar &input);
double min(const array::Scalar &input);
double max(const array::Scalar &input);
double absmax(const array::Scalar &input);

void apply_mask(const array::Scalar &M, double fill, array::Scalar &result);

void compute_magnitude(const array::Scalar &v_x,
                       const array::Scalar &v_y,
                       array::Scalar &result);

} // end of namespace pism

#endif /* PISM_ARRAY_SCALAR_H */
