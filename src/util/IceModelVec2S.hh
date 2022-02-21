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

#include "pism/util/iceModelVec.hh"

namespace pism {

/** A class for storing and accessing scalar 2D fields.
    IceModelVec2S is just IceModelVec2 with "dof == 1" */
class IceModelVec2S : public IceModelVec {
public:
  IceModelVec2S(IceGrid::ConstPtr grid, const std::string &name,
                IceModelVecKind ghostedp, int width = 1);

  typedef std::shared_ptr<IceModelVec2S> Ptr;
  typedef std::shared_ptr<const IceModelVec2S> ConstPtr;

  // does not need a copy constructor, because it does not add any new data members
  void copy_from(const IceModelVec2S &source);
  double** array();
  double const* const* array() const;
  void add(double alpha, const IceModelVec2S &x);
  void add(double alpha, const IceModelVec2S &x, IceModelVec2S &result) const;

  //! Provides access (both read and write) to the internal double array.
  /*!
    Note that i corresponds to the x direction and j to the y.
  */
  inline double& operator() (int i, int j);
  inline const double& operator()(int i, int j) const;
  inline stencils::Star<double> star(int i, int j) const;
  inline stencils::Box<double> box(int i, int j) const;
};

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

inline double& IceModelVec2S::operator() (int i, int j) {
#if (Pism_DEBUG==1)
  check_array_indices(i, j, 0);
#endif
  return static_cast<double**>(m_array)[j][i];
}

inline const double& IceModelVec2S::operator()(int i, int j) const {
#if (Pism_DEBUG==1)
  check_array_indices(i, j, 0);
#endif
  return static_cast<double**>(m_array)[j][i];
}

inline stencils::Star<double> IceModelVec2S::star(int i, int j) const {
  const IceModelVec2S &self = *this;

  stencils::Star<double> result;
  result.ij = self(i,j);
  result.e =  self(i+1,j);
  result.w =  self(i-1,j);
  result.n =  self(i,j+1);
  result.s =  self(i,j-1);

  return result;
}

inline stencils::Box<double> IceModelVec2S::box(int i, int j) const {
  const IceModelVec2S &x = *this;

  const int
      E = i + 1,
      W = i - 1,
      N = j + 1,
      S = j - 1;

  return {x(i, j), x(i, N), x(W, N), x(W, j), x(W, S), x(i, S), x(E, S), x(E, j), x(E, N)};
}

} // end of namespace pism

#endif /* PISM_ICEMODELVEC2S_H */
