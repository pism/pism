/* Copyright (C) 2015, 2016, 2017, 2019 PISM Authors
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

#ifndef _ICEMODELVEC_INLINE_H_
#define _ICEMODELVEC_INLINE_H_

// This header is included by iceModelVec.hh. Do not include it
// manually.

namespace pism {

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

inline stencils::Star<double> IceModelVec2Stag::star(int i, int j) const {
  const IceModelVec2Stag &self = *this;

  stencils::Star<double> result;

  result.ij = 0.0;             // has no meaning in this context
  result.e =  self(i, j, 0);
  result.w =  self(i-1, j, 0);
  result.n =  self(i, j, 1);
  result.s =  self(i, j-1, 1);

  return result;
}

inline int IceModelVec2Int::as_int(int i, int j) const {
#if (Pism_DEBUG==1)
  check_array_indices(i, j, 0);
#endif
  const double **a = (const double**) m_array;
  return static_cast<int>(floor(a[j][i] + 0.5));
}

inline stencils::Star<int> IceModelVec2Int::star(int i, int j) const {
  stencils::Star<int> result;

  result.ij = as_int(i,j);
  result.e =  as_int(i+1,j);
  result.w =  as_int(i-1,j);
  result.n =  as_int(i,j+1);
  result.s =  as_int(i,j-1);

  return result;
}

inline stencils::Box<int> IceModelVec2Int::box(int i, int j) const {
  const IceModelVec2Int &x = *this;

  const int
      E = i + 1,
      W = i - 1,
      N = j + 1,
      S = j - 1;

  return {x.as_int(i, j), x.as_int(i, N), x.as_int(W, N), x.as_int(W, j), x.as_int(W, S),
          x.as_int(i, S), x.as_int(E, S), x.as_int(E, j), x.as_int(E, N)};
}

inline double& IceModelVec3::operator() (int i, int j, int k) {
#if (Pism_DEBUG==1)
  check_array_indices(i, j, k);
#endif
  return static_cast<double***>(m_array)[j][i][k];
}

inline const double& IceModelVec3::operator() (int i, int j, int k) const {
#if (Pism_DEBUG==1)
  check_array_indices(i, j, k);
#endif
  return static_cast<double***>(m_array)[j][i][k];
}

} // end of namespace pism

#endif /* _ICEMODELVEC_INLINE_H_ */
