/* Copyright (C) 2015, 2016, 2017 PISM Authors
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

inline double& IceModelVec2::operator() (int i, int j, int k) {
#if (PISM_DEBUG==1)
  check_array_indices(i, j, k);
#endif
  return static_cast<double***>(m_array)[j][i][k];
}

inline const double& IceModelVec2::operator() (int i, int j, int k) const {
#if (PISM_DEBUG==1)
  check_array_indices(i, j, k);
#endif
  return static_cast<double***>(m_array)[j][i][k];
}

inline double& IceModelVec2S::operator() (int i, int j) {
#if (PISM_DEBUG==1)
  check_array_indices(i, j, 0);
#endif
  return static_cast<double**>(m_array)[j][i];
}

inline const double& IceModelVec2S::operator()(int i, int j) const {
#if (PISM_DEBUG==1)
  check_array_indices(i, j, 0);
#endif
  return static_cast<double**>(m_array)[j][i];
}

inline StarStencil<double> IceModelVec2S::star(int i, int j) const {
  const IceModelVec2S &self = *this;

  StarStencil<double> result;
  result.ij = self(i,j);
  result.e =  self(i+1,j);
  result.w =  self(i-1,j);
  result.n =  self(i,j+1);
  result.s =  self(i,j-1);

  return result;
}

inline StarStencil<double> IceModelVec2Stag::star(int i, int j) const {
  const IceModelVec2Stag &self = *this;

  StarStencil<double> result;

  result.ij = 0.0;             // has no meaning in this context
  result.e =  self(i, j, 0);
  result.w =  self(i-1, j, 0);
  result.n =  self(i, j, 1);
  result.s =  self(i, j-1, 1);

  return result;
}

inline int IceModelVec2Int::as_int(int i, int j) const {
#if (PISM_DEBUG==1)
  check_array_indices(i, j, 0);
#endif
  const double **a = (const double**) m_array;
  return static_cast<int>(floor(a[j][i] + 0.5));
}

inline StarStencil<int> IceModelVec2Int::int_star(int i, int j) const {
  StarStencil<int> result;

  result.ij = as_int(i,j);
  result.e =  as_int(i+1,j);
  result.w =  as_int(i-1,j);
  result.n =  as_int(i,j+1);
  result.s =  as_int(i,j-1);

  return result;
}


inline Vector2& IceModelVec2V::operator()(int i, int j) {
#if (PISM_DEBUG==1)
  check_array_indices(i, j, 0);
#endif
  return static_cast<Vector2**>(m_array)[j][i];
}

inline const Vector2& IceModelVec2V::operator()(int i, int j) const {
#if (PISM_DEBUG==1)
  check_array_indices(i, j, 0);
#endif
  return static_cast<Vector2**>(m_array)[j][i];
}

inline StarStencil<Vector2> IceModelVec2V::star(int i, int j) const {
  const IceModelVec2V &self = *this;

  StarStencil<Vector2> result;

  result.ij = self(i,j);
  result.e =  self(i+1,j);
  result.w =  self(i-1,j);
  result.n =  self(i,j+1);
  result.s =  self(i,j-1);

  return result;
}

inline double& IceModelVec3D::operator() (int i, int j, int k) {
#if (PISM_DEBUG==1)
  check_array_indices(i, j, k);
#endif
  return static_cast<double***>(m_array)[j][i][k];
}

inline const double& IceModelVec3D::operator() (int i, int j, int k) const {
#if (PISM_DEBUG==1)
  check_array_indices(i, j, k);
#endif
  return static_cast<double***>(m_array)[j][i][k];
}

} // end of namespace pism

#endif /* _ICEMODELVEC_INLINE_H_ */
