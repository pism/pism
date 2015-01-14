/* Copyright (C) 2015 PISM Authors
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
  return static_cast<double***>(array)[i][j][k];
}

inline const double& IceModelVec2::operator() (int i, int j, int k) const {
#if (PISM_DEBUG==1)
  check_array_indices(i, j, k);
#endif
  return static_cast<double***>(array)[i][j][k];
}

inline double& IceModelVec2S::operator() (int i, int j) {
#if (PISM_DEBUG==1)
  check_array_indices(i, j, 0);
#endif
  return static_cast<double**>(array)[i][j];
}

inline const double& IceModelVec2S::operator()(int i, int j) const {
#if (PISM_DEBUG==1)
  check_array_indices(i, j, 0);
#endif
  return static_cast<double**>(array)[i][j];
}

inline planeStar<double> IceModelVec2S::star(int i, int j) {
#if (PISM_DEBUG==1)
  check_array_indices(i, j, 0);
  check_array_indices(i+1, j, 0);
  check_array_indices(i-1, j, 0);
  check_array_indices(i, j+1, 0);
  check_array_indices(i, j-1, 0);
#endif
  planeStar<double> result;

  result.ij = operator()(i,j);
  result.e =  operator()(i+1,j);
  result.w =  operator()(i-1,j);
  result.n =  operator()(i,j+1);
  result.s =  operator()(i,j-1);

  return result;
}

inline planeStar<double> IceModelVec2Stag::star(int i, int j) const {
#if (PISM_DEBUG==1)
  check_array_indices(i, j, 0);
  check_array_indices(i+1, j, 0);
  check_array_indices(i-1, j, 0);
  check_array_indices(i, j+1, 0);
  check_array_indices(i, j-1, 0);
#endif
  planeStar<double> result;

  result.ij = 0.0;             // has no meaning in this context
  result.e =  operator()(i, j, 0);
  result.w =  operator()(i-1, j, 0);
  result.n =  operator()(i, j, 1);
  result.s =  operator()(i, j-1, 1);

  return result;
}

inline int IceModelVec2Int::as_int(int i, int j) const {
#if (PISM_DEBUG==1)
  check_array_indices(i, j, 0);
#endif
  const double **a = (const double**) array;
  return static_cast<int>(floor(a[i][j] + 0.5));
}

inline planeStar<int> IceModelVec2Int::int_star(int i, int j) const {
#if (PISM_DEBUG==1)
  check_array_indices(i, j, 0);
  check_array_indices(i+1, j, 0);
  check_array_indices(i-1, j, 0);
  check_array_indices(i, j+1, 0);
  check_array_indices(i, j-1, 0);
#endif

  planeStar<int> result;
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
  return static_cast<Vector2**>(array)[i][j];
}

inline const Vector2& IceModelVec2V::operator()(int i, int j) const {
#if (PISM_DEBUG==1)
  check_array_indices(i, j, 0);
#endif
  return static_cast<Vector2**>(array)[i][j];
}

inline planeStar<Vector2> IceModelVec2V::star(int i, int j) const {
#if (PISM_DEBUG==1)
  check_array_indices(i, j, 0);
  check_array_indices(i+1, j, 0);
  check_array_indices(i-1, j, 0);
  check_array_indices(i, j+1, 0);
  check_array_indices(i, j-1, 0);
#endif
  planeStar<Vector2> result;

  result.ij = operator()(i,j);
  result.e =  operator()(i+1,j);
  result.w =  operator()(i-1,j);
  result.n =  operator()(i,j+1);
  result.s =  operator()(i,j-1);

  return result;
}

inline double& IceModelVec3D::operator() (int i, int j, int k) {
#if (PISM_DEBUG==1)
  check_array_indices(i, j, k);
#endif
  return static_cast<double***>(array)[i][j][k];
}

inline const double& IceModelVec3D::operator() (int i, int j, int k) const {
#if (PISM_DEBUG==1)
  check_array_indices(i, j, k);
#endif
  return static_cast<double***>(array)[i][j][k];
}

} // end of namespace pism

#endif /* _ICEMODELVEC_INLINE_H_ */
