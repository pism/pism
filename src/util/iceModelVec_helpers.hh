// Copyright (C) 2011, 2013, 2014, 2016, 2017 PISM Authors
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

#ifndef _ICEMODELVEC_HELPERS_H_
#define _ICEMODELVEC_HELPERS_H_

namespace pism {

void compute_params(const IceModelVec* const x, const IceModelVec* const y,
                    const IceModelVec* const z, int &ghosts, bool &scatter);

//! \brief Computes result = x + alpha * y, where x, y, and z are 2D
//! IceModelVecs (scalar or vector).
/*!
 * This implementation tries to be smart about handling IceModelVecs with and
 * without ghosts and with different stencil widths.
 *
 * This template function was written to re-use this code for both
 * IceModelVec2S and IceModel2V.
 *
 * This cannot go into a protected member IceModelVec because
 * IceModelVec2S::operator() and IceModelVec2V::operator() return different
 * types.
 *
 * Note: this code uses overloaded operators (Vector2::operator*, etc).
 */
template<class V>
void add_2d(const IceModelVec* const x_in, double alpha,
            const IceModelVec* const y_in,
            IceModelVec* const result) {
  const V *x = dynamic_cast<const V*>(x_in),
    *y = dynamic_cast<const V*>(y_in);

  V *z = dynamic_cast<V*>(result);

  if (x == NULL || y == NULL || z == NULL) {
    throw RuntimeError(PISM_ERROR_LOCATION, "incompatible arguments");
  }

  int stencil = 0;
  bool scatter = false;
  compute_params(x, y, z, stencil, scatter);

  IceModelVec::AccessList list{x, y, z};
  for (PointsWithGhosts p(*z->grid(), stencil); p; p.next()) {
    const int i = p.i(), j = p.j();

    (*z)(i, j) = (*x)(i, j) + (*y)(i, j) * alpha;
  }

  if (scatter) {
    z->update_ghosts();
  }

  result->inc_state_counter();
}

template<class V>
void copy_2d(const IceModelVec* const source,
             IceModelVec* const destination) {
  const V *x = dynamic_cast<const V*>(source);

  V *z = dynamic_cast<V*>(destination);

  if (x == NULL || z == NULL) {
    throw RuntimeError(PISM_ERROR_LOCATION, "incompatible arguments");
  }

  int stencil = 0;
  bool scatter = false;
  compute_params(x, x, z, stencil, scatter);

  IceModelVec::AccessList list{x, z};

  for (PointsWithGhosts p(*z->grid(), stencil); p; p.next()) {
    const int i = p.i(), j = p.j();

    (*z)(i, j) = (*x)(i, j);
  }

  if (scatter) {
    z->update_ghosts();
  }

  destination->inc_state_counter();
}

} // end of namespace pism

#endif /* _ICEMODELVEC_HELPERS_H_ */
