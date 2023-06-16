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

#ifndef PISM_UNO_2D_H
#define PISM_UNO_2D_H

#include "pism/util/array/Scalar.hh"
#include "pism/util/array/Vector.hh"
#include "pism/util/array/Staggered.hh"

namespace pism {

class IceGrid;

namespace array {
class CellType1;
} // end of namespace array

enum UNOType {PISM_UNO_UPWIND1, PISM_UNO_LAX_WENDROFF, PISM_UNO_FROMM, PISM_UNO_2, PISM_UNO_3};

/*!
 * Implementation of Upwind Nonoscillatory advection schemes UNO2 and UNO3.
 *
 * See J.-G. Li, “Upstream Nonoscillatory Advection Schemes,” Monthly Weather Review, vol.
 * 136, Art. no. 12, Dec. 2008.
 *
 * Three more schemes (first-order upwinding, Lax-Wendroff and Fromm) can be stated in the
 * same framework and so we implement them as well since this requires very little effort
 * and simplifies testing and comparisons.
 */
class UNO {
public:
  UNO(std::shared_ptr<const IceGrid> grid, UNOType type);

  void update(double dt,
              const array::CellType1 &cell_type,
              const array::Scalar &x,
              const array::Vector &velocity,
              bool nonnegative = false);

  const array::Scalar& x() const;

private:
  void compute_interface_fluxes(const array::CellType1 &cell_type,
                                const array::Vector1 &velocity,
                                const array::Scalar2 &x_old,
                                double dt,
                                array::Staggered &result) const;

  typedef double (*MidFluxApproximation)(const double *, const double *, size_t, double, double, double);

  MidFluxApproximation m_approx;

  // interface fluxes (on the staggered grid; ghosted)
  array::Staggered1 m_q, m_q_limited;

  // ghosted copy of the input velocity field used to compute interface velocities
  array::Vector1 m_v_ghosted;

  // temporary storage for the old state of the advected quantity (ghosted, width=2)
  //
  // used to perform the step using flux divergence *and* to compute interface fluxes
  // (this requires a wide stencil)
  array::Scalar2 m_x_ghosted;

  // resulting advected quantity (not ghosted)
  array::Scalar m_x;
};

} // end of namespace pism

#endif /* PISM_UNO_2D_H */
