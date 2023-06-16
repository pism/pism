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

#ifndef PISM_MPDATA_2D_H
#define PISM_MPDATA_2D_H

#include "pism/util/array/Scalar.hh"
#include "pism/util/array/Vector.hh"
#include "pism/util/array/Staggered.hh"

namespace pism {

class IceGrid;

namespace array {
class CellType;
} // end of namespace array

class MPDATA2 {
public:
  MPDATA2(std::shared_ptr<const IceGrid> grid, int N);

  void update(double dt,
              const array::CellType &cell_type,
              const array::Scalar &x,
              const array::Vector &velocity,
              bool nonoscillatory = false);

  const array::Scalar& x() const;

private:
  // interface velocity (on the staggered grid; ghosted)
  array::Staggered1 m_v, m_v_old;

  array::Vector1 m_v_ghosted;

  // temporary storage for the result of the previous iteration
  array::Scalar2 m_x_previous;
  array::Scalar2 m_x_input;

  // advected quantity
  array::Scalar m_x;

  // number of iterations
  int m_N;
};

} // end of namespace pism

#endif /* PISM_MPDATA_2D_H */
