/* Copyright (C) 2026 PISM Authors
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

#ifndef PISM_TRANSPORT_SCHEME_H
#define PISM_TRANSPORT_SCHEME_H

#include <memory>
#include <string>

namespace pism {

class Grid;

namespace array {
class Array3D;
class CellType1;
class Scalar;
class Scalar1;
class Vector;
} // end of namespace array

/*!
 * Explicit, conservative advection of a map-plane scalar by a cell-centered velocity
 * field. Implementations wrap the schemes in `MPDATA2.hh` and `UNO.hh`.
 *
 * Faces between ice-covered and ice-free cells are treated as in `GeometryEvolution`: the
 * interface velocity is the velocity of the ice-covered side, so mass can leave the ice.
 */
class TransportScheme2D {
public:
  virtual ~TransportScheme2D() = default;

  //! Advance `x` by `dt`; the result is available through x().
  virtual void update(double dt, const array::CellType1 &cell_type, const array::Scalar &x,
                      const array::Vector &velocity) = 0;

  virtual const array::Scalar &x() const = 0;

  /*!
   * Create a scheme.
   *
   * @param[in] kind one of "upwind", "mpdata", "uno2", "uno3"
   * @param[in] N number of MPDATA passes (1 is first-order upwinding)
   * @param[in] nonoscillatory use the flux-corrected-transport limiter in MPDATA
   */
  static std::shared_ptr<TransportScheme2D> create(std::shared_ptr<const Grid> grid,
                                                   const std::string &kind, int N,
                                                   bool nonoscillatory);
};

/*!
 * Explicit, conservative advection of a quantity stored as *mass per unit area* in
 * control volumes around the levels of PISM's vertical grid (see
 * `debris::column::interfaces()`) by a 3D velocity field.
 *
 * All faces of the ice body are closed: horizontal faces are open only where both
 * columns contain ice at that level, the bed is closed, and so is the truncated top
 * volume. The total mass is therefore conserved to rounding error; ice removed at the
 * surface or the bed has to be handled by the caller.
 */
class TransportScheme3D {
public:
  virtual ~TransportScheme3D() = default;

  /*!
   * Advance `x` by `dt`.
   *
   * @param[in] ice_thickness column heights (ghosted)
   * @param[in] cell_type cell type mask (ghosted)
   * @param[in] x mass per unit area per level (no ghosts needed)
   * @param[in] u,v horizontal velocity components on the vertical grid (ghosted)
   * @param[in] w vertical velocity relative to the bed
   */
  virtual void update(double dt, const array::Scalar1 &ice_thickness,
                      const array::CellType1 &cell_type, const array::Array3D &x,
                      const array::Array3D &u, const array::Array3D &v,
                      const array::Array3D &w) = 0;

  virtual const array::Array3D &x() const = 0;

  //! @param[in] kind one of "upwind", "mpdata"
  static std::shared_ptr<TransportScheme3D> create(std::shared_ptr<const Grid> grid,
                                                   const std::string &kind, int N,
                                                   bool nonoscillatory);
};

} // end of namespace pism

#endif /* PISM_TRANSPORT_SCHEME_H */
