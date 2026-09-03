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

#ifndef PISM_MPDATA_3D_H
#define PISM_MPDATA_3D_H

#include <vector>

#include "pism/geometry/TransportScheme.hh"
#include "pism/util/array/Array3D.hh"
#include "pism/util/array/Scalar.hh"

namespace pism {

/*!
 * MPDATA (Smolarkiewicz, 1983; Smolarkiewicz and Grabowski, 1990) in three dimensions for
 * a quantity stored as mass per unit area in the control volumes around the levels of
 * PISM's vertical grid.
 *
 * Horizontal fluxes transport mass per unit area (the grid is uniform in x and y);
 * vertical fluxes transport the concentration (mass per unit area divided by the volume
 * thickness), so the non-uniform vertical spacing and the truncated top volume are
 * handled exactly by the finite-volume formulation. The first pass is first-order
 * upwinding with a positivity limiter for thin volumes; each following pass advects the
 * result with the anti-diffusive velocity of the previous pass. With `nonoscillatory` set
 * the flux-corrected-transport limiter keeps the solution within local bounds.
 */
class MPDATA3 : public TransportScheme3D {
public:
  MPDATA3(std::shared_ptr<const Grid> grid, int N, bool nonoscillatory);

  void update(double dt, const array::Scalar1 &ice_thickness, const array::CellType1 &cell_type,
              const array::Array3D &x, const array::Array3D &u, const array::Array3D &v,
              const array::Array3D &w);

  const array::Array3D &x() const;

private:
  void compute_interface_velocity(const array::Scalar1 &ice_thickness,
                                  const array::CellType1 &cell_type,
                                  const array::Array3D &u, const array::Array3D &v,
                                  const array::Array3D &w);

  void compute_corrective_velocity(double dt, const array::Scalar1 &ice_thickness);

  void limit(double dt, const array::Scalar1 &ice_thickness);

  void step(double dt, const array::Scalar1 &ice_thickness);

  //! thickness of volume `k` in a column of height `H`
  double dz(int k, double H) const;

  //! interfaces of the control volumes
  std::vector<double> m_zi;
  //! distances between the centers of adjacent volumes (size Mz - 1)
  std::vector<double> m_dz_face;

  //! face velocities: x-faces (east of a cell), y-faces (north of a cell), z-faces (top of
  //! a volume); index k of the z-face array is the face between volumes k and k + 1
  array::Array3D m_u_face, m_v_face, m_w_face;
  //! face velocities of the previous pass
  array::Array3D m_u_old, m_v_old, m_w_old;

  //! state at the start of a pass (ghosted, width 2)
  array::Array3D m_x_previous;
  //! input of the whole update (ghosted; needed by the limiter)
  array::Array3D m_x_input;
  //! result
  array::Array3D m_x;

  int m_N;
  bool m_nonoscillatory;
};

} // end of namespace pism

#endif /* PISM_MPDATA_3D_H */
