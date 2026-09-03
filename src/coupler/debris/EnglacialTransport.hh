// Copyright (C) 2026 Andy Aschwanden and Constantine Khroulev
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

#ifndef PISM_DEBRIS_ENGLACIAL_TRANSPORT_HH
#define PISM_DEBRIS_ENGLACIAL_TRANSPORT_HH

#include <memory>
#include <vector>

#include "pism/util/array/Array3D.hh"
#include "pism/util/array/Scalar.hh"

namespace pism {

class Grid;
class TransportScheme3D;

namespace array {
class CellType1;
}

namespace debris {

//! @brief Englacial debris: 3D advection, burial in the accumulation zone, melt-out in the
//! ablation zone (Section 2.1.3 of Verhaegen and Huybrechts, 2026).
/*!
  The state is the debris mass per unit area in the control volumes around the levels of
  PISM's vertical grid (see `column_kernels.hh`); the concentration (kg m^-3) is derived
  from it when needed.

  One step consists of

  1. advection with the 3D velocity field on the ice geometry the velocity was computed
     for (all faces of the ice body closed),
  2. truncation to the thickness after the flow step (mass above it is folded into the
     top volume, conserving mass),
  3. basal melt: the bottom of the column is removed and the rest shifts down,
  4. surface mass balance: ice removed at the top releases the debris it contains
     (`melt_out()`); ice added at the top is clean, except at debris input cells where it
     buries the prescribed input (`burial()`),
  5. columns that lost all their ice give up their debris (`ice_free_loss()`).
*/
class EnglacialTransport {
public:
  EnglacialTransport(std::shared_ptr<const Grid> grid, std::shared_ptr<TransportScheme3D> scheme,
                     double solid_density);

  /*!
   * Advance by `dt`.
   *
   * @param[in] H_old ice thickness the velocity field corresponds to (ghosted)
   * @param[in] cell_type_old cell type mask corresponding to `H_old` (ghosted)
   * @param[in] H_new ice thickness at the end of the step
   * @param[in] icy_new ice-covered cells at the end of the step
   * @param[in] top_smb ice thickness change at the surface during the step (m ice)
   * @param[in] bottom_smb ice thickness change at the base during the step (m ice)
   * @param[in] input_rate debris input rate (m of solid debris per second)
   */
  void step(double dt, const array::Scalar1 &H_old, const array::CellType1 &cell_type_old,
            const array::Scalar &H_new, const array::CellType1 &cell_type_new,
            const array::Scalar &top_smb, const array::Scalar &bottom_smb,
            const array::Scalar &input_rate,
            const array::Array3D &u, const array::Array3D &v, const array::Array3D &w);

  //! Debris mass per unit area per control volume (kg m^-2).
  array::Array3D &mass();
  const array::Array3D &mass() const;

  //! Debris released by surface melt during the last step (kg m^-2).
  const array::Scalar &melt_out() const;
  //! Debris buried by accumulation during the last step (kg m^-2).
  const array::Scalar &burial() const;
  //! Debris removed by basal melt during the last step (kg m^-2).
  const array::Scalar &basal_loss() const;
  //! Debris given up by columns that became ice-free during the last step (kg m^-2).
  const array::Scalar &ice_free_loss() const;

  //! Concentration (kg m^-3) corresponding to the current mass and the thickness `H`.
  void concentration(const array::Scalar &H, array::Array3D &result) const;
  //! Set the mass from a concentration field (kg m^-3) and the thickness `H`.
  void set_concentration(const array::Array3D &C, const array::Scalar &H);
  //! Mass per unit area integrated over the column (kg m^-2).
  void column_mass(array::Scalar &result) const;

  //! Total mass (kg) on this sub-domain (not reduced across processors).
  double local_mass() const;

private:
  std::shared_ptr<const Grid> m_grid;
  std::shared_ptr<TransportScheme3D> m_scheme;

  std::vector<double> m_zi;
  double m_solid_density;

  std::shared_ptr<array::Array3D> m_mass;
  array::Scalar m_melt_out, m_burial, m_basal_loss, m_ice_free_loss;
};

} // end of namespace debris
} // end of namespace pism

#endif /* PISM_DEBRIS_ENGLACIAL_TRANSPORT_HH */
