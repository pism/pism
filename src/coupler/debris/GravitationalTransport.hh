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

#ifndef PISM_DEBRIS_GRAVITATIONAL_TRANSPORT_HH
#define PISM_DEBRIS_GRAVITATIONAL_TRANSPORT_HH

#include <memory>

#include "pism/util/array/Scalar.hh"
#include "pism/util/array/Staggered.hh"

namespace pism {

class Grid;

namespace array {
class CellType1;
}

namespace debris {

//! @brief Gravitational (slope-driven) redistribution of supraglacial debris (equations
//! 21 and 22 in Verhaegen and Huybrechts, 2026).
/*!
  The debris flux is `F = -K grad(h_s + h_d)` with `K = mu (1 - phi) rho g h_d`, i.e. a
  nonlinear diffusion of the debris thickness down the slope of the debris surface. The
  flux is computed on the staggered grid with `K` taken from the upwind cell (compare
  `hydrology::Routing`), limited so that `h_d` stays non-negative, and the update is
  sub-cycled to satisfy the explicit stability conditions (diffusive and advective).
*/
class GravitationalTransport {
public:
  /*!
   * @param[in] mobility `mu` (Pa^-1 m^2 s^-1)
   * @param[in] solid_density `(1 - phi) rho` (kg m^-3)
   * @param[in] cfl_ratio fraction of the explicit stability limit used for sub-steps
   */
  GravitationalTransport(std::shared_ptr<const Grid> grid, double mobility, double solid_density,
                         double cfl_ratio);

  //! Advance `h_d` by `dt` (sub-cycling internally).
  void step(double dt, const array::CellType1 &cell_type, const array::Scalar &surface_elevation,
            array::Scalar1 &h_d);

  //! Debris flux on the staggered grid during the last sub-step (m^2 s^-1).
  const array::Staggered1 &flux() const;

  //! Number of sub-steps taken during the last call to step().
  int substeps() const;

  //! Largest stable time step for the diffusion coefficient `K_max`.
  double max_dt(double K_max) const;

  //! `K / h_d = mu (1 - phi) rho g` (m s^-1).
  double coefficient() const;

private:
  double compute_flux(const array::CellType1 &cell_type, const array::Scalar1 &h_d);
  void update_surface(const array::Scalar &surface_elevation, const array::Scalar1 &h_d);

  std::shared_ptr<const Grid> m_grid;
  array::Staggered1 m_flux, m_flux_limited;
  array::Scalar1 m_surface;
  double m_coefficient, m_cfl_ratio;
  int m_substeps;
};

} // end of namespace debris
} // end of namespace pism

#endif /* PISM_DEBRIS_GRAVITATIONAL_TRANSPORT_HH */
