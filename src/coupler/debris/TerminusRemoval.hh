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

#ifndef PISM_DEBRIS_TERMINUS_REMOVAL_HH
#define PISM_DEBRIS_TERMINUS_REMOVAL_HH

#include <memory>

#include "pism/util/array/Scalar.hh"

namespace pism {

class Grid;

namespace array {
class CellType2;
}

namespace debris {

//! @brief Removal of supraglacial debris into the foreland at the glacier margin
//! (Section 2.1.5 and equation 23 in Verhaegen and Huybrechts, 2026).
/*!
  At ice-covered cells with a lower ice-free neighbor the gravitational debris flux
  `K_d(h) |grad(h_s + h_d)|` directed into the foreland, divided by the marginal length
  scale `Gamma`, is a debris thickness loss rate. If `Gamma` exceeds the grid spacing the
  debris thickness and the slope are averaged over the margin cell and
  `ceil(Gamma / dx) - 1` cells up-glacier along the steepest ascent (at most 2, limited by
  the width of the ghost zone).
*/
class TerminusRemoval {
public:
  /*!
   * @param[in] gamma marginal length scale (m)
   * @param[in] mobility `mu` (Pa^-1 m^2 s^-1)
   * @param[in] solid_density `(1 - phi) rho` (kg m^-3)
   */
  TerminusRemoval(std::shared_ptr<const Grid> grid, double gamma, double mobility,
                  double solid_density);

  //! Remove debris from margin cells during `dt`.
  void step(double dt, const array::CellType2 &cell_type, const array::Scalar &surface_elevation,
            array::Scalar2 &h_d);

  //! Debris thickness removal rate during the last step (m s^-1).
  const array::Scalar &removal_rate() const;

  //! Number of cells up-glacier of the margin included in the averages.
  int upstream_cells() const;

private:
  std::shared_ptr<const Grid> m_grid;
  array::Scalar2 m_surface;
  array::Scalar m_removal_rate;
  double m_gamma, m_coefficient;
  int m_n;
};

} // end of namespace debris
} // end of namespace pism

#endif /* PISM_DEBRIS_TERMINUS_REMOVAL_HH */
