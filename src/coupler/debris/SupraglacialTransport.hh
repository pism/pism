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

#ifndef PISM_DEBRIS_SUPRAGLACIAL_TRANSPORT_HH
#define PISM_DEBRIS_SUPRAGLACIAL_TRANSPORT_HH

#include <memory>

#include "pism/util/array/Scalar.hh"
#include "pism/util/array/Vector.hh"

namespace pism {

class Grid;
class TransportScheme2D;

namespace array {
class Array3D;
class CellType1;
}

namespace debris {

//! @brief Advection of the supraglacial debris layer by the ice surface velocity (the
//! second term on the right-hand side of equation 17 in Verhaegen and Huybrechts, 2026).
class SupraglacialTransport {
public:
  SupraglacialTransport(std::shared_ptr<const Grid> grid, std::shared_ptr<TransportScheme2D> scheme);

  //! Extract the surface velocity from the 3D velocity field.
  void update_surface_velocity(const array::Array3D &u, const array::Array3D &v,
                               const array::Scalar &ice_thickness);

  //! Advance the debris thickness `h_d` by `dt` using the current surface velocity.
  void step(double dt, const array::CellType1 &cell_type, array::Scalar &h_d);

  const array::Vector &surface_velocity() const;

private:
  std::shared_ptr<const Grid> m_grid;
  std::shared_ptr<TransportScheme2D> m_scheme;
  array::Vector m_surface_velocity;
  array::Scalar m_u_surface, m_v_surface;
};

} // end of namespace debris
} // end of namespace pism

#endif /* PISM_DEBRIS_SUPRAGLACIAL_TRANSPORT_HH */
