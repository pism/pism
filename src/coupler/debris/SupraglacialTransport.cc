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

#include "pism/coupler/debris/SupraglacialTransport.hh"

#include "pism/geometry/TransportScheme.hh"
#include "pism/util/Grid.hh"
#include "pism/util/array/Array3D.hh"
#include "pism/util/array/CellType.hh"

namespace pism {
namespace debris {

SupraglacialTransport::SupraglacialTransport(std::shared_ptr<const Grid> grid,
                                             std::shared_ptr<TransportScheme2D> scheme)
  : m_grid(grid),
    m_scheme(scheme),
    m_surface_velocity(grid, "debris_surface_velocity"),
    m_u_surface(grid, "u_surface"),
    m_v_surface(grid, "v_surface") {

  m_surface_velocity.metadata(0)
      .long_name("x-component of the ice surface velocity advecting supraglacial debris")
      .units("m s^-1")
      .output_units("m year^-1");
  m_surface_velocity.metadata(1)
      .long_name("y-component of the ice surface velocity advecting supraglacial debris")
      .units("m s^-1")
      .output_units("m year^-1");
  m_surface_velocity.set(0.0);
}

void SupraglacialTransport::update_surface_velocity(const array::Array3D &u,
                                                    const array::Array3D &v,
                                                    const array::Scalar &ice_thickness) {
  array::extract_surface(u, ice_thickness, m_u_surface);
  array::extract_surface(v, ice_thickness, m_v_surface);

  array::AccessScope list{ &m_u_surface, &m_v_surface, &m_surface_velocity };

  for (auto p : m_grid->points()) {
    const int i = p.i(), j = p.j();
    m_surface_velocity(i, j) = { m_u_surface(i, j), m_v_surface(i, j) };
  }
}

void SupraglacialTransport::step(double dt, const array::CellType1 &cell_type, array::Scalar &h_d) {
  m_scheme->update(dt, cell_type, h_d, m_surface_velocity);
  h_d.copy_from(m_scheme->x());
}

const array::Vector &SupraglacialTransport::surface_velocity() const {
  return m_surface_velocity;
}

} // end of namespace debris
} // end of namespace pism
