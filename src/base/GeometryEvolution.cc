/* Copyright (C) 2016, 2017 PISM Authors
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

#include "GeometryEvolution.hh"

#include "base/util/iceModelVec.hh"

namespace pism {

struct GeometryEvolution::Impl {
  IceGrid::ConstPtr grid;

  /*! @brief Fluxes through cell interfaces (sides). */
  IceModelVec2Stag interface_fluxes;
  /*! @brief Flux divergence (used to track thickness changes due to flow). */
  IceModelVec2S flux_divergence;
  /*! @brief Conservation error due to enforcing non-negativity of ice thickness. */
  IceModelVec2S conservation_error;
  /*! @brief Effective surface mass balance. */
  IceModelVec2S effective_SMB;
  /*! @brief Effective basal mass balance. */
  IceModelVec2S effective_BMB;
  /*! @brief Change in ice thickness during the last time step. */
  IceModelVec2S thickness_change;
};

GeometryEvolution::GeometryEvolution(IceGrid::ConstPtr grid) {
  m_impl->grid = grid;
}

const IceModelVec2S& GeometryEvolution::flux_divergence() const {
  return m_impl->flux_divergence;
}

const IceModelVec2S& GeometryEvolution::top_surface_mass_balance() const {
  return m_impl->effective_SMB;
}

const IceModelVec2S& GeometryEvolution::bottom_surface_mass_balance() const {
  return m_impl->effective_BMB;
}

const IceModelVec2S& GeometryEvolution::thickness_change() const {
  return m_impl->thickness_change;
}

const IceModelVec2S& GeometryEvolution::conservation_error() const {
  return m_impl->conservation_error;
}

void GeometryEvolution::step(Geometry &ice_geometry, double dt,
                             const IceModelVec2V    &advective_velocity,
                             const IceModelVec2Stag &diffusive_flux,
                             const IceModelVec2Int  &velocity_bc_mask,
                             const IceModelVec2V    &velocity_bc_values,
                             const IceModelVec2Int  &thickness_bc_mask,
                             const IceModelVec2S    &surface_mass_balance_rate,
                             const IceModelVec2S    &basal_mass_balance_rate) {

  // uses B.C. data and includes modifications for regional runs
  compute_interface_velocity();

  // upwinding shows up here; includes modifications for regional runs
  compute_interface_fluxes();

  // simple finite differences
  compute_flux_divergence();

  // this is where part_grid should be implemented; uses thickness_bc_mask
  apply_flux_divergence();

  // computes the numerical conservation error
  ensure_nonnegativity();

  // ice extent may change due to the flow, so we need to update the mask in preparation for
  // applying the source terms (surface and basal mass balance)
  update_cell_type();

  // uses thickness_bc_mask
  apply_surface_and_basal_mass_balance();

  ice_geometry.ensure_consistency(m_impl->grid->ctx()->config()->get_double("geometry.ice_free_thickness_standard"));

  // calving is a separate issue
}

void GeometryEvolution::compute_interface_velocity() {

}

void GeometryEvolution::compute_interface_fluxes() {

}

void GeometryEvolution::compute_flux_divergence() {

}

void GeometryEvolution::apply_flux_divergence() {

}

void GeometryEvolution::ensure_nonnegativity() {

}

void GeometryEvolution::update_cell_type() {

}

void GeometryEvolution::apply_surface_and_basal_mass_balance() {

}


} // end of namespace pism
