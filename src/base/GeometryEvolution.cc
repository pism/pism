/* Copyright (C) 2016 PISM Authors
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

namespace pism {

struct Impl {
  IceModelVec2S flux_divergence;
  IceModelVec2S conservation_error;
  IceModelVec2S effective_SMB;
  IceModelVec2S effective_BMB;
  IceModelVec2S thickness_change;
};

void GeometryEvolution::step(Geometry &ice_geometry, double dt,
                             const IceModelVec2V    &advective_velocity,
                             const IceModelVec2Stag &diffusive_flux,
                             const IceModelVec2Int  &velocity_bc_mask,
                             const IceModelVec2V    &advective_bc_values,
                             const IceModelVec2Int  &thickness_bc_mask) {

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

  update_cell_type();

  // uses thickness_bc_mask
  apply_surface_and_basal_mass_balance();

  ice_geometry.ensure_consistency();

  // calving is a separate issue
}

} // end of namespace pism
