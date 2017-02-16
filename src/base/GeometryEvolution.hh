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

#ifndef GEOMETRYEVOLUTION_H
#define GEOMETRYEVOLUTION_H

#include "Geometry.hh"

namespace pism {

/*!
 * NB! Write this in a way that does not use ghosts of input fields (copy to temp. storage and
 * communicate).
 *
 * The promise:
 *
 * H_change + Href_change = dt * (SMB_rate + BMB_rate - flux_divergence) + conservation_error
 *
 * Href == 0 if H > 0
 */
class GeometryEvolution {
public:
  GeometryEvolution(IceGrid::ConstPtr grid);

  void step(Geometry &ice_geometry, double dt,
            const IceModelVec2V    &advective_velocity,
            const IceModelVec2Stag &diffusive_flux,
            const IceModelVec2Int  &velocity_bc_mask,
            const IceModelVec2V    &velocity_bc_values,
            const IceModelVec2Int  &thickness_bc_mask,
            const IceModelVec2S    &surface_mass_balance_rate,
            const IceModelVec2S    &basal_mass_balance_rate);

  const IceModelVec2S& flux_divergence() const;
  const IceModelVec2S& top_surface_mass_balance() const;
  const IceModelVec2S& bottom_surface_mass_balance() const;
  const IceModelVec2S& thickness_change() const;
  const IceModelVec2S& conservation_error() const;
protected:
  void compute_interface_velocity(const IceModelVec2CellType &cell_type,
                                  const IceModelVec2V &advective_velocity,
                                  const IceModelVec2Int &velocity_bc_mask,
                                  const IceModelVec2V &velocity_bc_values,
                                  IceModelVec2Stag &velocity_staggered);

  void compute_interface_fluxes(const IceModelVec2CellType &cell_type,
                                const IceModelVec2Stag &velocity_staggered,
                                const IceModelVec2S &ice_thickness,
                                const IceModelVec2Stag &diffusive_flux,
                                IceModelVec2Stag &flux_staggered);

  void compute_flux_divergence(const IceModelVec2Stag &flux_staggered,
                               const IceModelVec2Int &thickness_bc_mask,
                               IceModelVec2S &flux_fivergence);

  void apply_flux_divergence(const IceModelVec2CellType &cell_type,
                             const IceModelVec2S &ice_thickness,
                             const IceModelVec2S &ice_surface_elevation,
                             const IceModelVec2S &bed_elevation,
                             IceModelVec2S &thickness_change,
                             IceModelVec2S &area_specific_volume_change);

  void ensure_nonnegativity(const IceModelVec2S &ice_thickness,
                            const IceModelVec2S &area_specific_volume,
                            IceModelVec2S &thickness_change,
                            IceModelVec2S &area_specific_volume_change,
                            IceModelVec2S &conservation_error);

  // note: cells with area_specific_volume > 0 do not experience changes due to surface and basal
  // mass balance sources
  void apply_surface_and_basal_mass_balance(const IceModelVec2S &ice_thickness,
                                            const IceModelVec2Int &thickness_bc_mask,
                                            const IceModelVec2S &surface_mass_balance,
                                            const IceModelVec2S &basal_mass_balance,
                                            IceModelVec2S &thickness_change);
private:
  struct Impl;
  Impl *m_impl;
};

} // end of namespace pism

#endif /* GEOMETRYEVOLUTION_H */
