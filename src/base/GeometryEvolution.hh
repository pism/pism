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

#ifndef GEOMETRYEVOLUTION_H
#define GEOMETRYEVOLUTION_H

#include "IceGeometry.hh"

namespace pism {

/*!
 * NB! Write this in a way that does not use ghosts of input fields (copy to temp. storage and
 * communicate).
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
private:
  struct Impl;
  Impl *m_impl;
};

} // end of namespace pism

#endif /* GEOMETRYEVOLUTION_H */
