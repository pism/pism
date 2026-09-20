// Copyright (C) 2026 Andy Aschwanden
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 2 of the License, or (at your option) any later
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

#ifndef PISM_OCEAN_PLUME_H
#define PISM_OCEAN_PLUME_H

#include "pism/coupler/ocean/Picop.hh"

namespace pism {

namespace array {
class Forcing;
}

namespace ocean {

//! PICOP's buoyant-plume sub-shelf melt without the PICO box model.
/*!
 * The ambient ocean temperature `T_a` and salinity `S_a` driving the plume are read from
 * a file (`theta_ocean`, `salinity_ocean`) at every floating cell instead of being
 * averaged over an ocean basin and passed through PICO's boxes. Use it when the forcing
 * already resolves the near-glacier ocean state, e.g. thermal forcing extrapolated into
 * the fjords as in the ISMIP6/ISMIP7 Greenland protocol.
 *
 * With `ocean.plume.temperature_as_thermal_forcing` set, `theta_ocean` is thermal
 * forcing and `T_a = T_f(S_a, z_gl) + TF`, where `T_f` is the freezing point at the
 * grounding-line depth (Pelle et al. 2019, eqn. 4), so the plume sees exactly the
 * thermal forcing in the file.
 */
class Plume : public PlumeModel {
public:
  Plume(std::shared_ptr<const Grid> g);
  virtual ~Plume() = default;

protected:
  void init_impl(const Geometry &geometry);
  void update_impl(const Inputs &inputs, double t, double dt);
  MaxTimestep max_timestep_impl(double t, const CFLData *cfl_data) const;

  DiagnosticList spatial_diagnostics_impl() const;

private:
  std::shared_ptr<array::Forcing> m_theta_ocean, m_salinity_ocean;

  //! ambient temperature T_a (kelvin) and salinity S_a (g/kg) driving the plume on
  //! floating cells; named after PICOP's diagnostics so they are written under the
  //! same names
  array::Scalar m_ambient_temperature, m_ambient_salinity;

  //! whether theta_ocean holds thermal forcing rather than potential temperature
  bool m_thermal_forcing;

  void compute_ambient_temperature(const PicopPhysics &physics, const Geometry &geometry,
                                   array::Scalar &result) const;

  void compute_shelf_base_temperature(const PicopPhysics &physics, const Geometry &geometry,
                                      array::Scalar &result) const;
};

} // end of namespace ocean
} // end of namespace pism

#endif /* PISM_OCEAN_PLUME_H */
