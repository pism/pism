/* Copyright (C) 2016, 2018 PISM Authors
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

#ifndef UTILITIES_H
#define UTILITIES_H

namespace pism {

class IceModelVec2S;
class IceModelVec3;

namespace energy {

void compute_temperature(const IceModelVec3 &enthalpy,
                         const IceModelVec2S &ice_thickness,
                         IceModelVec3 &result);

void compute_enthalpy(const IceModelVec3 &temperature,
                      const IceModelVec3 &liquid_water_fraction,
                      const IceModelVec2S &ice_thickness,
                      IceModelVec3 &result);

void compute_enthalpy_cold(const IceModelVec3 &temperature,
                           const IceModelVec2S &ice_thickness,
                           IceModelVec3 &result);

void compute_liquid_water_fraction(const IceModelVec3 &enthalpy,
                                   const IceModelVec2S &ice_thickness,
                                   IceModelVec3 &result);

void compute_cts(const IceModelVec3 &enthalpy,
                 const IceModelVec2S &ice_thickness,
                 IceModelVec3 &result);

double total_ice_enthalpy(double thickness_threshold,
                          const IceModelVec3 &ice_enthalpy,
                          const IceModelVec2S &ice_thickness);

void bootstrap_ice_temperature(const IceModelVec2S &ice_thickness,
                               const IceModelVec2S &ice_surface_temp,
                               const IceModelVec2S &surface_mass_balance,
                               const IceModelVec2S &basal_heat_flux,
                               IceModelVec3 &result);

void bootstrap_ice_enthalpy(const IceModelVec2S &ice_thickness,
                            const IceModelVec2S &ice_surface_temp,
                            const IceModelVec2S &surface_mass_balance,
                            const IceModelVec2S &basal_heat_flux,
                            IceModelVec3 &result);

} // end of namespace energy
} // end of namespace pism

#endif /* UTILITIES_H */
