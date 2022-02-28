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

namespace array {
class Array3D;
class Scalar;
}

namespace energy {

void compute_temperature(const array::Array3D &enthalpy,
                         const array::Scalar &ice_thickness,
                         array::Array3D &result);

void compute_enthalpy(const array::Array3D &temperature,
                      const array::Array3D &liquid_water_fraction,
                      const array::Scalar &ice_thickness,
                      array::Array3D &result);

void compute_enthalpy_cold(const array::Array3D &temperature,
                           const array::Scalar &ice_thickness,
                           array::Array3D &result);

void compute_liquid_water_fraction(const array::Array3D &enthalpy,
                                   const array::Scalar &ice_thickness,
                                   array::Array3D &result);

void compute_cts(const array::Array3D &enthalpy,
                 const array::Scalar &ice_thickness,
                 array::Array3D &result);

double total_ice_enthalpy(double thickness_threshold,
                          const array::Array3D &ice_enthalpy,
                          const array::Scalar &ice_thickness);

void bootstrap_ice_temperature(const array::Scalar &ice_thickness,
                               const array::Scalar &ice_surface_temp,
                               const array::Scalar &surface_mass_balance,
                               const array::Scalar &basal_heat_flux,
                               array::Array3D &result);

void bootstrap_ice_enthalpy(const array::Scalar &ice_thickness,
                            const array::Scalar &ice_surface_temp,
                            const array::Scalar &surface_mass_balance,
                            const array::Scalar &basal_heat_flux,
                            array::Array3D &result);

} // end of namespace energy
} // end of namespace pism

#endif /* UTILITIES_H */
