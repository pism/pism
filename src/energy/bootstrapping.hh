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

#ifndef BOOTSTRAPPING_H
#define BOOTSTRAPPING_H

#include "pism/util/EnthalpyConverter.hh"

namespace pism {

namespace array {
class Array3D;
class Scalar;
}

namespace energy {

/*!
 * A heuristic formula for the temperature distribution within a column of ice. Used during
 * bootstrapping. A simple quartic guess.
 *
 * @param[in] EC enthalpy converter
 * @param[in] H ice thickness
 * @param[in] z height above the base of the ice, `0 <= z <= H`
 * @param[in] T_surface surface temperature, in kelvin
 * @param[in] G upward basal heat flux, in `W / meter^2`
 * @param[in] ice_k thermal conductivity of ice
 */
double ice_temperature_guess(EnthalpyConverter::Ptr EC,
                             double H, double z, double T_surface,
                             double G, double ice_k);

/*!
 * A heuristic formula for the temperature distribution within a column of ice. Used during
 * bootstrapping. Uses the surface mass balance to set the vertical velocity at the top surface.
 *
 * @param[in] EC enthalpy converter
 * @param[in] H ice thickness
 * @param[in] z height above the base of the ice, `0 <= z <= H`
 * @param[in] T_surface surface temperature, in kelvin
 * @param[in] G upward basal heat flux, in `W / meter^2`
 * @param[in] ice_k thermal conductivity of ice
 * @param[in] K temperature diffusivity in ice, ice_k / (ice_density * ice_c), where ice_c is the
 *              ice specific heat capacity
 * @param[in] SMB surface mass balance in `m / second`
 */
double ice_temperature_guess_smb(EnthalpyConverter::Ptr EC,
                                 double H, double z, double T_surface,
                                 double G, double ice_k, double K, double SMB);

} // end of namespace energy
} // end of namespace pism

#endif /* BOOTSTRAPPING_H */
