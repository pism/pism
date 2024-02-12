// Copyright (C) 2009--2024 PISM Authors
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

#ifndef PISM_DEBM_SIMPLE_POINTWISE_H
#define PISM_DEBM_SIMPLE_POINTWISE_H

#include <memory>
#include <array>

#include "pism/util/Mask.hh"
#include "pism/util/ScalarForcing.hh"

namespace pism {

class Context;
class Time;

namespace surface {

// The following three could be nested in DEBMSimplePointwise but SWIG does not appear to
// support nested classes, so we put them here to make them accessible from Python.

struct DEBMSimpleMelt {
  DEBMSimpleMelt();

  double temperature_melt;
  double insolation_melt;
  double offset_melt;
  double total_melt;
};

struct DEBMSimpleChanges {
  DEBMSimpleChanges();

  double snow_depth;
  double melt;
  double runoff;
  double smb;
};

struct DEBMSimpleOrbitalParameters {
  double declination;
  double distance_factor;
};

//! A dEBM-simple implementation
/*!
 * This class implements dEBM-simple, the simple diurnal energy balance model described in
 *
 * M. Zeitz, R. Reese, J. Beckmann, U. Krebs-Kanzow, and R. Winkelmann, “Impact of the
 * melt–albedo feedback on the future evolution of the Greenland Ice Sheet with
 * PISM-dEBM-simple,” The Cryosphere, vol. 15, Art. no. 12, Dec. 2021.
*/
class DEBMSimplePointwise {
public:
  DEBMSimplePointwise(const Context &ctx);

  double albedo(double melt_rate, MaskValue cell_type) const;

  DEBMSimpleOrbitalParameters orbital_parameters(double time) const;

  DEBMSimpleMelt melt(double declination,
                      double distance_factor,
                      double dt,
                      double T_std_deviation,
                      double T,
                      double surface_elevation,
                      double lat,
                      double albedo) const;

  DEBMSimpleChanges step(double ice_thickness, double max_melt, double snow_depth,
                         double accumulation) const;

  // public because it is a diagnostic field
  double atmosphere_transmissivity(double elevation) const;

  double insolation(double declination,
                    double distance_factor,
                    double latitude_degrees) const;

  // implementation details (exposed as "public" methods for testing)
  static double CalovGreveIntegrand(double sigma, double temperature);
  static double hour_angle(double phi, double latitude, double declination);
  static double solar_longitude(double year_fraction, double eccentricity,
                                double perihelion_longitude);
  static double distance_factor_present_day(double year_fraction);
  static double distance_factor_paleo(double eccentricity, double perihelion_longitude,
                                      double solar_longitude);
  static double solar_declination_present_day(double year_fraction);
  static double solar_declination_paleo(double obliquity,
                                        double solar_longitude);
  static double insolation(double solar_constant, double distance_factor, double hour_angle,
                           double latitude, double declination);

private:
  double eccentricity(double time) const;
  double obliquity(double time) const;
  double perihelion_longitude(double time) const;

  //! refreeze melted ice
  bool m_refreeze_ice_melt;
  //! refreeze fraction
  double m_refreeze_fraction;
  //! threshold temperature for the computation of temperature-driven melt
  double m_positive_threshold_temperature;

  double m_ice_density;
  double m_water_density;

  double m_albedo_max;
  double m_albedo_min;
  double m_albedo_ocean;

  //! slope used in the linear parameterization of the albedo as a function of melt
  double m_albedo_slope;

  //! slope used in the linear parameterization of transmissivity
  double m_transmissivity_slope;
  double m_transmissivity_intercept;

  // tuning parameters of the melt equation
  double m_melt_c1;
  double m_melt_c2;

  // threshold air temperature (no melt at temperatures below this)
  double m_melt_threshold_temp;

  //! latent heat of fusion
  double m_L;
  //! the solar constant
  double m_solar_constant;

  //! minimum solar elevation angle above which melt is possible
  double m_phi;

  std::unique_ptr<ScalarForcing> m_eccentricity;
  std::unique_ptr<ScalarForcing> m_obliquity;
  std::unique_ptr<ScalarForcing> m_perihelion_longitude;

  bool m_paleo;
  bool m_use_paleo_file;

  double m_constant_eccentricity;
  double m_constant_perihelion_longitude;
  double m_constant_obliquity;

  std::shared_ptr<const Time> m_time;
};

} // end of namespace surface
} // end of namespace pism

#endif  /* PISM_DEBM_SIMPLE_POINTWISE_H */
