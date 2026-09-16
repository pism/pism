// Copyright (C) 2009--2026 PISM Authors
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

class DEBMSimpleAtmosphereTransmissivity {
public:
  DEBMSimpleAtmosphereTransmissivity(const Config &config);
  double operator()(double surface_elevation) const;
private:
  double m_slope, m_intercept;
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

  DEBMSimpleMelt melt(double declination,
                      double distance_factor,
                      double dt,
                      double T_std_deviation,
                      double T,
                      double surface_elevation,
                      double lat,
                      double albedo) const;

  //! Insolation energy (J m^-2) reaching the surface over `dt`, computed from the
  //! analytic top-of-atmosphere insolation used by dEBM-simple.
  double insolation_energy(double declination,
                           double distance_factor,
                           double latitude,
                           double dt) const;

  //! Melt (and its components) given the insolation *energy* (J m^-2) reaching the
  //! surface over `dt`. dEBM-enhanced supplies a prescribed insolation field here instead
  //! of the analytic one. `declination` and `latitude` set the daily melt-period length
  //! that weights the temperature- and offset-driven melt terms.
  DEBMSimpleMelt melt_from_insolation(double declination,
                                      double latitude,
                                      double insolation_energy,
                                      double dt,
                                      double T_std_deviation,
                                      double T,
                                      double surface_elevation,
                                      double albedo) const;

  DEBMSimpleChanges step(double ice_thickness, double max_melt, double snow_depth,
                         double accumulation) const;

  // public because it is a diagnostic field
  double atmosphere_transmissivity(double elevation) const;

  double insolation_diagnostic(double declination,
                    double distance_factor,
                    double latitude_degrees) const;

  // implementation details (exposed as "public" methods for testing)
  static double CalovGreveIntegrand(double sigma, double temperature);
  static double hour_angle(double phi, double latitude, double declination);
  static double insolation_rate(double solar_constant, double distance_factor, double hour_angle,
                                double latitude, double declination);

private:
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

  // atmosphere transmissivity model
  DEBMSimpleAtmosphereTransmissivity m_transmissivity;
};

} // end of namespace surface
} // end of namespace pism

#endif  /* PISM_DEBM_SIMPLE_POINTWISE_H */
