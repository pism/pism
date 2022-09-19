// Copyright (C) 2009--2022 PISM Authors
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

#include "pism/util/Mask.hh"
#include "pism/util/ScalarForcing.hh"

namespace pism {

class Context;
class Time;

namespace surface {

//! A dEBM-simple implementation
/*!
  after Krebs-Kanzow et al. 2018
*/
class DEBMSimplePointwise {
public:
  DEBMSimplePointwise(const Context &ctx);

  unsigned int timeseries_length(double dt);

  double albedo(double melt_rate, MaskValue cell_type);

  class Melt {
  public:
    Melt();

    double temperature_melt;
    double insolation_melt;
    double background_melt;
    double total_melt;
    double transmissivity;
    double q_insol;
  };

  Melt calculate_melt(double time,
                      double dt,
                      double T_std_deviation,
                      double T,
                      double surface_elevation,
                      double lat,
                      double albedo);

  void get_snow_accumulation(const std::vector<double> &T, std::vector<double> &precip_rate);

  class Changes {
  public:
    Changes();

    double firn_depth;
    double snow_depth;
    double melt;
    double runoff;
    double smb;
  };

  Changes step(double ice_thickness, double max_melt, double firn_depth,
               double snow_depth, double accumulation);

  double eccentricity(double time);
  double obliquity(double time);
  double perihelion_longitude(double time);

  double distance_factor(double time);
  double distance_factor_paleo(double time);

  double solar_declination(double time);
  double solar_declination_paleo(double time);

  double atmosphere_transmissivity(double elevation);

private:
  double solar_longitude(double year_fraction,
                         double eccentricity,
                         double perihelion_longitude);

  double get_h_phi(double phi, double lat, double delta);

  double get_q_insol(double distance2, double h_phi, double lat, double delta);

  double CalovGreveIntegrand(double sigma, double TacC);
  //! interpret all the precipitation as snow (no rain)
  bool m_precip_as_snow;
  //! refreeze melted ice
  bool m_refreeze_ice_melt;
  //! refreeze fraction
  double m_refreeze_fraction;
  //! the temperature below which all precipitation is snow
  double m_Tmin;
  //! the temperature above which all precipitation is rain
  double m_Tmax;
  //! threshold temperature for the computation of temperature-driven melt
  double m_positive_threshold_temperature;

  //! year length used to compute the time series length required to get m_n_per_year
  //! evaluations
  double m_year_length;

  //! number of small time steps per year
  unsigned int m_n_per_year;

  double m_ice_density;
  double m_water_density;

  double m_albedo_snow;
  double m_albedo_ice;
  double m_albedo_land;
  double m_albedo_ocean;

  //! slope used in the linear parameterization of the albedo as a function of melt
  double m_albedo_slope;

  //! slope used in the linear parameterization of transmissivity
  double m_tau_a_slope;
  double m_tau_a_intercept;

  // tuning parameters of the melt equation
  double m_c1;
  double m_c2;

  double m_bm_temp;

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
