// Copyright (C) 2009-2018 Ricarda Winkelmann, Torsten Albrecht, Constantine Khrulev
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

// Implementation of the atmosphere model using constant-in-time precipitation
// and a cosine yearly cycle for near-surface air temperatures.

// This includes the PIK temperature parameterization.

#include "PIK.hh"

#include "pism/geometry/Geometry.hh"
#include "pism/util/ConfigInterface.hh"
#include "pism/util/IceGrid.hh"
#include "pism/util/MaxTimestep.hh"
#include "pism/util/error_handling.hh"

namespace pism {
namespace atmosphere {

PIK::PIK(IceGrid::ConstPtr g)
  : YearlyCycle(g) {

  auto parameterization = m_config->get_string("atmosphere.pik.parameterization");

  std::map<std::string, Parameterization>
    models = {{"martin",                    MARTIN},
              {"huybrechts_dewolde",        HUYBRECHTS_DEWOLDE},
              {"martin_huybrechts_dewolde", MARTIN_HUYBRECHTS_DEWOLDE},
              {"era_interim",               ERA_INTERIM},
              {"era_interim_sin",           ERA_INTERIM_SIN},
              {"era_interim_lon",           ERA_INTERIM_LON}};

  if (models.find(parameterization) == models.end()) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "invalid pik parameterization: %s",
                                  parameterization.c_str());
  }

  m_parameterization = models[parameterization];
}

PIK::~PIK() {
  // empty
}

void PIK::init_impl(const Geometry &geometry) {

  m_log->message(2,
                 "* Initializing PIK atmosphere model with air temperature parameterization based on \n"
                 "  Huybrechts & De Wolde (1999) or multiple regression analysis of ERA INTERIM data...\n"
                 "  Uses a time-independent precipitation field read from a file.");

  m_reference = "Winkelmann et al.";

  auto precip_file = m_config->get_string("atmosphere.pik.file");

  if (not precip_file.empty()) {
    YearlyCycle::init_internal(precip_file,
                               true, /* do regrid */
                               0 /* start (irrelevant) */);
  } else {
    // try to read precipitation from the input (-i) file
    YearlyCycle::init_impl(geometry);
  }

  switch (m_parameterization) {
  case HUYBRECHTS_DEWOLDE:
    m_log->message(2,
                   "    Parameterization based on: Huybrechts & De Wolde (1999).\n");
    break;
  case ERA_INTERIM:
    m_log->message(2,
                   "    Parameterization based on: multiple regression analysis of ERA INTERIM data.\n");
    break;
  case ERA_INTERIM_SIN:
    m_log->message(2,
                   "    Parameterization based on: multiple regression analysis of ERA INTERIM data with a sin(lat) dependence.\n");
    break;
  case ERA_INTERIM_LON:
    m_log->message(2,
                   "    Parameterization based on: multiple regression analysis of ERA INTERIM data with a cos(lon) dependence.\n");
    break;
  case MARTIN_HUYBRECHTS_DEWOLDE:
    m_log->message(2,
                   "    Mean annual temperature: as in Martin et al (2011).\n"
                   "    Mean summer temperature: anomaly to the parameterization used by Huybrechts & De Wolde (1999).\n");
    break;
  default:
  case MARTIN:
    m_log->message(2,
                   "    Mean annual temperature: as in Martin et al (2011).\n"
                   "    No seasonal variation in air temperature.\n");
  }
}

MaxTimestep PIK::max_timestep_impl(double t) const {
  (void) t;
  return MaxTimestep("atmosphere pik");
}

/*!
 * See equation C1 in HuybrechtsdeWolde.
 */
static double huybrechts_dewolde_mean_annual(double surface_elevation, double latitude) {
  double gamma_a = surface_elevation < 1500.0 ? -0.005102 : -0.014285;
  return 273.15 + 34.46 + gamma_a * surface_elevation - 0.68775 * latitude * (-1.0);
}

/*!
 * See equation C2 in HuybrechtsdeWolde.
 */
static double huybrechts_dewolde_mean_summer(double surface_elevation, double latitude) {
  return 273.15 + 16.81 - 0.00692 * surface_elevation - 0.27937 * latitude * (-1.0);
}

/*!
 * Parameterization of mean annual and mean summer near-surface temperature as in
 * Huybrechts & DeWolde (1999)
 */
static void huybrechts_dewolde(const Geometry &geometry, IceModelVec2S &T_ma, IceModelVec2S &T_ms) {
  IceGrid::ConstPtr grid = T_ma.grid();

  const IceModelVec2S
    &h   = geometry.ice_surface_elevation,
    &lat = geometry.latitude;

  IceModelVec::AccessList list{&h, &lat, &T_ma, &T_ms};

  for (Points p(*grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    T_ma(i, j) = huybrechts_dewolde_mean_annual(h(i, j), lat(i, j));
    T_ms(i, j) = huybrechts_dewolde_mean_summer(h(i, j), lat(i, j));
  }
}

/*!
 * Parametrization based on multiple regression analysis of ERA INTERIM data
 */
static void era_interim(const Geometry &geometry, IceModelVec2S &T_ma, IceModelVec2S &T_ms) {
  IceGrid::ConstPtr grid = T_ma.grid();

  const IceModelVec2S
    &h   = geometry.ice_surface_elevation,
    &lat = geometry.latitude;

  IceModelVec::AccessList list{&h, &lat, &T_ma, &T_ms};

  for (Points p(*grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    T_ma(i, j) = 273.15 + 29.2 - 0.0082 * h(i, j) - 0.576 * lat(i, j) * (-1.0);
    T_ms(i, j) = 273.15 + 16.5 - 0.0068 * h(i, j) - 0.248 * lat(i, j) * (-1.0);
  }
}

/*!
 * Parametrization based on multiple regression analysis of ERA INTERIM data with sin(lat)
 */
static void era_interim_sin(const Geometry &geometry, IceModelVec2S &T_ma, IceModelVec2S &T_ms) {
  IceGrid::ConstPtr grid = T_ma.grid();

  const IceModelVec2S
    &h   = geometry.ice_surface_elevation,
    &lat = geometry.latitude;

  IceModelVec::AccessList list{&h, &lat, &T_ma, &T_ms};

  for (Points p(*grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    T_ma(i, j) = 273.15 - 2.0 - 0.0082 * h(i, j) + 18.4 * (sin(3.1415 * lat(i, j) / 180.0) + 0.8910) / (1 - 0.8910);
    T_ms(i, j) = 273.15 + 3.2 - 0.0067 * h(i, j) + 8.3 * (sin(3.1415 * lat(i, j) / 180.0) + 0.8910) / (1 - 0.8910);
  }
}

/*!
 * Parametrization based on multiple regression analysis of ERA INTERIM data with cos(lon)
 */
static void era_interim_lon(const Geometry &geometry, IceModelVec2S &T_ma, IceModelVec2S &T_ms) {
  IceGrid::ConstPtr grid = T_ma.grid();

  const IceModelVec2S
    &h   = geometry.ice_surface_elevation,
    &lat = geometry.latitude,
    &lon = geometry.longitude;

  IceModelVec::AccessList list{&h, &lat, &lon, &T_ma, &T_ms};

  for (Points p(*grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    double hmod = std::max(1000.0, h(i, j));
    T_ma(i, j) = 273.15 + 37.49 - 0.0095 * hmod - 0.644 * lat(i, j) * (-1.0) + 2.146 * cos(3.1415 * (lon(i, j) + 110.0) / 180.0);
    T_ms(i, j) = 273.15 + 15.74 - 0.0083 * hmod - 0.196 * lat(i, j) * (-1.0) + 0.225 * cos(3.1415 * (lon(i, j) + 110.0) / 180.0);

  }
}

/*!
 * Parameterization of the mean annual near-surface air temperature, see equation (1) in
 * Martin et al, 2011.
 */
static double martin2011_mean_annual(double elevation, double latitude) {
  return 273.15 + 30 - 0.0075 * elevation - 0.68775 * latitude * (-1.0);
}

/*!
 * - annual mean temperature as in Martin et al. (2011)
 * - no seasonal variation of air temperature
 */
static void martin2011(const Geometry &geometry, IceModelVec2S &T_ma, IceModelVec2S &T_ms) {
  IceGrid::ConstPtr grid = T_ma.grid();

  const IceModelVec2S
    &h   = geometry.ice_surface_elevation,
    &lat = geometry.latitude;

  IceModelVec::AccessList list{&h, &lat, &T_ma, &T_ms};

  for (Points p(*grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    T_ma(i, j) = martin2011_mean_annual(h(i, j), lat(i, j));
    T_ms(i, j) = T_ma(i, j);
  }
}

/*!
 * - annual mean temperature as in Martin et al. (2011)
 * - summer mean temperature computed as an anomaly to Huybrechts & DeWolde (1999)
 */
static void martin_huybrechts_dewolde(const Geometry &geometry, IceModelVec2S &T_ma, IceModelVec2S &T_ms) {
  IceGrid::ConstPtr grid = T_ma.grid();

  const IceModelVec2S
    &h   = geometry.ice_surface_elevation,
    &lat = geometry.latitude;

  IceModelVec::AccessList list{&h, &lat, &T_ma, &T_ms};

  for (Points p(*grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    // mean annual surface temperature as in Martin et al. 2011, equation (1)
    T_ma(i, j) = 273.15 + 30 - 0.0075 * h(i, j) - 0.68775 * lat(i, j) * (-1.0);

    double
      TMA = huybrechts_dewolde_mean_annual(h(i, j), lat(i, j)),
      TMS = huybrechts_dewolde_mean_summer(h(i, j), lat(i, j));

    T_ms(i, j) = T_ma(i, j) + (TMS - TMA);
  }
}


/*!
 * Updates mean annual and mean summer (January) near-surface air temperatures. Note that
 * the precipitation rate is time-independent and does not need to be updated.
 */
void PIK::update_impl(const Geometry &geometry, double t, double dt) {
  (void) t;
  (void) dt;

  if (geometry.latitude.metadata().has_attribute("missing_at_bootstrap")) {
    throw RuntimeError(PISM_ERROR_LOCATION,
                       "latitude variable was missing at bootstrap;\n"
                       "PIK atmosphere model depends on latitude and would return nonsense!");
  }

  switch (m_parameterization) {
  case HUYBRECHTS_DEWOLDE:
    huybrechts_dewolde(geometry, m_air_temp_mean_annual, m_air_temp_mean_summer);
    break;
  case ERA_INTERIM:
    era_interim(geometry, m_air_temp_mean_annual, m_air_temp_mean_summer);
    break;
  case ERA_INTERIM_SIN:
    era_interim_sin(geometry, m_air_temp_mean_annual, m_air_temp_mean_summer);
    break;
  case ERA_INTERIM_LON:
    era_interim_lon(geometry, m_air_temp_mean_annual, m_air_temp_mean_summer);
    break;
  case MARTIN_HUYBRECHTS_DEWOLDE:
    martin_huybrechts_dewolde(geometry, m_air_temp_mean_annual, m_air_temp_mean_summer);
    break;
  default:
  case MARTIN:
    martin2011(geometry, m_air_temp_mean_annual, m_air_temp_mean_summer);
    break;
  }
}

} // end of namespace atmosphere
} // end of namespace pism
