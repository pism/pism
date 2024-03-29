// Copyright (C) 2008-2018, 2022, 2023 Ed Bueler, Constantine Khroulev, Ricarda Winkelmann,
// Gudfinna Adalgeirsdottir and Andy Aschwanden
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

// Implementation of the atmosphere model using constant-in-time precipitation
// and a cosine yearly cycle for near-surface air temperatures.

// This includes the SeaRISE Greenland parameterization.

#include "pism/coupler/atmosphere/SeariseGreenland.hh"
#include "pism/geometry/Geometry.hh"
#include "pism/util/ConfigInterface.hh"
#include "pism/util/Grid.hh"
#include "pism/util/MaxTimestep.hh"
#include "pism/util/error_handling.hh"

namespace pism {
namespace atmosphere {

///// SeaRISEGreenland

SeaRISEGreenland::SeaRISEGreenland(std::shared_ptr<const Grid> g)
  : YearlyCycle(g) {
  // empty
}

SeaRISEGreenland::~SeaRISEGreenland() {
}

void SeaRISEGreenland::init_impl(const Geometry &geometry) {

  m_log->message(2,
                 "* Initializing SeaRISE-Greenland atmosphere model based on the Fausto et al (2009)\n"
                 "  air temperature parameterization and using stored time-independent precipitation...\n");

  m_reference =
    "R. S. Fausto, A. P. Ahlstrom, D. V. As, C. E. Boggild, and S. J. Johnsen, 2009. "
    "A new present-day temperature parameterization for Greenland. J. Glaciol. 55 (189), 95-105.";

  auto precip_file = m_config->get_string("atmosphere.searise_greenland.file");

  if (not precip_file.empty()) {
    m_log->message(2,
                   "  * Reading precipitation from '%s'...\n",
                   precip_file.c_str());

    YearlyCycle::init_internal(precip_file,
                               true, /* do regrid */
                               0 /* start (irrelevant) */);
  } else {
    YearlyCycle::init_impl(geometry);
  }
}

void SeaRISEGreenland::precip_time_series_impl(int i, int j, std::vector<double> &result) const {

  for (unsigned int k = 0; k < m_ts_times.size(); k++) {
    result[k] = m_precipitation(i,j);
  }
}

MaxTimestep SeaRISEGreenland::max_timestep_impl(double t) const {
  (void) t;
  return MaxTimestep("atmosphere searise_greenland");
}

//! \brief Updates mean annual and mean summer (July) near-surface air temperatures.
//! Note that the precipitation rate is time-independent and does not need
//! to be updated.
void SeaRISEGreenland::update_impl(const Geometry &geometry, double t, double dt) {
  (void) t;
  (void) dt;

  const double
    d_ma     = m_config->get_number("atmosphere.fausto_air_temp.d_ma"),      // K
    gamma_ma = m_config->get_number("atmosphere.fausto_air_temp.gamma_ma"),  // K m-1
    c_ma     = m_config->get_number("atmosphere.fausto_air_temp.c_ma"),      // K (degN)-1
    kappa_ma = m_config->get_number("atmosphere.fausto_air_temp.kappa_ma"),  // K (degW)-1
    d_mj     = m_config->get_number("atmosphere.fausto_air_temp.d_mj"),      // SAME UNITS as for _ma ...
    gamma_mj = m_config->get_number("atmosphere.fausto_air_temp.gamma_mj"),
    c_mj     = m_config->get_number("atmosphere.fausto_air_temp.c_mj"),
    kappa_mj = m_config->get_number("atmosphere.fausto_air_temp.kappa_mj");


  // initialize pointers to fields the parameterization depends on:
  const array::Scalar
    &h        = geometry.ice_surface_elevation,
    &lat_degN = geometry.latitude,
    &lon_degE = geometry.longitude;

  if (lat_degN.metadata().has_attribute("missing_at_bootstrap")) {
    throw RuntimeError(PISM_ERROR_LOCATION, "latitude variable was missing at bootstrap;\n"
                       "SeaRISE-Greenland atmosphere model depends on latitude and would return nonsense!");
  }

  if (lon_degE.metadata().has_attribute("missing_at_bootstrap")) {
    throw RuntimeError(PISM_ERROR_LOCATION, "longitude variable was missing at bootstrap;\n"
                       "SeaRISE-Greenland atmosphere model depends on longitude and would return nonsense!");
  }

  array::AccessScope list{&h, &lat_degN, &lon_degE, &m_air_temp_mean_annual, &m_air_temp_mean_summer};

  for (auto p = m_grid->points(); p; p.next()) {
    const int i = p.i(), j = p.j();
    m_air_temp_mean_annual(i,j) = d_ma + gamma_ma * h(i,j) + c_ma * lat_degN(i,j) + kappa_ma * (-lon_degE(i,j));
    m_air_temp_mean_summer(i,j)   = d_mj + gamma_mj * h(i,j) + c_mj * lat_degN(i,j) + kappa_mj * (-lon_degE(i,j));
  }
}

} // end of namespace atmosphere
} // end of namespace pism
