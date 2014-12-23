// Copyright (C) 2008-2014 Ed Bueler, Constantine Khroulev, Ricarda Winkelmann,
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

#include "PASeariseGreenland.hh"
#include "PISMVars.hh"
#include "IceGrid.hh"
#include "pism_options.hh"
#include "PISMTime.hh"
#include <assert.h>
#include "PISMConfig.hh"

#include "error_handling.hh"

namespace pism {

///// PA_SeaRISE_Greenland

PA_SeaRISE_Greenland::PA_SeaRISE_Greenland(const IceGrid &g)
  : PAYearlyCycle(g) {
  // empty
}

PA_SeaRISE_Greenland::~PA_SeaRISE_Greenland() {
}

void PA_SeaRISE_Greenland::init() {

  m_t = m_dt = GSL_NAN;  // every re-init restarts the clock

  verbPrintf(2, m_grid.com,
             "* Initializing SeaRISE-Greenland atmosphere model based on the Fausto et al (2009)\n"
             "  air temperature parameterization and using stored time-independent precipitation...\n");

  m_reference =
    "R. S. Fausto, A. P. Ahlstrom, D. V. As, C. E. Boggild, and S. J. Johnsen, 2009. "
    "A new present-day temperature parameterization for Greenland. J. Glaciol. 55 (189), 95-105.";

  std::string option_prefix = "-atmosphere_searise_greenland";
  options::String precip_file(option_prefix + "_file",
                              "Specifies a file with boundary conditions");
  m_precip_filename = precip_file;

  if (precip_file.is_set()) {
    verbPrintf(2, m_grid.com,
               "  * Option '-atmosphere_searise_greenland %s' is set...\n",
               m_precip_filename.c_str());

    PAYearlyCycle::init_internal(m_precip_filename,
                                 true, /* do regrid */
                                 0 /* start (irrelevant) */);
  } else {
    PAYearlyCycle::init();
  }

  // initialize pointers to fields the parameterization depends on:
  m_surfelev = m_grid.variables().get_2d_scalar("surface_altitude");
  m_lat      = m_grid.variables().get_2d_scalar("latitude");
  m_lon      = m_grid.variables().get_2d_scalar("longitude");
}

void PA_SeaRISE_Greenland::precip_time_series(int i, int j, std::vector<double> &result) {

  for (unsigned int k = 0; k < m_ts_times.size(); k++) {
    result[k] = m_precipitation(i,j);
  }
}

//! \brief Updates mean annual and mean July near-surface air temperatures.
//! Note that the precipitation rate is time-independent and does not need
//! to be updated.
void PA_SeaRISE_Greenland::update(double my_t, double my_dt) {

  if (m_lat->metadata().has_attribute("missing_at_bootstrap")) {
    throw RuntimeError("latitude variable was missing at bootstrap;\n"
                       "SeaRISE-Greenland atmosphere model depends on latitude and would return nonsense!");
  }

  if (m_lon->metadata().has_attribute("missing_at_bootstrap")) {
    throw RuntimeError("longitude variable was missing at bootstrap;\n"
                       "SeaRISE-Greenland atmosphere model depends on longitude and would return nonsense!");
  }

  if ((fabs(my_t - m_t) < 1e-12) &&
      (fabs(my_dt - m_dt) < 1e-12)) {
    return;
  }

  m_t  = my_t;
  m_dt = my_dt;

  const double 
    d_ma     = m_config.get("snow_temp_fausto_d_ma"),      // K
    gamma_ma = m_config.get("snow_temp_fausto_gamma_ma"),  // K m-1
    c_ma     = m_config.get("snow_temp_fausto_c_ma"),      // K (degN)-1
    kappa_ma = m_config.get("snow_temp_fausto_kappa_ma"),  // K (degW)-1
    d_mj     = m_config.get("snow_temp_fausto_d_mj"),      // SAME UNITS as for _ma ...
    gamma_mj = m_config.get("snow_temp_fausto_gamma_mj"),
    c_mj     = m_config.get("snow_temp_fausto_c_mj"),
    kappa_mj = m_config.get("snow_temp_fausto_kappa_mj");

  IceModelVec2S &h = *m_surfelev, &lat_degN = *m_lat, &lon_degE = *m_lon;

  IceModelVec::AccessList list;
  list.add(h);
  list.add(lat_degN);
  list.add(lon_degE);
  list.add(m_air_temp_mean_annual);
  list.add(m_air_temp_mean_july);

  for (Points p(m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();
    m_air_temp_mean_annual(i,j) = d_ma + gamma_ma * h(i,j) + c_ma * lat_degN(i,j) + kappa_ma * (-lon_degE(i,j));
    m_air_temp_mean_july(i,j)   = d_mj + gamma_mj * h(i,j) + c_mj * lat_degN(i,j) + kappa_mj * (-lon_degE(i,j));
  }
}

} // end of namespace pism
