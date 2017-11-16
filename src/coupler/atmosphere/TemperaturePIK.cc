// Copyright (C)  2009-2017 Ricarda Winkelmann, Torsten Albrecht
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

#include <gsl/gsl_math.h>

#include "TemperaturePIK.hh"
#include "pism/util/Vars.hh"
#include "pism/util/IceGrid.hh"
#include "pism/util/pism_options.hh"
#include "pism/util/Time.hh"
#include <cassert>
#include "pism/util/ConfigInterface.hh"
#include "pism/util/MaxTimestep.hh"
#include "pism/util/error_handling.hh"


namespace pism {
namespace atmosphere {


///// TemperaturePIK
TemperaturePIK::TemperaturePIK(IceGrid::ConstPtr g)
  : YearlyCycle(g) 
{
    //empty
}

TemperaturePIK::~TemperaturePIK() 
{
    //empty
}



void TemperaturePIK::init_impl() {

  m_t = m_dt = GSL_NAN;  // every re-init restarts the clock

  m_log->message(2,
             "* Initializing PIK atmosphere model with air temperature parameterization based on \n"
             "  Huybrechts & De Wolde (1999) and/or multiple regression analysis of ERA INTERIM data...\n"
             "  Precipitation is per default time-independent, use modifiers paleo_precip, delta_P or lapse_rate!\n");

  m_reference =
    "Winkelmann et al.";


  std::string option_prefix = "-atmosphere_pik_temp";
  options::String precip_file(option_prefix + "_file",
                              "Specifies a file with boundary conditions");

  if (precip_file.is_set()) {
    m_log->message(2,
                   "  * Option '-atmosphere_pik_temp %s' is set...\n",
                   precip_file->c_str());

    YearlyCycle::init_internal(precip_file,
                               true, /* do regrid */
                               0 /* start (irrelevant) */);
  } else {
    YearlyCycle::init_impl();
  }



  /// Surface (annual mean and summer mean) temperature parametrization: 
  temp_huybrechts_dewolde99_set = options::Bool("-temp_huybrechts_dewolde99",
    "Near-surface air temperature is parameterized as in Huybrechts & De Wolde (1999)");

  temp_era_interim_set = options::Bool("-temp_era_interim",
    "Near-surface air temperature is parameterized based on ERA INTERIM data");

  temp_era_interim_sin_set = options::Bool("-temp_era_interim_sin",
    "Near-surface air temperature is parameterized based on ERA INTERIM data with a sin(lat) dependence");

  temp_era_interim_lon_set = options::Bool("-temp_era_interim_lon",
    "Near-surface air temperature is parameterized based on ERA INTERIM data with a cos(lon) dependence");


  if (temp_huybrechts_dewolde99_set) {
      m_log->message(2,
             "    Near-surface air temperature is parameterized as in Huybrechts & De Wolde (1999).\n");
  } else if (temp_era_interim_set) {
      m_log->message(2,
             "    Near-surface air temperature is parameterized based on ERA INTERIM data.\n");
  } else if (temp_era_interim_sin_set) {
      m_log->message(2,
             "    Near-surface air temperature is parameterized based on ERA INTERIM data with a sin(lat) dependence.\n");
  } else if (temp_era_interim_lon_set) {
      m_log->message(2,
             "    Near-surface air temperature is parameterized based on ERA INTERIM data with a cos(lon) dependence.\n");
  }else {
      m_log->message(2,
             "    Near-surface annual mean air temperature is parameterized as in Martin et al. (2011),\n"
             "    and near-surface summer mean air temperature is computed as anomaly to the Huybrechts & De Wolde (1999) - temperature.\n");
  }
}


MaxTimestep TemperaturePIK::max_timestep_impl(double t) const {
  (void) t;
  //return MaxTimestep();
  return MaxTimestep("atmosphere pik_temp");
}

void TemperaturePIK::precip_time_series_impl(int i, int j, std::vector<double> &result) const {

  for (unsigned int k = 0; k < m_ts_times.size(); k++) {
    result[k] = m_precipitation(i,j);
  }
}


//! \brief Updates mean annual and mean July near-surface air temperatures.
//! Note that the precipitation rate is time-independent and does not need
//! to be updated.
void TemperaturePIK::update_impl(double my_t, double my_dt) {

  if ((fabs(my_t - m_t) < 1e-12) &&
      (fabs(my_dt - m_dt) < 1e-12)) {
    return;
  }

  m_t  = my_t;
  m_dt = my_dt;

  // initialize pointers to fields the parameterization depends on:
  const IceModelVec2S
    &h        = *m_grid->variables().get_2d_scalar("surface_altitude"),
    &lat_degN = *m_grid->variables().get_2d_scalar("latitude"),
    &lon_degE = *m_grid->variables().get_2d_scalar("longitude");

  if (lat_degN.metadata().has_attribute("missing_at_bootstrap")) {
    throw RuntimeError(PISM_ERROR_LOCATION, "latitude variable was missing at bootstrap;\n"
                       "TemperaturePIK atmosphere model depends on latitude and would return nonsense!");
  }

  IceModelVec::AccessList list;
  list.add(h);
  list.add(lat_degN);
  list.add(lon_degE);
  list.add(m_air_temp_mean_annual);
  list.add(m_air_temp_mean_july); //FIXME: change name to summer, since it is january in Antarctica!

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    double gamma_a;
    if (h(i,j) < 1500.0) {
      gamma_a = -0.005102;
    }else{
      gamma_a = -0.014285;
    }
  
    if (temp_huybrechts_dewolde99_set){
      //m_air_temp_mean_annual(i,j) = d_ma + gamma_ma * h(i,j) + c_ma * lat_degN(i,j) + kappa_ma * (-lon_degE(i,j));
      //m_air_temp_mean_july(i,j)   = d_mj + gamma_mj * h(i,j) + c_mj * lat_degN(i,j) + kappa_mj * (-lon_degE(i,j));
      m_air_temp_mean_annual(i,j) = 273.15 + 34.46 + gamma_a * h(i,j) - 0.68775 * lat_degN(i,j)*(-1.0);
      m_air_temp_mean_july(i,j)   = 273.15 + 16.81 - 0.00692 * h(i,j) - 0.27937 * lat_degN(i,j)*(-1.0);
    }

    else if (temp_era_interim_set){  // parametrization based on multiple regression analysis of ERA INTERIM data
      m_air_temp_mean_annual(i,j) = 273.15 + 29.2 - 0.0082 * h(i,j) - 0.576 * lat_degN(i,j)*(-1.0);
      m_air_temp_mean_july(i,j)   = 273.15 + 16.5 - 0.0068 * h(i,j) - 0.248 * lat_degN(i,j)*(-1.0);
    }
    else if (temp_era_interim_sin_set){  // parametrization based on multiple regression analysis of ERA INTERIM data with sin(lat)
      m_air_temp_mean_annual(i,j) = 273.15 - 2.0 -0.0082*h(i,j) + 18.4 * (sin(3.1415*lat_degN(i,j)/180.0)+0.8910)/(1-0.8910);
      m_air_temp_mean_july(i,j)   = 273.15 + 3.2 -0.0067*h(i,j) +  8.3 * (sin(3.1415*lat_degN(i,j)/180.0)+0.8910)/(1-0.8910);
    }
    else if (temp_era_interim_lon_set){  // parametrization based on multiple regression analysis of ERA INTERIM data with cos(lon)
      double hmod = std::max(500.0,h(i,j));  //FIXME: if icefree ocean, hmod=0
      m_air_temp_mean_annual(i,j) = 273.15 + 36.81 -0.00797*hmod -0.688*lat_degN(i,j)*(-1.0) + 2.574*cos(3.1415*(lon_degE(i,j)-110.0)/180.0);
      m_air_temp_mean_july(i,j)   = 273.15 + 22.58 -0.00940*hmod -0.234*lat_degN(i,j)*(-1.0) + 0.828*cos(3.1415*(lon_degE(i,j)-110.0)/180.0);
    }
    else {

    // annual mean temperature = Martin et al. (2011) parametrization
    // summer mean temperature = anomaly to Huybrechts & DeWolde (1999)
      m_air_temp_mean_annual(i,j) = 273.15 + 30 - 0.0075 * h(i,j) - 0.68775 * lat_degN(i,j)*(-1.0);  // surface temperature parameterization as in Martin et al. 2011, Eqn. 2.0.2

      double TMA = 273.15 + 34.46 + gamma_a * h(i,j) - 0.68775 * lat_degN(i,j)*(-1.0); // = TMA, mean annual temperature in Huybrechts & DeWolde (1999)
      double TMS = 273.15 + 16.81 - 0.00692 * h(i,j) - 0.27937 * lat_degN(i,j)*(-1.0); // = TMS, mean summer temperature in Huybrechts & DeWolde (1999)

      m_air_temp_mean_july(i,j) = m_air_temp_mean_annual(i,j) + (TMS - TMA);
    }

  }
}

} // end of namespace atmosphere
} // end of namespace pism

