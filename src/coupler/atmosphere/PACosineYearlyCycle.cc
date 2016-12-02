// Copyright (C) 2012, 2013, 2014, 2015, 2016 PISM Authors
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

#include <gsl/gsl_math.h>       // M_PI, GSL_NAN

#include "PACosineYearlyCycle.hh"
#include "base/util/Timeseries.hh"
#include "base/util/PISMTime.hh"
#include "base/util/pism_options.hh"
#include "base/util/PISMConfigInterface.hh"
#include "base/util/IceGrid.hh"
#include "base/util/io/PIO.hh"

#include "base/util/error_handling.hh"
#include "base/util/MaxTimestep.hh"

namespace pism {
namespace atmosphere {

CosineYearlyCycle::CosineYearlyCycle(IceGrid::ConstPtr g)
  : YearlyCycle(g), m_A(NULL) {
}

CosineYearlyCycle::~CosineYearlyCycle() {
  if (m_A != NULL) {
    delete m_A;
  }
}

void CosineYearlyCycle::init_impl() {

  m_t = m_dt = GSL_NAN;  // every re-init restarts the clock

  m_log->message(2,
             "* Initializing the 'cosine yearly cycle' atmosphere model (-atmosphere yearly_cycle)...\n");


  options::String input_file("-atmosphere_yearly_cycle_file",
                             "CosineYearlyCycle input file name");
  options::String scaling_file("-atmosphere_yearly_cycle_scaling_file",
                               "CosineYearlyCycle amplitude scaling input file name");

  if (not input_file.is_set()) {
    throw RuntimeError(PISM_ERROR_LOCATION, "Please specify an '-atmosphere yearly_cycle' input file\n"
                       "using the -atmosphere_yearly_cycle_file option.");
  }

  m_log->message(2,
             "  Reading mean annual air temperature, mean July air temperature, and\n"
             "  precipitation fields from '%s'...\n", input_file->c_str());

  m_air_temp_mean_annual.regrid(input_file, CRITICAL);
  m_air_temp_mean_july.regrid(input_file, CRITICAL);
  m_precipitation.regrid(input_file, CRITICAL);

  if (scaling_file.is_set()) {

    if (m_A == NULL) {
      m_A = new Timeseries(*m_grid, "amplitude_scaling",
                           m_config->get_string("time.dimension_name"));
      m_A->metadata().set_string("units", "1");
      m_A->metadata().set_string("long_name", "cosine yearly cycle amplitude scaling");
      m_A->dimension_metadata().set_string("units", m_grid->ctx()->time()->units_string());
    }

    m_log->message(2,
               "  Reading cosine yearly cycle amplitude scaling from '%s'...\n",
               scaling_file->c_str());

    PIO nc(m_grid->com, "netcdf3", scaling_file, PISM_READONLY);    // OK to use netcdf3
    {
      m_A->read(nc, *m_grid->ctx()->time(), *m_grid->ctx()->log());
    }
    nc.close();

  } else {
    if (m_A != NULL) {
      delete m_A;
    }
    m_A = NULL;
  }
}

MaxTimestep CosineYearlyCycle::max_timestep_impl(double t) const {
  (void) t;
  return MaxTimestep("atmosphere cosine_yearly_cycle");
}

void CosineYearlyCycle::update_impl(double my_t, double my_dt) {
  m_t = my_t;
  m_dt = my_dt;
}

void CosineYearlyCycle::init_timeseries_impl(const std::vector<double> &ts) const {

  YearlyCycle::init_timeseries_impl(ts);

  if (m_A != NULL) {
    for (unsigned int k = 0; k < ts.size(); ++k) {
      m_cosine_cycle[k] *= (*m_A)(ts[k]);
    }
  }
}

} // end of namespace atmosphere
} // end of namespace pism
