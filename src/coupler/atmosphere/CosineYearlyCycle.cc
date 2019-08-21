// Copyright (C) 2012, 2013, 2014, 2015, 2016, 2017, 2018 PISM Authors
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

#include "CosineYearlyCycle.hh"
#include "pism/util/Timeseries.hh"
#include "pism/util/Time.hh"
#include "pism/util/pism_options.hh"
#include "pism/util/ConfigInterface.hh"
#include "pism/util/IceGrid.hh"
#include "pism/util/io/PIO.hh"

#include "pism/util/error_handling.hh"
#include "pism/util/MaxTimestep.hh"

namespace pism {
namespace atmosphere {

CosineYearlyCycle::CosineYearlyCycle(IceGrid::ConstPtr g)
  : YearlyCycle(g) {
  
}

CosineYearlyCycle::~CosineYearlyCycle() {
  // empty
}

void CosineYearlyCycle::init_impl(const Geometry &geometry) {
  (void) geometry;

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
             "  Reading mean annual air temperature, mean summer air temperature (NH:July, SH:January), and\n"
             "  precipitation fields from '%s'...\n", input_file->c_str());

  m_air_temp_mean_annual.regrid(input_file, CRITICAL);
  m_air_temp_mean_summer.regrid(input_file, CRITICAL);
  m_precipitation.regrid(input_file, CRITICAL);

  if (scaling_file.is_set()) {

    m_A.reset(new Timeseries(*m_grid, "amplitude_scaling",
                             m_config->get_string("time.dimension_name")));
    m_A->variable().set_string("units", "1");
    m_A->variable().set_string("long_name", "cosine yearly cycle amplitude scaling");
    m_A->dimension().set_string("units", m_grid->ctx()->time()->units_string());

    m_log->message(2,
                   "  Reading cosine yearly cycle amplitude scaling from '%s'...\n",
                   scaling_file->c_str());

    PIO nc(m_grid->com, "netcdf3", scaling_file, PISM_READONLY);    // OK to use netcdf3
    {
      m_A->read(nc, *m_grid->ctx()->time(), *m_grid->ctx()->log());
    }
    nc.close();

  } else {
    m_A = NULL;
  }
}

MaxTimestep CosineYearlyCycle::max_timestep_impl(double t) const {
  (void) t;
  return MaxTimestep("atmosphere cosine_yearly_cycle");
}

void CosineYearlyCycle::update_impl(const Geometry &geometry, double t, double dt) {
  (void) geometry;
  (void) t;
  (void) dt;
  // nothing to do here, but an implementation is required
}

void CosineYearlyCycle::init_timeseries_impl(const std::vector<double> &ts) const {

  YearlyCycle::init_timeseries_impl(ts);

  if (m_A) {
    for (unsigned int k = 0; k < ts.size(); ++k) {
      m_cosine_cycle[k] *= (*m_A)(ts[k]);
    }
  }
}

} // end of namespace atmosphere
} // end of namespace pism
