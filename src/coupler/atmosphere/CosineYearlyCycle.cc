// Copyright (C) 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2021, 2023, 2025 PISM Authors
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

#include "pism/coupler/atmosphere/CosineYearlyCycle.hh"
#include "pism/util/Time.hh"
#include "pism/util/pism_options.hh"
#include "pism/util/Config.hh"
#include "pism/util/Grid.hh"
#include "pism/util/io/File.hh"

#include "pism/util/error_handling.hh"
#include "pism/util/MaxTimestep.hh"

#include "pism/util/ScalarForcing.hh"
#include "pism/util/Logger.hh"

namespace pism {
namespace atmosphere {

CosineYearlyCycle::CosineYearlyCycle(std::shared_ptr<const Grid> grid)
  : YearlyCycle(grid) {

  auto scaling_file = m_config->get_string("atmosphere.yearly_cycle.scaling.file");

  if (not scaling_file.empty()) {
    m_A.reset(new ScalarForcing(*grid->ctx(),
                                "atmosphere.yearly_cycle.scaling",
                                "amplitude_scaling",
                                "1", "1",
                                "temperature amplitude scaling"));
  }
}

void CosineYearlyCycle::init_impl(const Geometry &geometry) {
  (void) geometry;

  m_log->message(2,
                 "* Initializing the 'cosine yearly cycle' atmosphere model...\n");

  auto input_file   = m_config->get_string("atmosphere.yearly_cycle.file");

  if (input_file.empty()) {
    throw RuntimeError(PISM_ERROR_LOCATION,
                       "Please specify an '-atmosphere yearly_cycle' input file\n"
                       "using atmosphere.yearly_cycle.file or a command-line option.");
  }

  m_log->message(2,
                 "  Reading mean annual air temperature, mean July air temperature, and\n"
                 "  precipitation fields from '%s'...\n", input_file.c_str());

  m_air_temp_mean_annual.regrid(input_file, io::Default::Nil());
  m_air_temp_mean_summer.regrid(input_file, io::Default::Nil());
  m_precipitation.regrid(input_file, io::Default::Nil());
}

MaxTimestep CosineYearlyCycle::max_timestep_impl(double t) const {
  (void) t;
  return MaxTimestep("atmosphere cosine_yearly_cycle");
}

void CosineYearlyCycle::update_impl(const Geometry &geometry, double t, double dt) {
  (void) geometry;
  (void) t;
  (void) dt;
  // an implementation is necessary because the base class does not define this
}

void CosineYearlyCycle::init_timeseries_impl(const std::vector<double> &ts) const {

  YearlyCycle::init_timeseries_impl(ts);

  if (m_A) {
    for (unsigned int k = 0; k < ts.size(); ++k) {
      m_cosine_cycle[k] *= m_A->value(ts[k]);
    }
  }
}

} // end of namespace atmosphere
} // end of namespace pism
