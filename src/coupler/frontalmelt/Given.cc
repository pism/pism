// Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019 PISM Authors
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

#include "Given.hh"

#include "pism/util/IceGrid.hh"
#include "pism/util/Time.hh"

#include "pism/coupler/util/options.hh"

namespace pism {
namespace frontalmelt {

Given::Given(IceGrid::ConstPtr g)
  : FrontalMelt(g, nullptr) {

  m_frontal_melt_rate = allocate_frontal_melt_rate(g);
}

Given::~Given() {
  // empty
}

void Given::init_impl(const Geometry &geometry) {
  (void) geometry;

  m_log->message(2,
                 "* Initializing the frontal melt model reading melt rates\n"
                 "  from a file...\n");

  ForcingOptions opt(*m_grid->ctx(), "frontal_melt.given");

  {
    unsigned int buffer_size = m_config->get_double("climate_forcing.buffer_size");
    unsigned int evaluations_per_year = m_config->get_double("climate_forcing.evaluations_per_year");
    bool periodic = opt.period > 0;

    PIO file(m_grid->com, "netcdf3", opt.filename, PISM_READONLY);

    m_frontalmeltrate = IceModelVec2T::ForcingField(m_grid,
                                                    file,
                                                    "frontal_melt_rate",
                                                    "", // no standard name
                                                    buffer_size,
                                                    evaluations_per_year,
                                                    periodic);
  }

  m_frontalmeltrate->set_attrs("climate_forcing", "frontal melt rate", "m s-1", "");
  m_frontalmeltrate->metadata().set_string("glaciological_units", "m year-1");

  m_frontalmeltrate->init(opt.filename, opt.period, opt.reference_time);
}

void Given::bootstrap_impl(const Geometry &geometry) {
  (void) geometry;
}
  
void Given::update_impl(const FrontalMeltInputs &inputs, double t, double dt) {
  (void) inputs;

  m_frontalmeltrate->update(t, dt);

  m_frontalmeltrate->average(t, dt);

  m_frontal_melt_rate->copy_from(*m_frontalmeltrate);
}

MaxTimestep Given::max_timestep_impl(double t) const {
  (void) t;
  // FIXME: get time-step restriction from the input data
  return MaxTimestep("frontal_melt given");
}


const IceModelVec2S& Given::frontal_melt_rate_impl() const {
  return *m_frontal_melt_rate;
}

} // end of namespace frontalmelt
} // end of namespace pism
