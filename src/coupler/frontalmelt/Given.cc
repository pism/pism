// Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022 PISM Authors
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

#include "pism/util/Grid.hh"
#include "pism/util/Time.hh"

#include "pism/coupler/util/options.hh"
#include "pism/geometry/Geometry.hh"

namespace pism {
namespace frontalmelt {

Given::Given(std::shared_ptr<const Grid> grid)
  : FrontalMelt(grid, std::shared_ptr<FrontalMelt>()) {

  m_frontal_melt_rate = array::Forcing::Constant(grid, "frontal_melt_rate", 0.0);
}

void Given::init_impl(const Geometry &geometry) {
  (void) geometry;

  m_log->message(2,
                 "* Initializing the frontal melt model reading melt rates\n"
                 "  from a file...\n");

  ForcingOptions opt(*m_grid->ctx(), "frontal_melt.given");

  {
    unsigned int buffer_size = m_config->get_number("input.forcing.buffer_size");

    File file(m_grid->com, opt.filename, io::PISM_NETCDF3, io::PISM_READONLY);

    m_frontal_melt_rate = std::make_shared<array::Forcing>(m_grid,
                                                      file,
                                                      "frontal_melt_rate",
                                                      "", // no standard name
                                                      buffer_size,
                                                      opt.periodic);
  }

  m_frontal_melt_rate->set_attrs("climate_forcing", "frontal melt rate",
                                 "m s-1", "m year-1", "", 0);

  m_frontal_melt_rate->init(opt.filename, opt.periodic);
}

void Given::update_impl(const FrontalMeltInputs &inputs, double t, double dt) {

  const auto &cell_type = inputs.geometry->cell_type;

  // fill m_frontal_melt_rate with values read from an file
  m_frontal_melt_rate->update(t, dt);
  m_frontal_melt_rate->average(t, dt);

  // post-processing: keep values at grounded (or grounded and floating) margins and in
  // the interior, filling the rest with zeros

  array::AccessScope list{&cell_type, m_frontal_melt_rate.get()};

  for (auto p = m_grid->points(); p; p.next()) {
    const int i = p.i(), j = p.j();

    auto R = (*m_frontal_melt_rate)(i, j);

    if (apply(cell_type, i, j)) {
      (*m_frontal_melt_rate)(i, j) = R;
    } else {
      (*m_frontal_melt_rate)(i, j) = 0.0;
    }
  }
}

MaxTimestep Given::max_timestep_impl(double t) const {

  auto dt = m_frontal_melt_rate->max_timestep(t);

  if (dt.finite()) {
    return MaxTimestep(dt.value(), "frontal_melt given");
  } else {
    return MaxTimestep("frontal_melt given");
  }
}


const array::Scalar& Given::frontal_melt_rate_impl() const {
  return *m_frontal_melt_rate;
}

} // end of namespace frontalmelt
} // end of namespace pism
