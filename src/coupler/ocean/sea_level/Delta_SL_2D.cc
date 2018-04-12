// Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018 PISM Authors
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

#include "Delta_SL_2D.hh"

#include "pism/util/IceGrid.hh"
#include "pism/coupler/util/options.hh"

namespace pism {
namespace ocean {
namespace sea_level {

Delta_SL_2D::Delta_SL_2D(IceGrid::ConstPtr grid, std::shared_ptr<SeaLevel> in)
  : SeaLevel(grid, in) {

  ForcingOptions opt(*m_grid->ctx(), "ocean.delta_sl_2d");

  {
    unsigned int buffer_size = m_config->get_double("climate_forcing.buffer_size");
    unsigned int evaluations_per_year = m_config->get_double("climate_forcing.evaluations_per_year");
    bool periodic = opt.period > 0;

    PIO file(m_grid->com, "netcdf3", opt.filename, PISM_READONLY);


    m_forcing = IceModelVec2T::ForcingField(m_grid,
                                            file,
                                            "delta_SL",
                                            "", // no standard name
                                            buffer_size,
                                            evaluations_per_year,
                                            periodic);
    m_forcing->set_attrs("climate_forcing",
                         "two-dimensional sea level offsets",
                         "meters", "");
  }
}

Delta_SL_2D::~Delta_SL_2D() {
  // empty
}

void Delta_SL_2D::init_impl(const Geometry &geometry) {

  m_input_model->init(geometry);

  ForcingOptions opt(*m_grid->ctx(), "ocean.delta_sl_2d");

  m_log->message(2, "* Initializing 2D sea level forcing...\n");
  m_log->message(2, "    reading anomalies from %s ...\n", opt.filename.c_str());

  m_forcing->init(opt.filename, opt.period, opt.reference_time);
}

void Delta_SL_2D::update_impl(const Geometry &geometry, double t, double dt) {
  m_input_model->update(geometry, t, dt);

  m_forcing->update(t, dt);
  m_forcing->average(t, dt);

  m_input_model->elevation().add(1.0, *m_forcing, m_sea_level);
}

} // end of namespace sea_level
} // end of namespace ocean
} // end of namespace pism
