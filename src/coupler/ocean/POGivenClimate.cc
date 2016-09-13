// Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016 PISM Authors
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

#include <gsl/gsl_math.h>

#include "POGivenClimate.hh"
#include "base/util/IceGrid.hh"

namespace pism {
namespace ocean {
Given::Given(IceGrid::ConstPtr g)
  : PGivenClimate<OceanModifier,OceanModel>(g, NULL) {

  m_option_prefix   = "-ocean_given";

  // will be de-allocated by the parent's destructor
  m_shelfbtemp     = new IceModelVec2T;
  m_shelfbmassflux = new IceModelVec2T;

  m_fields["shelfbtemp"]     = m_shelfbtemp;
  m_fields["shelfbmassflux"] = m_shelfbmassflux;

  process_options();

  std::map<std::string, std::string> standard_names;
  set_vec_parameters(standard_names);

  m_shelfbtemp->create(m_grid, "shelfbtemp");
  m_shelfbmassflux->create(m_grid, "shelfbmassflux");

  m_shelfbtemp->set_attrs("climate_forcing",
                        "absolute temperature at ice shelf base",
                        "Kelvin", "");
  m_shelfbmassflux->set_attrs("climate_forcing",
                            "ice mass flux from ice shelf base (positive flux is loss from ice shelf)",
                            "kg m-2 s-1", "");
  m_shelfbmassflux->metadata().set_string("glaciological_units", "kg m-2 year-1");
  m_shelfbmassflux->write_in_glaciological_units = true;
}

Given::~Given() {
  // empty
}

void Given::init_impl() {

  m_t = m_dt = GSL_NAN;  // every re-init restarts the clock

  m_log->message(2,
             "* Initializing the ocean model reading base of the shelf temperature\n"
             "  and sub-shelf mass flux from a file...\n");

  m_shelfbtemp->init(m_filename, m_bc_period, m_bc_reference_time);
  m_shelfbmassflux->init(m_filename, m_bc_period, m_bc_reference_time);

  // read time-independent data right away:
  if (m_shelfbtemp->get_n_records() == 1 && m_shelfbmassflux->get_n_records() == 1) {
    update(m_grid->ctx()->time()->current(), 0); // dt is irrelevant
  }
}

void Given::update_impl(double my_t, double my_dt) {
  update_internal(my_t, my_dt);

  m_shelfbmassflux->average(m_t, m_dt);
  m_shelfbtemp->average(m_t, m_dt);
}

void Given::sea_level_elevation_impl(double &result) const {
  result = m_sea_level;
}

void Given::shelf_base_temperature_impl(IceModelVec2S &result) const {
  result.copy_from(*m_shelfbtemp);
}


void Given::shelf_base_mass_flux_impl(IceModelVec2S &result) const {
  result.copy_from(*m_shelfbmassflux);
}

void Given::melange_back_pressure_fraction_impl(IceModelVec2S &result) const {
  result.set(0.0);
}
} // end of namespace ocean
} // end of namespace pism
