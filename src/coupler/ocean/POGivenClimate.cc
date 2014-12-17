// Copyright (C) 2011, 2012, 2013, 2014 PISM Authors
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

#include "POGivenClimate.hh"
#include "IceGrid.hh"

namespace pism {

POGiven::POGiven(IceGrid &g)
  : PGivenClimate<POModifier,OceanModel>(g, NULL) {

  option_prefix   = "-ocean_given";

  // will be de-allocated by the parent's destructor
  shelfbtemp     = new IceModelVec2T;
  shelfbmassflux = new IceModelVec2T;

  m_fields["shelfbtemp"]     = shelfbtemp;
  m_fields["shelfbmassflux"] = shelfbmassflux;

  process_options();

  std::map<std::string, std::string> standard_names;
  set_vec_parameters(standard_names);

  shelfbtemp->create(m_grid, "shelfbtemp", false);
  shelfbmassflux->create(m_grid, "shelfbmassflux", false);

  shelfbtemp->set_attrs("climate_forcing",
                        "absolute temperature at ice shelf base",
                        "Kelvin", "");
  shelfbmassflux->set_attrs("climate_forcing",
                            "ice mass flux from ice shelf base (positive flux is loss from ice shelf)",
                            "kg m-2 s-1", "");
  shelfbmassflux->set_glaciological_units("kg m-2 year-1");
  shelfbmassflux->write_in_glaciological_units = true;
}

POGiven::~POGiven() {
  // empty
}

void POGiven::init(Vars &) {

  m_t = m_dt = GSL_NAN;  // every re-init restarts the clock

  verbPrintf(2, m_grid.com,
             "* Initializing the ocean model reading base of the shelf temperature\n"
             "  and sub-shelf mass flux from a file...\n");

  shelfbtemp->init(filename, bc_period, bc_reference_time);
  shelfbmassflux->init(filename, bc_period, bc_reference_time);

  // read time-independent data right away:
  if (shelfbtemp->get_n_records() == 1 && shelfbmassflux->get_n_records() == 1) {
    update(m_grid.time->current(), 0); // dt is irrelevant
  }
}

void POGiven::update(double my_t, double my_dt) {
  update_internal(my_t, my_dt);

  shelfbmassflux->average(m_t, m_dt);
  shelfbtemp->average(m_t, m_dt);
}

void POGiven::sea_level_elevation(double &result) {
  result = sea_level;
}

void POGiven::shelf_base_temperature(IceModelVec2S &result) {
  shelfbtemp->copy_to(result);
}


void POGiven::shelf_base_mass_flux(IceModelVec2S &result) {
  shelfbmassflux->copy_to(result);
}

void POGiven::melange_back_pressure_fraction(IceModelVec2S &result) {
  result.set(0.0);
}

} // end of namespace pism
