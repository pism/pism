/* Copyright (C) 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2021, 2022 PISM Authors
 *
 * This file is part of PISM.
 *
 * PISM is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option) any later
 * version.
 *
 * PISM is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PISM; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include "GivenRate.hh"

#include "pism/util/Mask.hh"
#include "pism/util/IceGrid.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/coupler/util/options.hh"
#include "pism/util/iceModelVec2T.hh"
#include "pism/util/io/File.hh"

namespace pism {

//! @brief Calving and iceberg removal code.
namespace calving {

GivenRate::GivenRate(IceGrid::ConstPtr g)
  : Component(g) {

  ForcingOptions opt(*m_grid->ctx(), "calving.given_calving");
  {
    unsigned int buffer_size = m_config->get_number("input.forcing.buffer_size");

    File file(m_grid->com, opt.filename, PISM_NETCDF3, PISM_READONLY);

    m_calving_rate = IceModelVec2T::ForcingField(m_grid,
                                                      file,
                                                      "calving_rate",
                                                      "", // no standard name
                                                      buffer_size,
                                                      opt.periodic,
                                                      LINEAR);

    m_calving_rate->set_attrs("diagnostic",
                                   "'calving rate' as used in given_calving method",
                                   "m year-1", "m year-1",
                                   "", 0); // no standard name
    m_calving_rate->metadata()["valid_min"] = {0.0};
  }
}

void GivenRate::init() {

  m_log->message(2, "* Initializing the 'given calving rate' mechanism...\n");

  ForcingOptions opt(*m_grid->ctx(), "calving.thickness_calving");

  File file(m_grid->com, opt.filename, PISM_GUESS, PISM_READONLY);
  auto variable_exists = file.find_variable(m_calving_rate->get_name());
  if (variable_exists) {
    m_log->message(2,
                   "  Reading calving rate from file '%s'...\n",
                   opt.filename.c_str());

    m_calving_rate->init(opt.filename, opt.periodic);
  } else {
    double calving_rate = m_config->get_number("calving.given_calving.rate");

    SpatialVariableMetadata attributes = m_calving_rate->metadata();
    // replace with a constant IceModelVec2T
    m_calving_rate = IceModelVec2T::Constant(m_grid,
                                                 "given_calving_rate",
                                                  calving_rate);
    // restore metadata
    m_calving_rate->metadata() = attributes;

    m_log->message(2,
                   "  Calving rate: %3.3f meters year-1.\n", calving_rate);
  }
}

/**
 * Updates calving rate according to given calving rate
 *
 * @param[in] t beginning of the time step
 * @param[in] dt length of the time step
 *
 * @return 0 on success
 */
void GivenRate::update(double t,
                       double dt) {

  m_calving_rate->update(t, dt);
  m_calving_rate->average(t, dt);

}

const IceModelVec2S& GivenRate::calvingrate() const {
  return *m_calving_rate;
}

DiagnosticList GivenRate::diagnostics_impl() const {
  return {{"given_calving_rate", Diagnostic::wrap(*m_calving_rate)}};
}

} // end of namespace calving
} // end of namespace pism
