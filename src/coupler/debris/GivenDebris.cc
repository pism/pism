// Copyright (C) 2026 Constantine Khroulev and Andy Aschwanden
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

#include "pism/coupler/debris/GivenDebris.hh"

#include "pism/util/Time.hh"
#include "pism/util/Grid.hh"
#include "pism/coupler/util/options.hh"
#include "pism/util/array/Forcing.hh"
#include "pism/util/Logger.hh"
#include "pism/util/io/IO_Flags.hh"

namespace pism {
namespace debris {

Given::Given(std::shared_ptr<const Grid> grid)
  : DebrisModel(grid)
{

  ForcingOptions opt(*m_grid->ctx(), "debris.given");

  {
    unsigned int buffer_size = m_config->get_number("input.forcing.buffer_size");

    File file(m_grid->com, opt.filename, io::PISM_GUESS, io::PISM_READONLY);

    m_debris_thickness = std::make_shared<array::Forcing>(m_grid,
                                                file,
                                                "debris_thickness",
                                                "", // no standard name
                                                buffer_size,
                                                opt.periodic,
                                                LINEAR);
  }

  m_debris_thickness->metadata(0)
      .long_name("debris thickness")
      .units("m");

}

void Given::init_impl(const Geometry &geometry) {

  m_log->message(2,
                 "* Initializing the debris model reading debris thickness\n"
                 "  from a file...\n");

  ForcingOptions opt(*m_grid->ctx(), "debris.given");

  m_debris_thickness->init(opt.filename, opt.periodic);

  // read time-independent data right away:
  if (m_debris_thickness->buffer_size() == 1) {
    Inputs inputs;
    inputs.geometry = &geometry;
    update(inputs, time().current(), 0); // dt is irrelevant
  }
}

void Given::update_impl(const Inputs &inputs, double t, double dt) {
  (void) inputs;

  m_debris_thickness->update(t, dt);

  m_debris_thickness->average(t, dt);

}

const array::Scalar &Given::debris_impl() const {
  return *m_debris_thickness;
}

std::set<VariableMetadata> Given::state_impl() const {
  return array::metadata({ m_debris_thickness.get()});
}

void Given::write_state_impl(const OutputFile &output) const {
  m_debris_thickness->write(output);
}

} // end of namespace debris
} // end of namespace pism
