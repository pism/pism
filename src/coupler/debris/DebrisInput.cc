// Copyright (C) 2026 Andy Aschwanden and Constantine Khroulev
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

#include "pism/coupler/debris/DebrisInput.hh"

#include "pism/util/Config.hh"
#include "pism/util/Context.hh"
#include "pism/util/Grid.hh"
#include "pism/util/Logger.hh"
#include "pism/util/array/Forcing.hh"
#include "pism/util/io/File.hh"
#include "pism/util/io/IO_Flags.hh"

namespace pism {
namespace debris {

DebrisInput::DebrisInput(std::shared_ptr<const Grid> grid)
  : m_grid(grid), m_zero(grid, "debris_input_rate") {

  const Config &config = *grid->ctx()->config();

  m_filename = config.get_string("debris.transport.input.file");
  m_periodic = config.get_flag("debris.transport.input.periodic");

  m_zero.metadata(0)
      .long_name("rate of debris input from the surrounding terrain (solid debris thickness)")
      .units("m s^-1")
      .output_units("m year^-1");
  m_zero.set(0.0);

  if (not m_filename.empty()) {
    unsigned int buffer_size = config.get_number("input.forcing.buffer_size");

    File file(m_grid->com, m_filename, io::PISM_GUESS, io::PISM_READONLY);

    // a rate: piecewise constant over the time bounds of the records (a step-like input,
    // e.g. a source switched on at some time, is represented exactly)
    m_rate = std::make_shared<array::Forcing>(m_grid, file, "debris_input_rate",
                                              "", // no standard name
                                              buffer_size, m_periodic, PIECEWISE_CONSTANT);
    m_rate->metadata(0)
        .long_name("rate of debris input from the surrounding terrain (solid debris thickness)")
        .units("m s^-1")
        .output_units("m year^-1");
  }
}

void DebrisInput::init() {
  const Logger &log = *m_grid->ctx()->log();

  if (m_rate) {
    log.message(2, "  - Reading the debris input rate from '%s'...\n", m_filename.c_str());
    m_rate->init(m_filename, m_periodic);
  } else {
    log.message(2, "  - No debris input file given (debris.transport.input.file):"
                   " the debris input rate is zero.\n");
  }
}

void DebrisInput::update(double t, double dt) {
  if (m_rate) {
    m_rate->update(t, dt);
    m_rate->average(t, dt);
  }
}

const array::Scalar &DebrisInput::rate() const {
  if (m_rate) {
    return *m_rate;
  }
  return m_zero;
}

bool DebrisInput::enabled() const {
  return (bool)m_rate;
}

} // end of namespace debris
} // end of namespace pism
