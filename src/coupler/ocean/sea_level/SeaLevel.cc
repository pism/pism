/* Copyright (C) 2018, 2019, 2021, 2023 PISM Authors
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

#include "pism/coupler/SeaLevel.hh"

#include "pism/util/MaxTimestep.hh"

#include "pism/util/pism_utilities.hh" // combine

namespace pism {
namespace ocean {
namespace sea_level {

// "Modifier" constructor.
SeaLevel::SeaLevel(std::shared_ptr<const Grid> grid, std::shared_ptr<SeaLevel> input)
  : Component(grid),
    m_input_model(input),
    m_sea_level(grid, "sea_level") {

  m_sea_level.metadata(0)
      .long_name("sea level elevation, relative to the geoid")
      .units("meter");
}

// "Model" constructor (returns sea level is zero).
SeaLevel::SeaLevel(std::shared_ptr<const Grid> g)
  : SeaLevel(g, std::shared_ptr<SeaLevel>()) {
  // empty
}

void SeaLevel::init(const Geometry &geometry) {
  init_impl(geometry);
}

void SeaLevel::init_impl(const Geometry &geometry) {
  if (m_input_model) {
    m_input_model->init(geometry);
  } else {
    double z_s = m_config->get_number("sea_level.constant.value");
    m_sea_level.set(z_s);
    m_log->message(2, "* Using constant sea level at %f meters...\n", z_s);
  }
}

void SeaLevel::update(const Geometry &geometry, double t, double dt) {
  update_impl(geometry, t, dt);
}

void SeaLevel::update_impl(const Geometry &geometry, double t, double dt) {
  if (m_input_model) {
    m_input_model->update(geometry, t, dt);
  } else {
    double z_s = m_config->get_number("sea_level.constant.value");
    m_sea_level.set(z_s);
  }
}

const array::Scalar& SeaLevel::elevation() const {
  return m_sea_level;
}

MaxTimestep SeaLevel::max_timestep_impl(double t) const {
  if (m_input_model) {
    return m_input_model->max_timestep(t);
  }
  return MaxTimestep("sea level forcing");
}

void SeaLevel::define_model_state_impl(const File &output) const {
  if (m_input_model) {
    m_input_model->define_model_state(output);
  }
}

void SeaLevel::write_model_state_impl(const File &output) const {
  if (m_input_model) {
    m_input_model->write_model_state(output);
  }
}

namespace diagnostics {

/*! @brief Sea level elevation. */
class SL : public Diag<SeaLevel> {
public:
  SL(const SeaLevel *m)
    : Diag<SeaLevel>(m) {
    /* set metadata: */
    m_vars = {SpatialVariableMetadata(m_sys, "sea_level")};

    set_attrs("sea level elevation, relative to the geoid", "", "meters", "meters", 0);
  }

protected:
  array::Array::Ptr compute_impl() const {

    array::Scalar::Ptr result(new array::Scalar(m_grid, "sea_level"));
    result->metadata(0) = m_vars[0];

    result->copy_from(model->elevation());

    return result;
  }
};

} // end of namespace diagnostics

DiagnosticList SeaLevel::diagnostics_impl() const {
  DiagnosticList result = {
    {"sea_level", Diagnostic::Ptr(new diagnostics::SL(this))},
  };

  if (m_input_model) {
    return combine(result, m_input_model->diagnostics());
  } else {
    return result;
  }
}

TSDiagnosticList SeaLevel::ts_diagnostics_impl() const {
  if (m_input_model) {
    return m_input_model->ts_diagnostics();
  } else {
    return {};
  }
}

} // end of namespace sea_level
} // end of namespace ocean
} // end of namespace pism
