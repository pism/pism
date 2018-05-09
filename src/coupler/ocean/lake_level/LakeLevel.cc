/* Copyright (C) 2018 PISM Authors
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

#include "pism/coupler/LakeLevel.hh"

#include "pism/util/MaxTimestep.hh"

#include "pism/util/pism_utilities.hh" // combine

namespace pism {
namespace ocean {
namespace lake_level {

// "Modifier" constructor.
LakeLevel::LakeLevel(IceGrid::ConstPtr grid, std::shared_ptr<LakeLevel> input)
  : Component(grid),
    m_input_model(input),
    m_lake_level(grid, "lake_level", WITHOUT_GHOSTS),
    m_fill_value(m_config->get_double("output.fill_value")) {

  m_lake_level.set_attrs("diagnostic",
                         "lake level elevation, relative to the geoid",
                         "meter", "");
}

// "Model" constructor (returns lake level is missing).
LakeLevel::LakeLevel(IceGrid::ConstPtr g)
  : LakeLevel(g, nullptr) {
  // empty
}

LakeLevel::~LakeLevel() {
  // empty
}

void LakeLevel::init(const Geometry &geometry) {
  init_impl(geometry);
}

void LakeLevel::init_impl(const Geometry &geometry) {
  if (m_input_model) {
    m_input_model->init(geometry);
  } else {
    // set the default value
    m_lake_level.set(m_fill_value);
    m_log->message(2, "* Using no lake level model...\n");
  }
}

void LakeLevel::update(const Geometry &geometry, double t, double dt) {
  update_impl(geometry, t, dt);
}

void LakeLevel::update_impl(const Geometry &geometry, double t, double dt) {
  if (m_input_model) {
    m_input_model->update(geometry, t, dt);
  } else {
    m_lake_level.set(m_fill_value);
  }
}

const IceModelVec2S& LakeLevel::elevation() const {
  return m_lake_level;
}

MaxTimestep LakeLevel::max_timestep_impl(double t) const {
  if (m_input_model) {
    return m_input_model->max_timestep(t);
  }
  return MaxTimestep("lake level forcing");
}

void LakeLevel::define_model_state_impl(const PIO &output) const {
  if (m_input_model) {
    m_input_model->define_model_state(output);
  }
}

void LakeLevel::write_model_state_impl(const PIO &output) const {
  if (m_input_model) {
    m_input_model->write_model_state(output);
  }
}

namespace diagnostics {

/*! @brief Lake level elevation. */
class LL : public Diag<LakeLevel> {
public:
  LL(const LakeLevel *m)
    : Diag<LakeLevel>(m) {
    /* set metadata: */
    m_vars = {SpatialVariableMetadata(m_sys, "lake_level")};

    set_attrs("lake level elevation, relative to the geoid", "", "meters", "meters", 0);
    metadata().set_double("_FillValue", m_fill_value);
  }

protected:
  IceModelVec::Ptr compute_impl() const {

    IceModelVec2S::Ptr result(new IceModelVec2S(m_grid, "lake_level", WITHOUT_GHOSTS));
    result->metadata(0) = m_vars[0];

    result->copy_from(model->elevation());

    return result;
  }
};

} // end of namespace diagnostics

DiagnosticList LakeLevel::diagnostics_impl() const {
  DiagnosticList result = {
    {"lake_level", Diagnostic::Ptr(new diagnostics::LL(this))},
  };

  if (m_input_model) {
    return combine(result, m_input_model->diagnostics());
  } else {
    return result;
  }
}

TSDiagnosticList LakeLevel::ts_diagnostics_impl() const {
  if (m_input_model) {
    return m_input_model->ts_diagnostics();
  } else {
    return {};
  }
}

} // end of namespace lake_level
} // end of namespace ocean
} // end of namespace pism
