/* Copyright (C) 2018, 2023 PISM Authors
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
#include "pism/util/Mask.hh"
#include "pism/util/Vars.hh"

#include "pism/util/pism_utilities.hh" // combine

namespace pism {
namespace ocean {
namespace lake_level {

// "Modifier" constructor.
LakeLevel::LakeLevel(IceGrid::ConstPtr grid, std::shared_ptr<LakeLevel> input)
  : Component(grid),
    m_input_model(input),
    m_lake_level(grid, "lake_level", WITHOUT_GHOSTS),
    m_fill_value(m_config->get_number("output.fill_value")) {

  m_lake_level.set_attrs("diagnostic",
                         "lake level elevation, relative to the geoid",
                         "meter", "meter", "", 0);
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
    // do nothing
  }
}

const IceModelVec2S& LakeLevel::elevation() const {
  return m_lake_level;
}

bool LakeLevel::expandMargins() const {
  return expandMargins_impl();
}

bool LakeLevel::expandMargins_impl() const {
  if (m_input_model) {
    return m_input_model->expandMargins();
  }
  return false;
}

MaxTimestep LakeLevel::max_timestep_impl(double t) const {
  if (m_input_model) {
    return m_input_model->max_timestep(t);
  }
  return MaxTimestep("lake level forcing");
}

void LakeLevel::define_model_state_impl(const File &output) const {
  if (m_input_model) {
    m_input_model->define_model_state(output);
  }
}

void LakeLevel::write_model_state_impl(const File &output) const {
  if (m_input_model) {
    m_input_model->write_model_state(output);
  }
}

namespace diagnostics {

/*! @brief Lake level elevation. */
class LakeLevelReal : public Diag<LakeLevel> {
public:
  LakeLevelReal(const LakeLevel *m)
    : Diag<LakeLevel>(m) {
    /* set metadata: */
    m_vars = {SpatialVariableMetadata(m_sys, "lake_level_real")};

    set_attrs("lake level elevation, relative to the geoid", "", "meters", "meters", 0);
    metadata().set_number("_FillValue", m_fill_value);
  }

protected:
  IceModelVec::Ptr compute_impl() const {

    IceModelVec2S::Ptr result(new IceModelVec2S(m_grid, "lake_level", WITHOUT_GHOSTS)),
                       ll(new IceModelVec2S(m_grid, "lake_level", WITHOUT_GHOSTS));
    result->metadata(0) = m_vars[0];
    ll->copy_from(model->elevation());

    result->set(m_fill_value);

    const IceModelVec2S *bed = m_grid->variables().get_2d_scalar("bedrock_altitude"),
                        *thk = m_grid->variables().get_2d_scalar("land_ice_thickness");

    IceModelVec::AccessList list{ bed, thk };
    list.add(*ll);
    list.add(*result);

    GeometryCalculator gc(*m_config);

    ParallelSection ParSec(m_grid->com);
    try {
      for (Points p(*m_grid); p; p.next()) {
        const int i = p.i(), j = p.j();
        const bool isLake = mask::ocean(gc.mask(m_fill_value, (*bed)(i, j), (*thk)(i, j), (*ll)(i, j)));

        if (isLake) {
          (*result)(i, j) = (*ll)(i, j);
        }
      }
    } catch (...) {
      ParSec.failed();
    }
    ParSec.check();

    return result;
  }
};

/*! @brief Lake level elevation. */
class LakeDepth : public Diag<LakeLevel> {
public:
  LakeDepth(const LakeLevel *m)
    : Diag<LakeLevel>(m) {
    /* set metadata: */
    m_vars = {SpatialVariableMetadata(m_sys, "lake_depth")};

    set_attrs("lake depth", "", "meters", "meters", 0);
    metadata().set_number("_FillValue", m_fill_value);
  }

protected:
  IceModelVec::Ptr compute_impl() const {

    IceModelVec2S::Ptr result(new IceModelVec2S(m_grid, "lake_depth", WITHOUT_GHOSTS)),
                       ll(new IceModelVec2S(m_grid, "lake_level", WITHOUT_GHOSTS));
    result->metadata(0) = m_vars[0];
    ll->copy_from(model->elevation());

    result->set(m_fill_value);

    const IceModelVec2S *bed = m_grid->variables().get_2d_scalar("bedrock_altitude"),
                        *thk = m_grid->variables().get_2d_scalar("land_ice_thickness");

    IceModelVec::AccessList list{ bed, thk };
    list.add(*ll);
    list.add(*result);

    GeometryCalculator gc(*m_config);

    ParallelSection ParSec(m_grid->com);
    try {
      for (Points p(*m_grid); p; p.next()) {
        const int i = p.i(), j = p.j();
        const bool isLake = mask::ocean(gc.mask(m_fill_value, (*bed)(i, j), (*thk)(i, j), (*ll)(i, j)));

        if (isLake) {
          (*result)(i, j) = (*ll)(i, j) - (*bed)(i, j);
        }
      }
    } catch (...) {
      ParSec.failed();
    }
    ParSec.check();

    return result;
  }
};


} // end of namespace diagnostics

DiagnosticList LakeLevel::diagnostics_impl() const {
  DiagnosticList result = {
    {"lake_level_real", Diagnostic::Ptr(new diagnostics::LakeLevelReal(this))},
    {"lake_depth", Diagnostic::Ptr(new diagnostics::LakeDepth(this))},
    {"lake_level", Diagnostic::wrap(m_lake_level)}
  };

  if (m_input_model) {
    return combine(result, m_input_model->diagnostics());
  }

  return result;
}

TSDiagnosticList LakeLevel::ts_diagnostics_impl() const {
  if (m_input_model) {
    return m_input_model->ts_diagnostics();
  }

  return {};
}

} // end of namespace lake_level
} // end of namespace ocean
} // end of namespace pism
