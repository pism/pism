/* Copyright (C) 2018, 2019 PISM Authors
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
#include "pism/geometry/Geometry.hh"
#include "pism/util/Vars.hh"

#include "pism/util/MaxTimestep.hh"

#include "pism/util/pism_utilities.hh" // combine
#include "pism/util/pism_options.hh"

namespace pism {
namespace ocean {
namespace sea_level {

// "Modifier" constructor.
SeaLevel::SeaLevel(IceGrid::ConstPtr grid, std::shared_ptr<SeaLevel> input)
  : Component(grid),
    m_input_model(input),
    m_sea_level(grid, "sea_level", WITHOUT_GHOSTS),
    m_fill_value(m_config->get_number("output.fill_value")) {

  m_sea_level.set_attrs("diagnostic",
                        "sea level elevation, relative to the geoid",
                        "meter", "meter", "", 0);

  m_const_sl = 0.0;
}

// "Model" constructor (returns sea level is zero).
SeaLevel::SeaLevel(IceGrid::ConstPtr g)
  : SeaLevel(g, std::shared_ptr<SeaLevel>()) {
  // empty
}

SeaLevel::~SeaLevel() {
  // empty
}

void SeaLevel::init(const Geometry &geometry) {
  init_impl(geometry);
}

void SeaLevel::init_impl(const Geometry &geometry) {
  if (m_input_model) {
    m_input_model->init(geometry);
  } else {

    double const_sl = m_const_sl;
    const_sl = options::Real("-ocean_const_sl",
                             "Constant sea level elevation (m)", const_sl);
    m_const_sl = const_sl;

    // set the sea level
    m_sea_level.set(m_const_sl);
    m_log->message(2, "* Using constant sea level (%gm)...\n", m_const_sl);
  }
}

void SeaLevel::update(const Geometry &geometry, double t, double dt) {
  update_impl(geometry, t, dt);
}

void SeaLevel::update_impl(const Geometry &geometry, double t, double dt) {
  if (m_input_model) {
    m_input_model->update(geometry, t, dt);
  } else {
    m_sea_level.set(m_const_sl);
  }
}

const IceModelVec2S& SeaLevel::elevation() const {
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

bool SeaLevel::expandMargins() const {
  return expandMargins_impl();
}

bool SeaLevel::expandMargins_impl() const {
  if (m_input_model) {
    return m_input_model->expandMargins();
  }
  return false;
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
    metadata().set_number("_FillValue", m_fill_value);
  }

protected:
  IceModelVec::Ptr compute_impl() const {

    IceModelVec2S::Ptr result(new IceModelVec2S(m_grid, "sea_level", WITHOUT_GHOSTS));
    result->metadata(0) = m_vars[0];

    result->copy_from(model->elevation());

    return result;
  }
};

/*! @brief Sea level elevation. */
class OceanDepth : public Diag<SeaLevel> {
public:
  OceanDepth(const SeaLevel *m)
    : Diag<SeaLevel>(m) {
    /* set metadata: */
    m_vars = {SpatialVariableMetadata(m_sys, "ocean_depths")};

    set_attrs("ocean depth", "", "meters", "meters", 0);
    metadata().set_number("_FillValue", m_fill_value);
  }

protected:
  IceModelVec::Ptr compute_impl() const {

    IceModelVec2S::Ptr result(new IceModelVec2S(m_grid, "ocean_depth", WITHOUT_GHOSTS));
    result->metadata(0) = m_vars[0];

    const IceModelVec2S &sl  = model->elevation(),
                        &bed = *m_grid->variables().get_2d_scalar("bedrock_altitude"),
                        &thk = *m_grid->variables().get_2d_scalar("land_ice_thickness");
    IceModelVec::AccessList list{ &sl, &bed, &thk };
    list.add(*result);

    GeometryCalculator gc(*m_config);

    for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (mask::ocean(gc.mask(sl(i, j), bed(i,j), thk(i, j)))) {
      (*result)(i, j) = sl(i, j) - bed(i, j);
    } else {
      (*result)(i, j) = m_fill_value;
    }
  } // end loop grid

    return result;
  }
};

} // end of namespace diagnostics

DiagnosticList SeaLevel::diagnostics_impl() const {
  DiagnosticList result = {
    {"sea_level", Diagnostic::Ptr(new diagnostics::SL(this))},
    {"ocean_depth", Diagnostic::Ptr(new diagnostics::OceanDepth(this))},
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
