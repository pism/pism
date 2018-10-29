/* Copyright (C) 2013, 2014, 2015, 2016, 2017, 2018 PISM Authors
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

#include <gsl/gsl_math.h>       // GSL_NAN

#include "pism/coupler/FrontalMeltModel.hh"
#include "pism/util/iceModelVec.hh"
#include "pism/util/MaxTimestep.hh"
#include "pism/util/pism_utilities.hh"

namespace pism {

FrontalMeltInputs::FrontalMeltInputs() {
  geometry = nullptr;

  subglacial_discharge_at_grounding_line = nullptr;

}

  
namespace frontalmelt {

IceModelVec2S::Ptr FrontalMeltModel::allocate_frontal_melt_rate(IceGrid::ConstPtr g) {
  IceModelVec2S::Ptr result(new IceModelVec2S(g, "frontalmeltrate", WITHOUT_GHOSTS));

  result->set_attrs("diagnostic", "frontal melt rate",
                    "m s-1", "");
  result->metadata().set_string("glaciological_units", "m year-1");

  return result;
}


// "modifier" constructor
FrontalMeltModel::FrontalMeltModel(IceGrid::ConstPtr g, std::shared_ptr<FrontalMeltModel> input)
  : Component(g), m_input_model(input) {

  if (not input) {
  }
}

// "model" constructor
FrontalMeltModel::FrontalMeltModel(IceGrid::ConstPtr g)
  : FrontalMeltModel(g, nullptr) {
  // empty
}


FrontalMeltModel::~FrontalMeltModel() {
  // empty
}

void FrontalMeltModel::init(const Geometry &geometry) {
  this->init_impl(geometry);
}

void FrontalMeltModel::bootstrap(const Geometry &geometry) {
  this->bootstrap_impl(geometry);
}

void FrontalMeltModel::init_impl(const Geometry &geometry) {
  if (m_input_model) {
    m_input_model->init(geometry);
  }
}

void FrontalMeltModel::bootstrap_impl(const Geometry &geometry) {
  if (m_input_model) {
    m_input_model->bootstrap(geometry);
  }
}

  
void FrontalMeltModel::update(const Geometry &geometry, double t, double dt) {
  this->update_impl(geometry, t, dt);
}


const IceModelVec2S& FrontalMeltModel::frontal_melt_rate() const {
  return frontal_melt_rate_impl();
}

// pass-through default implementations for "modifiers"

void FrontalMeltModel::update_impl(const Geometry &geometry, double t, double dt) {
  if (m_input_model) {
    m_input_model->update(geometry, t, dt);
  } else {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "no input model");
  }
}

MaxTimestep FrontalMeltModel::max_timestep_impl(double t) const {
  if (m_input_model) {
    return m_input_model->max_timestep(t);
  } else {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "no input model");
  }
}

void FrontalMeltModel::define_model_state_impl(const PIO &output) const {
  if (m_input_model) {
    return m_input_model->define_model_state(output);
  } else {
    // no state to define
  }
}

void FrontalMeltModel::write_model_state_impl(const PIO &output) const {
  if (m_input_model) {
    return m_input_model->write_model_state(output);
  } else {
    // no state to write
  }
}

const IceModelVec2S& FrontalMeltModel::frontal_melt_rate_impl() const {
  if (m_input_model) {
    return m_input_model->frontal_melt_rate();
  } else {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "no input model");
  }
}

  
namespace diagnostics {


/*! @brief Frontal melt rate. */
class PFM_frontal_melt_rate : public Diag<FrontalMeltModel>
{
public:
  PFM_frontal_melt_rate(const FrontalMeltModel *m)
    : Diag<FrontalMeltModel>(m) {

    /* set metadata: */
    m_vars = {SpatialVariableMetadata(m_sys, "frontal_melt_rate")};

    set_attrs("frontal melt rate", "",
              "m s-1", "m year-1", 0);
  }
protected:
  IceModelVec::Ptr compute_impl() const {

    IceModelVec2S::Ptr result(new IceModelVec2S(m_grid, "frontalmeltrate", WITHOUT_GHOSTS));
    result->metadata(0) = m_vars[0];

    result->copy_from(model->frontal_melt_rate());

    return result;
  }
};


} // end of namespace diagnostics

DiagnosticList FrontalMeltModel::diagnostics_impl() const {
  using namespace diagnostics;
  DiagnosticList result = {
    {"frontalmeltrate",                 Diagnostic::Ptr(new PFM_frontal_melt_rate(this))},
  };

  if (m_input_model) {
    return combine(m_input_model->diagnostics(), result);
  } else {
    return result;
  }
}

TSDiagnosticList FrontalMeltModel::ts_diagnostics_impl() const {
  if (m_input_model) {
    return m_input_model->ts_diagnostics();
  } else {
    return {};
  }
}


} // end of namespace frontalmelt
} // end of namespace pism
