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

#include <gsl/gsl_math.h>       // GSL_NAN
#include <memory>

#include "pism/coupler/DebrisModel.hh"
#include "pism/util/Time.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/MaxTimestep.hh"
#include "pism/util/Context.hh"

namespace pism {
namespace debris {

std::shared_ptr<array::Scalar> DebrisModel::allocate_debris(std::shared_ptr<const Grid> grid) {
  auto result = std::make_shared<array::Scalar>(grid, "debris_thickness");

  result->metadata(0)
      .long_name("debris thickness")
      .units("m");

  return result;
}


DebrisModel::DebrisModel(std::shared_ptr<const Grid> g)
  : Component(g) {
  // empty
}

DebrisModel::DebrisModel(std::shared_ptr<const Grid> g,
                                 std::shared_ptr<DebrisModel> input)
  :Component(g), m_input_model(input) {
  // empty
}

void DebrisModel::init(const Geometry &geometry) {
  this->init_impl(geometry);
}

Inputs::Inputs()
  : geometry(nullptr),
    u3(nullptr), v3(nullptr), w3(nullptr),
    top_surface_mass_balance(nullptr),
    bottom_surface_mass_balance(nullptr),
    no_model_mask(nullptr) {
  // empty
}

void Inputs::check() const {
  if (geometry == nullptr) {
    throw RuntimeError(PISM_ERROR_LOCATION, "debris::Inputs: geometry is not set");
  }
}

void Inputs::check_transport() const {
  check();

  const char *missing = nullptr;
  if (u3 == nullptr) {
    missing = "u3";
  } else if (v3 == nullptr) {
    missing = "v3";
  } else if (w3 == nullptr) {
    missing = "w3";
  } else if (top_surface_mass_balance == nullptr) {
    missing = "top_surface_mass_balance";
  } else if (bottom_surface_mass_balance == nullptr) {
    missing = "bottom_surface_mass_balance";
  }

  if (missing != nullptr) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "debris::Inputs: '%s' is required by the debris transport model"
                                  " but is not set", missing);
  }
}

void DebrisModel::update(const Inputs &inputs, double t, double dt) {
  inputs.check();
  this->update_impl(inputs, t, dt);
}

const array::Scalar& DebrisModel::debris() const {
  return this->debris_impl();
}

void DebrisModel::begin_pointwise_access() const {
  this->begin_pointwise_access_impl();
}

void DebrisModel::end_pointwise_access() const {
  this->end_pointwise_access_impl();
}

void DebrisModel::init_timeseries(const std::vector<double> &ts) const {
  this->init_timeseries_impl(ts);
}

void DebrisModel::debris_time_series(int i, int j, std::vector<double> &result) const {
  result.resize(m_ts_times.size());
  this->debris_time_series_impl(i, j, result);
}

namespace diagnostics {

/*! @brief Instantaneous debris thickness. */
class DebrisSnapshot : public Diag<DebrisModel> {
public:
  DebrisSnapshot(const DebrisModel *m) : Diag<DebrisModel>(m) {
    m_vars = { { m_sys, "debris_thickness snapshot", *m_grid } };
    m_vars[0].long_name("instantaneous value of the debris thickness").units("m");
  }

protected:
  std::shared_ptr<array::Array> compute_impl(const Geometry &/*geometry*/) const {

    auto result = allocate<array::Scalar>("debris_snapshot");

    std::vector<double> current_time = { m_grid->ctx()->time()->current() };
    std::vector<double> debris  = { 0.0 };

    model->init_timeseries(current_time);

    model->begin_pointwise_access();

    array::AccessScope list(*result);
    ParallelSection loop(m_grid->com);
    try {
      for (auto p : m_grid->points()) {
        const int i = p.i(), j = p.j();

        model->debris_time_series(i, j, debris);

        (*result)(i, j) = debris[0];
      }
    } catch (...) {
      loop.failed();
    }
    loop.check();

    model->end_pointwise_access();

    return result;
  }
};

/*! @brief Effective near-surface mean-annual debris thickness. */
class Debris : public Diag<DebrisModel> {
public:
  Debris(const DebrisModel *m) : Diag<DebrisModel>(m) {
    m_vars = { { m_sys, "effective_debris_thickness", *m_grid } };
    m_vars[0].long_name("effective mean-annual debris thickness").units("m");
  }

protected:
  std::shared_ptr<array::Array> compute_impl(const Geometry &/*geometry*/) const {
    auto result = allocate<array::Scalar>("effective_debris");

    result->copy_from(model->debris());

    return result;
  }
};


} // end of namespace diagnostics

void DebrisModel::update_impl(const Inputs &inputs, double t, double dt) {
  if (m_input_model) {
    m_input_model->update(inputs, t, dt);
  }
}

MaxTimestep DebrisModel::max_timestep_impl(double t, const CFLData *cfl_data) const {
  if (m_input_model) {
    return m_input_model->max_timestep(t, cfl_data);
  }
  return MaxTimestep("debris model");
}

DiagnosticList DebrisModel::spatial_diagnostics_impl() const {
  using namespace diagnostics;

  DiagnosticList result = {
    {"debris_thickness_snapshot",       Diagnostic::Ptr(new DebrisSnapshot(this))},
    {"effective_debris_thickness",      Diagnostic::Ptr(new Debris(this))},
  };

  if (m_input_model) {
    result = combine(result, m_input_model->spatial_diagnostics());
  }

  return result;
}

TSDiagnosticList DebrisModel::scalar_diagnostics_impl() const {
  if (m_input_model) {
    return m_input_model->scalar_diagnostics();
  }

  return {};
}

std::set<VariableMetadata> DebrisModel::state_impl() const {
  if (m_input_model) {
    return m_input_model->state();
  }
  return {};
}

void DebrisModel::write_state_impl(const OutputFile &output) const {
  if (m_input_model) {
    m_input_model->write_state(output);
  }
}

const array::Scalar& DebrisModel::debris_impl() const {
  if (m_input_model) {
    return m_input_model->debris();
  }
  throw RuntimeError::formatted(PISM_ERROR_LOCATION, "no input model");
}


void DebrisModel::begin_pointwise_access_impl() const {
  if (m_input_model) {
    m_input_model->begin_pointwise_access();
  } else {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "no input model");
  }
}

void DebrisModel::end_pointwise_access_impl() const {
  if (m_input_model) {
    m_input_model->end_pointwise_access();
  } else {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "no input model");
  }
}

void DebrisModel::debris_time_series_impl(int i, int j, std::vector<double> &result) const {
  if (m_input_model) {
    m_input_model->debris_time_series(i, j, result);
  } else {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "no input model");
  }
}

void DebrisModel::init_timeseries_impl(const std::vector<double> &ts) const {
  if (m_input_model) {
    m_input_model->init_timeseries(ts);
    m_ts_times = ts;
  } else {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "no input model");
  }
}

} // end of namespace debris
} // end of namespace pism
