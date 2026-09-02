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

#include <algorithm>            // std::max
#include <cmath>                // std::pow

#include "pism/coupler/debris/IceMeltEnhancement.hh"

#include "pism/coupler/DebrisModel.hh"
#include "pism/coupler/util/options.hh"
#include "pism/util/Config.hh"
#include "pism/util/Context.hh"
#include "pism/util/Grid.hh"
#include "pism/util/Logger.hh"
#include "pism/util/MaxTimestep.hh"
#include "pism/util/Time.hh"
#include "pism/util/array/Forcing.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/io/IO_Flags.hh"

namespace pism {
namespace debris {

static IceMeltEnhancement::Model model_from_config(const Config &config) {
  auto model = config.get_string("debris.ice_melt_enhancement.model");

  if (model == "none") {
    return IceMeltEnhancement::NONE;
  }

  if (model == "given") {
    return IceMeltEnhancement::GIVEN;
  }

  if (model == "verhaegen") {
    return IceMeltEnhancement::VERHAEGEN;
  }

  throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                "invalid debris.ice_melt_enhancement.model: '%s'"
                                " (allowed choices: %s)",
                                model.c_str(),
                                config.choices("debris.ice_melt_enhancement.model").c_str());
}

IceMeltEnhancement::IceMeltEnhancement(std::shared_ptr<const Grid> grid)
  : IceMeltEnhancement(grid, nullptr) {
  // empty
}

IceMeltEnhancement::IceMeltEnhancement(std::shared_ptr<const Grid> grid,
                                       std::shared_ptr<DebrisModel> debris_model)
  : Component(grid),
    m_model(model_from_config(*grid->ctx()->config())),
    m_ice_melt_enhancement(grid, "ice_melt_enhancement"),
    m_debris_model(debris_model) {

  m_ice_melt_enhancement.metadata(0)
      .long_name("ratio of the sub-debris ice melt rate to the clean ice melt rate")
      .units("1");

  // the "none" model never updates this field, so this is the value used for the whole
  // run
  m_ice_melt_enhancement.set(1.0);

  switch (m_model) {
  case GIVEN:
    {
      ForcingOptions opt(*m_grid->ctx(), "debris.ice_melt_enhancement");

      unsigned int buffer_size = m_config->get_number("input.forcing.buffer_size");

      File file(m_grid->com, opt.filename, io::PISM_GUESS, io::PISM_READONLY);

      m_melt_factor = std::make_shared<array::Forcing>(m_grid,
                                                       file,
                                                       "debris_melt_factor",
                                                       "", // no standard name
                                                       buffer_size,
                                                       opt.periodic,
                                                       LINEAR);

      m_melt_factor->metadata(0)
          .long_name("sub-debris melt enhancement factor")
          .units("1");
      break;
    }
  case VERHAEGEN:
    if (not m_debris_model) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                    "debris.ice_melt_enhancement.model = \"verhaegen\""
                                    " requires a debris model supplying the debris thickness");
    }
    break;
  case NONE:
    break;
  }
}

void IceMeltEnhancement::init(const Geometry &geometry) {

  switch (m_model) {
  case NONE:
    m_log->message(2, "* Debris-related ice melt enhancement is disabled.\n");
    break;
  case GIVEN:
    {
      m_log->message(2,
                     "* Initializing the debris ice melt enhancement model reading\n"
                     "  'debris_melt_factor' from a file...\n");

      ForcingOptions opt(*m_grid->ctx(), "debris.ice_melt_enhancement");

      m_melt_factor->init(opt.filename, opt.periodic);

      // read time-independent data right away:
      if (m_melt_factor->buffer_size() == 1) {
        update(geometry, time().current(), 0); // dt is irrelevant
      }
      break;
    }
  case VERHAEGEN:
    m_log->message(2,
                   "* Initializing the debris ice melt enhancement model using\n"
                   "  equation (14) of Verhaegen and Huybrechts (2026)...\n");

    m_debris_model->init(geometry);
    // the debris thickness may be time-independent, so compute the factor right away:
    update(geometry, time().current(), 0);
    break;
  }
}

void IceMeltEnhancement::update(const Geometry &geometry, double t, double dt) {

  switch (m_model) {
  case NONE:
    // m_ice_melt_enhancement is 1.0 everywhere: nothing to do
    break;
  case GIVEN:
    m_melt_factor->update(t, dt);
    m_melt_factor->average(t, dt);

    m_ice_melt_enhancement.copy_from(*m_melt_factor);
    break;
  case VERHAEGEN:
    {
      m_debris_model->update(geometry, t, dt);

      const array::Scalar &debris_thickness = m_debris_model->debris();

      array::AccessScope list{ &debris_thickness, &m_ice_melt_enhancement };

      for (auto p : m_grid->points()) {
        const int i = p.i(), j = p.j();

        m_ice_melt_enhancement(i, j) = verhaegen_melt_factor(debris_thickness(i, j));
      }
      break;
    }
  }
}

/*!
 * Equation (14) in Verhaegen and Huybrechts (2026),
 * https://doi.org/10.1029/2025JF008748.
 *
 * An empirical fit to the "Ostrem curve" averaged over the 21 debris-covered glacier
 * studies in their Table 1: melt is *enhanced* by a thin, dispersed cover, reaching a
 * maximum of about 1.4 times the clean-ice rate at the effective debris thickness
 * `h_e = 1.5` cm, drops back to the clean-ice rate near the critical thickness
 * `h_c = 4.1` cm, and is suppressed by thicker debris, becoming negligible beyond
 * 1.5 m (Popovnin and Rozova, 2002).
 *
 * @param[in] debris_thickness debris thickness, meters
 * @return the ratio of the sub-debris melt rate to the clean ice melt rate
 */
double IceMeltEnhancement::verhaegen_melt_factor(double debris_thickness) {

  // effective debris thickness (meters): the factor peaks here
  const double h_e = 0.015;
  // critical debris thickness (meters): the factor crosses 1 here
  const double h_c = 0.041;
  // lower bound corresponding to "no significant melt"
  const double f_min = 1e-3;

  const double h = debris_thickness;

  // Clean ice. Written this way so that NaNs (and negative values produced by
  // interpolation) are treated as "no debris" instead of poisoning the melt rate.
  if (not (h > 0.0)) {
    return 1.0;
  }

  if (h <= h_e) {
    return 26.667 * h + 1.0;
  }

  if (h <= h_c) {
    return -16.0 * h + 1.64;
  }

  return std::max(f_min, 0.1061 * std::pow(h, -0.7205) - 0.07922);
}

const array::Scalar &IceMeltEnhancement::ice_melt_enhancement() const {
  return m_ice_melt_enhancement;
}

IceMeltEnhancement::Model IceMeltEnhancement::model() const {
  return m_model;
}

MaxTimestep IceMeltEnhancement::max_timestep_impl(double t, const CFLData *cfl_data) const {
  switch (m_model) {
  case GIVEN:
    return m_melt_factor->max_timestep(t);
  case VERHAEGEN:
    return m_debris_model->max_timestep(t, cfl_data);
  case NONE:
  default:
    return MaxTimestep("debris ice melt enhancement");
  }
}

namespace diagnostics {

/*! @brief Debris-related ice melt enhancement factor. */
class MeltEnhancement : public Diag<IceMeltEnhancement> {
public:
  MeltEnhancement(const IceMeltEnhancement *m) : Diag<IceMeltEnhancement>(m) {
    m_vars = { { m_sys, "ice_melt_enhancement", *m_grid } };
    m_vars[0]
        .long_name("ratio of the sub-debris ice melt rate to the clean ice melt rate")
        .units("1");
  }

protected:
  std::shared_ptr<array::Array> compute_impl(const Geometry & /*geometry*/) const {
    auto result = allocate<array::Scalar>("ice_melt_enhancement");

    result->copy_from(model->ice_melt_enhancement());

    return result;
  }
};

} // end of namespace diagnostics

DiagnosticList IceMeltEnhancement::spatial_diagnostics_impl() const {
  using namespace diagnostics;

  DiagnosticList result = {
    { "ice_melt_enhancement", Diagnostic::Ptr(new MeltEnhancement(this)) },
  };

  if (m_debris_model) {
    result = combine(result, m_debris_model->spatial_diagnostics());
  }

  return result;
}

} // end of namespace debris
} // end of namespace pism
