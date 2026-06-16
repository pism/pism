// Copyright (C) 2026 PISM Authors
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

#include <utility> // std::move

#include "pism/coupler/surface/DEBMEnhanced.hh"

#include "pism/coupler/AtmosphereModel.hh"
#include "pism/geometry/Geometry.hh"
#include "pism/util/Diagnostic.hh"
#include "pism/util/Grid.hh"
#include "pism/util/Logger.hh"
#include "pism/util/Units.hh"
#include "pism/util/array/Array.hh"
#include "pism/util/array/Forcing.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/io/File.hh"
#include "pism/util/io/IO_Flags.hh"

namespace pism {
namespace surface {

///// dEBM-enhanced: dEBM-simple driven by a prescribed daily surface-insolation field.

DEBMEnhanced::DEBMEnhanced(std::shared_ptr<const Grid> g,
                           std::shared_ptr<atmosphere::AtmosphereModel> input)
  : DEBMSimple(g, std::move(input)) {

  // dEBM-enhanced reuses the analytic orbit only for the temperature/offset melt-period
  // weighting. The prescribed insolation is computed for a fixed (present-day) orbit, so it
  // is inconsistent with the paleo orbit parameterization.
  if (m_config->get_flag("surface.debm_simple.paleo.enabled")) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "the 'debm_enhanced' surface model is incompatible with "
                                  "surface.debm_simple.paleo.enabled");
  }

  m_insolation_file = m_config->get_string("surface.debm_enhanced.file");
  if (m_insolation_file.empty()) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "the 'debm_enhanced' surface model requires "
                                  "surface.debm_enhanced.file (create it with "
                                  "util/pism_compute_insolation)");
  }

  m_insolation_periodic = m_config->get_flag("surface.debm_enhanced.periodic");

  m_seconds_per_day = units::convert(m_sys, 1.0, "days", "seconds");

  m_log->message(2, "  Reading daily surface insolation from '%s'...\n",
                 m_insolation_file.c_str());

  int buffer_size = static_cast<int>(m_config->get_number("input.forcing.buffer_size"));

  File file(m_grid->com, m_insolation_file, io::PISM_GUESS, io::PISM_READONLY);

  m_input_insolation = std::make_shared<array::Forcing>(m_grid,
                                                        file,
                                                        "insolation",
                                                        "", // no standard name
                                                        buffer_size,
                                                        m_insolation_periodic,
                                                        LINEAR);
  m_input_insolation->metadata(0)
      .long_name("prescribed daily surface insolation energy")
      .units("J m^-2");
}

void DEBMEnhanced::init_impl(const Geometry &geometry) {
  // initialize the dEBM-simple machinery first
  DEBMSimple::init_impl(geometry);

  m_log->message(2,
                 "* dEBM-enhanced: the insolation-driven melt uses prescribed daily surface\n"
                 "  insolation read from '%s' (instead of the analytic top-of-atmosphere\n"
                 "  parameterization).\n",
                 m_insolation_file.c_str());

  m_input_insolation->init(m_insolation_file, m_insolation_periodic);
}

void DEBMEnhanced::update_insolation_input(double t, double dt,
                                           const std::vector<double> &ts,
                                           array::AccessScope &list) {
  m_input_insolation->update(t, dt);
  m_input_insolation->init_interpolation(ts);
  list.add(*m_input_insolation);
}

void DEBMEnhanced::insolation_energy_series(int i, int j,
                                            const std::vector<DEBMSimpleOrbitalParameters> &orbital,
                                            const std::vector<double> &ts,
                                            double dt_sub,
                                            double latitude,
                                            std::vector<double> &result) const {
  (void)ts;       // sampling times are already set via init_interpolation()
  (void)latitude; // spatial variability is encoded in the prescribed field

  // interpolate the prescribed *daily* insolation energy (J m^-2 day^-1) at this cell to
  // each sub-step time
  m_input_insolation->interp(i, j, result);

  // convert the daily energy to the energy reaching the surface during one sub-step of
  // length dt_sub
  double scale = dt_sub / m_seconds_per_day;
  for (size_t k = 0; k < orbital.size(); ++k) {
    result[k] *= scale;
  }
}

DiagnosticList DEBMEnhanced::spatial_diagnostics_impl() const {
  DiagnosticList result = DEBMSimple::spatial_diagnostics_impl();
  // expose the prescribed insolation actually driving the melt
  result["prescribed_insolation"] = Diagnostic::wrap(*m_input_insolation);
  return result;
}

} // end of namespace surface
} // end of namespace pism
