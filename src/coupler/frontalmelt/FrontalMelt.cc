/* Copyright (C) 2018, 2019, 2020, 2022 Constantine Khroulev and Andy Aschwanden
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

#include "pism/coupler/FrontalMelt.hh"
#include "pism/util/MaxTimestep.hh"
#include "pism/util/pism_utilities.hh" // combine()
#include "pism/geometry/Geometry.hh"
#include "pism/geometry/part_grid_threshold_thickness.hh"
#include "pism/util/Mask.hh"         // GeometryCalculator

namespace pism {

FrontalMeltInputs::FrontalMeltInputs() {
  geometry = nullptr;

  subglacial_water_flux = nullptr;
}

namespace frontalmelt {

/*!
 * Compute retreat rate corresponding to a given frontal melt rate.
 *
 * This code computes the fraction of the front submerged in ocean water and uses it to
 * scale the provided melt rate.
 */
void FrontalMelt::compute_retreat_rate(const Geometry &geometry,
                                       const IceModelVec2S &frontal_melt_rate,
                                       IceModelVec2S &result) const {

  GeometryCalculator gc(*m_config);

  const IceModelVec2S
    &bed_elevation       = geometry.bed_elevation,
    &surface_elevation   = geometry.ice_surface_elevation,
    &ice_thickness       = geometry.ice_thickness,
    &sea_level_elevation = geometry.sea_level_elevation;
  const auto &cell_type = geometry.cell_type;

  const double
    ice_density = m_config->get_number("constants.ice.density"),
    alpha       = ice_density / m_config->get_number("constants.sea_water.density");

  IceModelVec::AccessList list{&cell_type, &frontal_melt_rate, &sea_level_elevation,
                               &bed_elevation, &surface_elevation, &ice_thickness, &result};

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      if (cell_type.ice_free_ocean(i, j) and cell_type.next_to_ice(i, j)) {
        const double
          bed       = bed_elevation(i, j),
          sea_level = sea_level_elevation(i, j);

        auto H = ice_thickness.star(i, j);
        auto h = surface_elevation.star(i, j);
        auto M = cell_type.star_int(i, j);

        double H_threshold = part_grid_threshold_thickness(M, H, h, bed);

        int m = gc.mask(sea_level, bed, H_threshold);

        double H_submerged = (mask::grounded(m) ?
                              std::max(sea_level - bed, 0.0) :
                              alpha * H_threshold);

        result(i, j) = (H_submerged / H_threshold) * frontal_melt_rate(i, j);
      } else {
        result(i, j) = 0.0;
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();
}

// "modifier" constructor
FrontalMelt::FrontalMelt(IceGrid::ConstPtr g, std::shared_ptr<FrontalMelt> input)
  : Component(g),
    m_input_model(input),
    m_retreat_rate(m_grid, "retreat_rate_due_to_frontal_melt")
{
  m_retreat_rate.set_attrs("diagnostic", "retreat rate due to frontal melt",
                           "m s-1", "m day-1", "", 0);

  m_include_floating_ice = m_config->get_flag("frontal_melt.include_floating_ice");
}

// "model" constructor
FrontalMelt::FrontalMelt(IceGrid::ConstPtr g)
  : FrontalMelt(g, nullptr) {
  // empty
}

void FrontalMelt::init(const Geometry &geometry) {
  this->init_impl(geometry);
}

void FrontalMelt::init_impl(const Geometry &geometry) {
  if (m_input_model) {
    m_input_model->init(geometry);
  }
}

void FrontalMelt::update(const FrontalMeltInputs &inputs, double t, double dt) {
  this->update_impl(inputs, t, dt);

  compute_retreat_rate(*inputs.geometry, frontal_melt_rate(), m_retreat_rate);
}

const IceModelVec2S& FrontalMelt::frontal_melt_rate() const {
  return frontal_melt_rate_impl();
}

const IceModelVec2S& FrontalMelt::retreat_rate() const {
  return m_retreat_rate;
}

// pass-through default implementations for "modifiers"

void FrontalMelt::update_impl(const FrontalMeltInputs &inputs, double t, double dt) {
  if (m_input_model) {
    m_input_model->update(inputs, t, dt);
  } else {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "no input model");
  }
}

MaxTimestep FrontalMelt::max_timestep_impl(double t) const {
  if (m_input_model) {
    return m_input_model->max_timestep(t);
  } else {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "no input model");
  }
}

void FrontalMelt::define_model_state_impl(const File &output) const {
  if (m_input_model) {
    return m_input_model->define_model_state(output);
  } else {
    // no state to define
  }
}

void FrontalMelt::write_model_state_impl(const File &output) const {
  if (m_input_model) {
    return m_input_model->write_model_state(output);
  } else {
    // no state to write
  }
}

namespace diagnostics {

/*! @brief Report frontal melt rate. */
class FrontalMeltRate : public DiagAverageRate<FrontalMelt>
{
public:
  FrontalMeltRate(const FrontalMelt *m)
    : DiagAverageRate<FrontalMelt>(m, "frontal_melt_rate", RATE) {

    m_vars = {SpatialVariableMetadata(m_sys, "frontal_melt_rate")};
    m_accumulator.metadata()["units"] = "m";

    set_attrs("frontal melt rate", "",
              "m second-1", "m day-1", 0);
    m_vars[0]["cell_methods"] = "time: mean";

    m_vars[0]["_FillValue"] = {to_internal(m_fill_value)};
  }

protected:
  const IceModelVec2S& model_input() {
    return model->frontal_melt_rate();
  }
};

/*! @brief Report retreat rate due to frontal melt. */
class FrontalMeltRetreatRate : public DiagAverageRate<FrontalMelt>
{
public:
  FrontalMeltRetreatRate(const FrontalMelt *m)
    : DiagAverageRate<FrontalMelt>(m, "frontal_melt_retreat_rate", RATE) {

    m_vars = {SpatialVariableMetadata(m_sys, "frontal_melt_retreat_rate")};
    m_accumulator.metadata()["units"] = "m";

    set_attrs("retreat rate due to frontal melt", "",
              "m second-1", "m year-1", 0);
    m_vars[0]["cell_methods"] = "time: mean";

    m_vars[0]["_FillValue"] = {to_internal(m_fill_value)};
    m_vars[0]["comment"] = "takes into account what part of the front is submerged";
  }

protected:
  const IceModelVec2S& model_input() {
    return model->retreat_rate();
  }
};

} // end of namespace diagnostics

DiagnosticList FrontalMelt::diagnostics_impl() const {
  using namespace diagnostics;
  DiagnosticList result = {
    {"frontal_melt_rate", Diagnostic::Ptr(new FrontalMeltRate(this))},
    {"frontal_melt_retreat_rate", Diagnostic::Ptr(new FrontalMeltRetreatRate(this))}
  };

  if (m_input_model) {
    return combine(m_input_model->diagnostics(), result);
  } else {
    return result;
  }
}

TSDiagnosticList FrontalMelt::ts_diagnostics_impl() const {
  if (m_input_model) {
    return m_input_model->ts_diagnostics();
  } else {
    return {};
  }
}

bool FrontalMelt::apply(const CellTypeArray0 &M, int i, int j) {
  // icy and grounded_ice cells are included for visualization only (values at these
  // locations have no effect)
  if (m_include_floating_ice) {
    return (M.ice_free_ocean(i, j) and M.next_to_ice(i, j)) or M.icy(i, j);
  } else {
    return (M.ice_free_ocean(i, j) and M.next_to_grounded_ice(i, j)) or M.grounded_ice(i, j);
  }
}


} // end of namespace frontalmelt
} // end of namespace pism
