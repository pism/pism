// Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016 Ed Bueler and Constantine Khroulev
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

#include <gsl/gsl_math.h>

#include "BedThermalUnit.hh"
#include "base/util/io/PIO.hh"
#include "base/util/PISMVars.hh"
#include "base/util/IceGrid.hh"
#include "base/util/pism_options.hh"
#include "base/util/PISMConfigInterface.hh"
#include "base/util/error_handling.hh"
#include "base/util/MaxTimestep.hh"
#include "base/util/pism_utilities.hh"

#include "BTU_Full.hh"
#include "BTU_Minimal.hh"

namespace pism {
namespace energy {

BTUGrid::BTUGrid(Context::ConstPtr ctx) {
  Mbz = (unsigned int) ctx->config()->get_double("grid.Mbz");
  Lbz = ctx->config()->get_double("grid.Lbz");
}


BTUGrid BTUGrid::FromOptions(Context::ConstPtr ctx) {
  BTUGrid result(ctx);

  InputOptions opts = process_input_options(ctx->com());

  const Logger &log = *ctx->log();

  if (opts.type == INIT_RESTART) {
    options::ignored(log, "-Mbz");
    options::ignored(log, "-Lbz");

    // If we're initializing from a file we need to get the number of bedrock
    // levels and the depth of the bed thermal layer from it:
    PIO input_file(ctx->com(), "guess_mode", opts.filename, PISM_READONLY);

    if (input_file.inq_var("litho_temp")) {
      grid_info info(input_file, "litho_temp", ctx->unit_system(),
                     NOT_PERIODIC); // periodicity is irrelevant

      result.Mbz = info.z_len;
      result.Lbz = -info.z_min;
    } else {
      // override values we got using config.get_double() in the constructor
      result.Mbz = 1;
      result.Lbz = 0;
    }

    input_file.close();
  } else {
    // Bootstrapping or initializing without an input file.
    options::Integer Mbz("-Mbz", "number of levels in bedrock thermal layer",
                         result.Mbz);

    options::Real Lbz("-Lbz", "depth (thickness) of bedrock thermal layer, in meters",
                      result.Lbz);

    if (Mbz.is_set() and Mbz == 1) {
      options::ignored(log, "-Lbz");
      result.Lbz = 0;
      result.Mbz = 1;
    } else {
      if (Mbz.is_set() ^ Lbz.is_set()) {
        throw RuntimeError(PISM_ERROR_LOCATION, "please specify both -Mbz and -Lbz");
      }

      result.Lbz = Lbz;
      result.Mbz = Mbz;
    }
  }

  return result;
}

/*! Allocate a complete or minimal bedrock thermal unit depending on the number of bedrock levels.
 *
 */
BedThermalUnit* BedThermalUnit::FromOptions(IceGrid::ConstPtr grid,
                                            Context::ConstPtr ctx) {

  BTUGrid bedrock_grid = BTUGrid::FromOptions(ctx);

  if (bedrock_grid.Mbz > 1) {
    return new energy::BTU_Full(grid, bedrock_grid);
  } else {
    return new energy::BTU_Minimal(grid);
  }
}


BedThermalUnit::BedThermalUnit(IceGrid::ConstPtr g)
  : Component_TS(g) {

  {
    m_top_surface_flux.create(m_grid, "heat_flux_from_bedrock", WITHOUT_GHOSTS);
    m_top_surface_flux.set_attrs("diagnostic", "upward geothermal flux at the top bedrock surface",
                                 "W m-2", "");
    m_top_surface_flux.metadata().set_string("glaciological_units", "mW m-2");
    m_top_surface_flux.metadata().set_string("comment", "positive values correspond to an upward flux");
    m_top_surface_flux.write_in_glaciological_units = true;
  }
  {
    m_bottom_surface_flux.create(m_grid, "bheatflx", WITHOUT_GHOSTS);
    // PROPOSED standard_name = lithosphere_upward_heat_flux
    m_bottom_surface_flux.set_attrs("model_state",
                                    "upward geothermal flux at the bottom bedrock surface",
                                    "W m-2", "");
    m_bottom_surface_flux.metadata().set_string("glaciological_units", "mW m-2");
    m_bottom_surface_flux.metadata().set_string("comment", "positive values correspond to an upward flux");
    m_bottom_surface_flux.write_in_glaciological_units = true;
    m_bottom_surface_flux.set_time_independent(true);
  }
}

BedThermalUnit::~BedThermalUnit() {
  // empty
}

void BedThermalUnit::init(const InputOptions &opts) {
  m_t = m_dt = GSL_NAN;  // every re-init restarts the clock

  this->init_impl(opts);
}

//! \brief Initialize the bedrock thermal unit.
void BedThermalUnit::init_impl(const InputOptions &opts) {
  switch (opts.type) {
  case INIT_RESTART:
    m_bottom_surface_flux.read(opts.filename, opts.record);
    break;
  case INIT_BOOTSTRAP:
    m_bottom_surface_flux.regrid(opts.filename, OPTIONAL,
                                 m_config->get_double("bootstrapping.defaults.geothermal_flux"));
    break;
  case INIT_OTHER:
  default:
    initialize_bottom_surface_flux();
  }

  regrid("bedrock thermal layer", m_bottom_surface_flux, REGRID_WITHOUT_REGRID_VARS);
}

void BedThermalUnit::initialize_bottom_surface_flux() {
  const double heat_flux = m_config->get_double("bootstrapping.defaults.geothermal_flux");

  m_log->message(2, "  using constant geothermal flux %f W m-2 ...\n",
                 heat_flux);

  m_bottom_surface_flux.set(heat_flux);
}

/** Returns the vertical spacing used by the bedrock grid.
 *
 */
double BedThermalUnit::vertical_spacing() const {
  return this->vertical_spacing_impl();
}

/*!
 * Return the depth of the bedrock thermal layer.
 */
double BedThermalUnit::depth() const {
  return this->depth_impl();
}

/*!
 * Return the number of levels in the bedrock thermal layer.
 */
unsigned int BedThermalUnit::Mz() const {
  return this->Mz_impl();
}

void BedThermalUnit::define_model_state_impl(const PIO &output) const {
  m_bottom_surface_flux.define(output);
}

void BedThermalUnit::write_model_state_impl(const PIO &output) const {
  m_bottom_surface_flux.write(output);
}

std::map<std::string, Diagnostic::Ptr> BedThermalUnit::diagnostics_impl() const {
  return {{"hfgeoubed", Diagnostic::Ptr(new BTU_geothermal_flux_at_ground_level(this))}};
}

void BedThermalUnit::update_impl(double t, double dt) {
  // CHECK: has the desired time-interval already been dealt with?
  if ((fabs(t - m_t) < 1e-12) and (fabs(dt - m_dt) < 1e-12)) {
    return;
  }

  const IceModelVec2S
    &bedtoptemp = *m_grid->variables().get_2d_scalar("bedtoptemp");

  this->update_impl(bedtoptemp, t, dt);
}

void BedThermalUnit::update(const IceModelVec2S &bedrock_top_temperature,
                            double t, double dt) {
  this->update_impl(bedrock_top_temperature, t, dt);
}

const IceModelVec2S& BedThermalUnit::flux_through_top_surface() const {
  return m_top_surface_flux;
}

const IceModelVec2S& BedThermalUnit::flux_through_bottom_surface() const {
  return m_bottom_surface_flux;
}

BTU_geothermal_flux_at_ground_level::BTU_geothermal_flux_at_ground_level(const BedThermalUnit *m)
  : Diag<BedThermalUnit>(m) {
  m_vars = {SpatialVariableMetadata(m_sys, "hfgeoubed")};
  set_attrs("upward geothermal flux at ground (top of the bedrock) level",
            "",                 // no standard name
            "W m-2", "W m-2", 0);
}

IceModelVec::Ptr BTU_geothermal_flux_at_ground_level::compute_impl() {
  IceModelVec2S::Ptr result(new IceModelVec2S);
  result->create(m_grid, "hfgeoubed", WITHOUT_GHOSTS);
  result->metadata() = m_vars[0];

  result->copy_from(model->flux_through_top_surface());

  return result;
}

} // end of namespace energy
} // end of namespace pism
