// Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019 Ed Bueler and Constantine Khroulev
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

#include "BedThermalUnit.hh"
#include "pism/util/io/PIO.hh"
#include "pism/util/Vars.hh"
#include "pism/util/IceGrid.hh"
#include "pism/util/pism_options.hh"
#include "pism/util/ConfigInterface.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/MaxTimestep.hh"
#include "pism/util/pism_utilities.hh"

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

  InputOptions opts = process_input_options(ctx->com(), ctx->config());

  const Logger &log = *ctx->log();

  if (opts.type == INIT_RESTART) {
    options::ignored(log, "-Mbz");
    options::ignored(log, "-Lbz");

    // If we're initializing from a file we need to get the number of bedrock
    // levels and the depth of the bed thermal layer from it:
    PIO input_file(ctx->com(), "guess_mode", opts.filename, PISM_READONLY);

    if (input_file.inq_var("litho_temp")) {
      grid_info info(input_file, "litho_temp", ctx->unit_system(),
                     CELL_CENTER); // grid registration is irrelevant

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
    options::Integer M("-Mbz", "number of levels in bedrock thermal layer",
                       result.Mbz);

    options::Real L("-Lbz", "depth (thickness) of bedrock thermal layer, in meters",
                    result.Lbz);

    if (M.is_set() and M == 1) {
      options::ignored(log, "-Lbz");
      result.Lbz = 0;
      result.Mbz = 1;
    } else {
      if (M.is_set() ^ L.is_set()) {
        throw RuntimeError(PISM_ERROR_LOCATION, "please specify both -Mbz and -Lbz");
      }

      result.Lbz = L;
      result.Mbz = M;
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
  : Component(g) {

  {
    m_top_surface_flux.create(m_grid, "heat_flux_from_bedrock", WITHOUT_GHOSTS);
    m_top_surface_flux.set_attrs("diagnostic", "upward geothermal flux at the top bedrock surface",
                                 "W m-2",
                                 "upward_geothermal_heat_flux_at_ground_level"); // InitMIP "standard" name
    m_top_surface_flux.metadata().set_string("glaciological_units", "mW m-2");
    m_top_surface_flux.metadata().set_string("comment", "positive values correspond to an upward flux");
  }
  {
    m_bottom_surface_flux.create(m_grid, "bheatflx", WITHOUT_GHOSTS);
    // PROPOSED standard_name = lithosphere_upward_heat_flux
    m_bottom_surface_flux.set_attrs("model_state",
                                    "upward geothermal flux at the bottom bedrock surface",
                                    "W m-2", "");
    m_bottom_surface_flux.metadata().set_string("glaciological_units", "mW m-2");
    m_bottom_surface_flux.metadata().set_string("comment", "positive values correspond to an upward flux");
    m_bottom_surface_flux.set_time_independent(true);
  }
}

BedThermalUnit::~BedThermalUnit() {
  // empty
}

void BedThermalUnit::init(const InputOptions &opts) {
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

DiagnosticList BedThermalUnit::diagnostics_impl() const {
  return {
    {"bheatflx",   Diagnostic::wrap(m_bottom_surface_flux)},
    {"hfgeoubed", Diagnostic::Ptr(new BTU_geothermal_flux_at_ground_level(this))}};
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
  m_vars = {model->flux_through_top_surface().metadata()};
}

IceModelVec::Ptr BTU_geothermal_flux_at_ground_level::compute_impl() const {
  IceModelVec2S::Ptr result(new IceModelVec2S(m_grid, "hfgeoubed", WITHOUT_GHOSTS));
  result->metadata() = m_vars[0];

  result->copy_from(model->flux_through_top_surface());

  return result;
}

} // end of namespace energy
} // end of namespace pism
