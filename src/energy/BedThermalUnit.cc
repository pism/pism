// Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2023 Ed Bueler and Constantine Khroulev
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

#include <memory>

#include "pism/energy/BedThermalUnit.hh"
#include "pism/util/ConfigInterface.hh"
#include "pism/util/Grid.hh"
#include "pism/util/io/File.hh"

#include "pism/energy/BTU_Full.hh"
#include "pism/energy/BTU_Minimal.hh"

namespace pism {
namespace energy {

BTUGrid::BTUGrid(std::shared_ptr<const Context> ctx) {
  Mbz = (unsigned int) ctx->config()->get_number("grid.Mbz");
  Lbz = ctx->config()->get_number("grid.Lbz");
}


BTUGrid BTUGrid::FromOptions(std::shared_ptr<const Context> ctx) {
  BTUGrid result(ctx);

  Config::ConstPtr config = ctx->config();
  InputOptions opts = process_input_options(ctx->com(), config);

  if (opts.type == INIT_RESTART) {
    // If we're initializing from a file we need to get the number of bedrock
    // levels and the depth of the bed thermal layer from it:
    File input_file(ctx->com(), opts.filename, io::PISM_NETCDF3, io::PISM_READONLY);

    if (input_file.find_variable("litho_temp")) {
      grid::InputGridInfo info(input_file, "litho_temp", ctx->unit_system(),
                               grid::CELL_CENTER); // grid registration is irrelevant

      result.Mbz = info.z_len;
      result.Lbz = -info.z_min;
    } else {
      // override values we got using config.get_number() in the constructor
      result.Mbz = 1;
      result.Lbz = 0;
    }

    input_file.close();
  } else {
    // Bootstrapping or initializing without an input file.
    result.Mbz = config->get_number("grid.Mbz");
    result.Lbz = config->get_number("grid.Lbz");

    if (result.Mbz == 1) {
      result.Lbz = 0;
      result.Mbz = 1;
    }
  }

  return result;
}

/*! Allocate a complete or minimal bedrock thermal unit depending on the number of bedrock levels.
 *
 */
std::shared_ptr<BedThermalUnit> BedThermalUnit::FromOptions(std::shared_ptr<const Grid> grid,
                                                            std::shared_ptr<const Context> ctx) {

  auto bedrock_grid = BTUGrid::FromOptions(ctx);

  if (bedrock_grid.Mbz > 1) {
    return std::make_shared<energy::BTU_Full>(grid, bedrock_grid);
  }

  return std::make_shared<energy::BTU_Minimal>(grid);
}


BedThermalUnit::BedThermalUnit(std::shared_ptr<const Grid> g)
  : Component(g),
    m_bottom_surface_flux(m_grid, "bheatflx"),
    m_top_surface_flux(m_grid, "heat_flux_from_bedrock") {

  {
    m_top_surface_flux.metadata(0)
        .long_name("upward geothermal flux at the top bedrock surface")
        .units("W m-2")
        .output_units("mW m-2")
        .standard_name("upward_geothermal_heat_flux_at_ground_level_in_land_ice");
    m_top_surface_flux.metadata()["comment"] = "positive values correspond to an upward flux";
  }
  {
    // PROPOSED standard_name = lithosphere_upward_heat_flux
    m_bottom_surface_flux.metadata(0)
        .long_name("upward geothermal flux at the bottom bedrock surface")
        .units("W m-2"); // note: don't convert to
                         // "mW m-2" when saving

    m_bottom_surface_flux.metadata()["comment"] = "positive values correspond to an upward flux";
    m_bottom_surface_flux.set_time_independent(true);
  }
}

void BedThermalUnit::init(const InputOptions &opts) {
  this->init_impl(opts);
}

//! \brief Initialize the bedrock thermal unit.
void BedThermalUnit::init_impl(const InputOptions &opts) {
  auto input_file = m_config->get_string("energy.bedrock_thermal.file");

  if (not input_file.empty()) {
    m_log->message(2, "  - Reading geothermal flux from '%s' ...\n",
                   input_file.c_str());

    m_bottom_surface_flux.regrid(input_file, io::CRITICAL);
  } else {
    m_log->message(2,
                   "  - Parameter %s is not set. Reading geothermal flux from '%s'...\n",
                   "energy.bedrock_thermal.file",
                   opts.filename.c_str());

    switch (opts.type) {
    case INIT_RESTART:
      m_bottom_surface_flux.read(opts.filename, opts.record);
      break;
    case INIT_BOOTSTRAP:
      m_bottom_surface_flux.regrid(opts.filename, io::OPTIONAL,
                                   m_config->get_number("bootstrapping.defaults.geothermal_flux"));
      break;
    case INIT_OTHER:
    default:
      initialize_bottom_surface_flux();
    }
  }

  regrid("bedrock thermal layer", m_bottom_surface_flux, REGRID_WITHOUT_REGRID_VARS);
}

void BedThermalUnit::initialize_bottom_surface_flux() {
  const double heat_flux = m_config->get_number("bootstrapping.defaults.geothermal_flux");

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

void BedThermalUnit::define_model_state_impl(const File &output) const {
  m_bottom_surface_flux.define(output, io::PISM_DOUBLE);
}

void BedThermalUnit::write_model_state_impl(const File &output) const {
  m_bottom_surface_flux.write(output);
}

DiagnosticList BedThermalUnit::diagnostics_impl() const {
  DiagnosticList result = {
    {"bheatflx",   Diagnostic::wrap(m_bottom_surface_flux)},
    {"heat_flux_from_bedrock", Diagnostic::Ptr(new BTU_geothermal_flux_at_ground_level(this))}};

  if (m_config->get_flag("output.ISMIP6")) {
    result["hfgeoubed"] = Diagnostic::Ptr(new BTU_geothermal_flux_at_ground_level(this));
  }
  return result;
}

void BedThermalUnit::update(const array::Scalar &bedrock_top_temperature,
                            double t, double dt) {
  this->update_impl(bedrock_top_temperature, t, dt);
}

const array::Scalar& BedThermalUnit::flux_through_top_surface() const {
  return m_top_surface_flux;
}

const array::Scalar& BedThermalUnit::flux_through_bottom_surface() const {
  return m_bottom_surface_flux;
}

BTU_geothermal_flux_at_ground_level::BTU_geothermal_flux_at_ground_level(const BedThermalUnit *m)
  : Diag<BedThermalUnit>(m) {

  auto ismip6 = m_config->get_flag("output.ISMIP6");

  // set metadata:
  m_vars = {SpatialVariableMetadata(m_sys, ismip6 ? "hfgeoubed" : "heat_flux_from_bedrock")};

  set_attrs("upward geothermal flux at the top bedrock surface",
            (ismip6 ?
             "upward_geothermal_heat_flux_in_land_ice" :
             "upward_geothermal_heat_flux_at_ground_level_in_land_ice"),
            "W m-2", "mW m-2", 0);
  m_vars[0]["comment"] = "positive values correspond to an upward flux";
}

array::Array::Ptr BTU_geothermal_flux_at_ground_level::compute_impl() const {
  auto result = std::make_shared<array::Scalar>(m_grid, "hfgeoubed");
  result->metadata() = m_vars[0];

  result->copy_from(model->flux_through_top_surface());

  return result;
}

} // end of namespace energy
} // end of namespace pism
