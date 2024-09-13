// Copyright (C) 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2022, 2023, 2024 Constantine Khroulev
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

#include "pism/earth/BedDef.hh"
#include "pism/util/ConfigInterface.hh"
#include "pism/util/Context.hh"
#include "pism/util/Grid.hh"
#include "pism/util/MaxTimestep.hh"
#include "pism/util/Time.hh"
#include <string>

namespace pism {
namespace bed {

BedDef::BedDef(std::shared_ptr<const Grid> grid, const std::string &model_name)
  : Component(grid),
    m_topg(m_grid, "topg"),
    m_topg_last(m_grid, "topg"),
    m_load(grid, "bed_def_load"),
    m_load_accumulator(grid, "bed_def_load_accumulator"),
    m_uplift(m_grid, "dbdt"),
    m_model_name(model_name)
{

  m_time_name = m_config->get_string("time.dimension_name") + "_bed_deformation";
  m_t_last = time().current();
  m_update_interval = m_config->get_number("bed_deformation.update_interval", "seconds");
  m_t_eps = m_config->get_number("time_stepping.resolution", "seconds");

  m_topg.metadata(0)
      .long_name("bedrock surface elevation")
      .units("m")
      .standard_name("bedrock_altitude");

  m_topg_last.metadata(0)
      .long_name("bedrock surface elevation after the previous update (used to compute uplift)")
      .units("m")
      .standard_name("bedrock_altitude");

  m_load.metadata(0)
    .long_name("load on the bed expressed as ice-equivalent thickness")
    .units("m");

  m_load_accumulator.metadata(0)
    .long_name("accumulated load on the bed expressed as a time integral of ice-equivalent thickness")
    .units("m s");

  m_uplift.metadata(0)
      .long_name("bedrock uplift rate")
      .units("m s^-1")
      .standard_name("tendency_of_bedrock_altitude");
}

const array::Scalar &BedDef::bed_elevation() const {
  return m_topg;
}

const array::Scalar &BedDef::uplift() const {
  return m_uplift;
}

void BedDef::define_model_state_impl(const File &output) const {
  m_uplift.define(output, io::PISM_DOUBLE);
  m_topg.define(output, io::PISM_DOUBLE);

  if (not output.variable_exists(m_time_name)) {
    output.define_variable(m_time_name, io::PISM_DOUBLE, {});

    output.write_attribute(m_time_name, "long_name",
                        "time of the last update of the Lingle-Clark bed deformation model");
    output.write_attribute(m_time_name, "calendar", time().calendar());
    output.write_attribute(m_time_name, "units", time().units_string());
  }
}

void BedDef::write_model_state_impl(const File &output) const {
  m_uplift.write(output);
  m_topg.write(output);

  output.write_variable(m_time_name, {0}, {1}, &m_t_last);
}

DiagnosticList BedDef::diagnostics_impl() const {
  DiagnosticList result;
  result = { { "dbdt", Diagnostic::wrap(m_uplift) }, { "topg", Diagnostic::wrap(m_topg) } };

  return result;
}

void BedDef::init(const InputOptions &opts, const array::Scalar &ice_thickness,
                  const array::Scalar &sea_level_elevation) {

  m_log->message(2, "* Initializing the %s bed deformation model...\n",
                 m_model_name.c_str());

  if (opts.type == INIT_RESTART or opts.type == INIT_BOOTSTRAP) {
    File input_file(m_grid->com, opts.filename, io::PISM_NETCDF3, io::PISM_READONLY);

    if (input_file.variable_exists(m_time_name)) {
      input_file.read_variable(m_time_name, {0}, {1}, &m_t_last);
    } else {
      m_t_last = time().current();
    }
  } else {
    m_t_last = time().current();
  }

  {
    switch (opts.type) {
    case INIT_RESTART:
      // read bed elevation and uplift rate from file
      m_log->message(2, "    reading bed topography and uplift from %s ... \n",
                     opts.filename.c_str());
      // re-starting
      m_topg.read(opts.filename, opts.record);   // fails if not found!
      m_uplift.read(opts.filename, opts.record); // fails if not found!
      break;
    case INIT_BOOTSTRAP:
      // bootstrapping
      m_topg.regrid(opts.filename, io::Default(m_config->get_number("bootstrapping.defaults.bed")));
      m_uplift.regrid(opts.filename,
                      io::Default(m_config->get_number("bootstrapping.defaults.uplift")));
      break;
    case INIT_OTHER:
    default: {
      // do nothing
    }
    }

    // process -regrid_file and -regrid_vars
    regrid("bed deformation", m_topg);
    // uplift is not a part of the model state, but the user may want to take it from a -regrid_file
    // during bootstrapping
    regrid("bed deformation", m_uplift);

    auto uplift_file = m_config->get_string("bed_deformation.bed_uplift_file");
    if (not uplift_file.empty()) {
      m_log->message(2, "    reading bed uplift from %s ... \n", uplift_file.c_str());
      m_uplift.regrid(uplift_file, io::Default::Nil());
    }

    auto correction_file = m_config->get_string("bed_deformation.bed_topography_delta_file");
    if (not correction_file.empty()) {
      m_log->message(2, "  Adding a bed topography correction read in from '%s'...\n",
                     correction_file.c_str());
      apply_topg_offset(correction_file, m_topg);
    }
  }

  this->init_impl(opts, ice_thickness, sea_level_elevation);

  // this should be the last thing we do
  m_topg_last.copy_from(m_topg);
}

//! Initialize using provided bed elevation and uplift.
void BedDef::bootstrap(const array::Scalar &bed_elevation, const array::Scalar &bed_uplift,
                       const array::Scalar &ice_thickness,
                       const array::Scalar &sea_level_elevation) {
  m_t_last = time().current();

  m_topg.copy_from(bed_elevation);
  m_uplift.copy_from(bed_uplift);
  m_topg_last.copy_from(bed_elevation);

  this->bootstrap_impl(bed_elevation, bed_uplift, ice_thickness, sea_level_elevation);
}

void BedDef::bootstrap_impl(const array::Scalar & /*bed_elevation*/,
                           const array::Scalar & /*bed_uplift*/,
                           const array::Scalar & /*ice_thickness*/,
                           const array::Scalar & /*sea_level_elevation*/) {
  // empty
}

void BedDef::update(const array::Scalar &ice_thickness, const array::Scalar &sea_level_elevation,
                    double t, double dt) {

  double t_final = t + dt;

  if (t_final < m_t_last) {
    throw RuntimeError(PISM_ERROR_LOCATION, "cannot go back in time");
  }

  compute_load(m_topg_last, ice_thickness, sea_level_elevation, m_load);
  m_load_accumulator.add(dt, m_load);

  double t_next = m_update_interval > 0.0 ? m_t_last + m_update_interval : t_final;
  if (std::abs(t_next - t_final) < m_t_eps) { // reached the next update time

    double dt_beddef = t_final - m_t_last;

    // compute time-averaged load and reset the accumulator
    {
      m_load.copy_from(m_load_accumulator);
      m_load.scale(1.0 / dt_beddef);

      m_load_accumulator.set(0.0);
    }

    this->update_impl(m_load, m_t_last, dt_beddef);
    // note: we don't know if a derived class modified m_topg in update_impl(), so we *do
    // not* call m_topg.inc_state_counter() here -- it should be done in update_impl(), if
    // necessary

    m_t_last = t_final;

    // Update m_uplift and m_topg_last
    {
      m_topg.add(-1, m_topg_last, m_uplift);
      //! uplift = (m_topg - m_topg_last) / dt
      m_uplift.scale(1.0 / dt_beddef);
    }
    m_topg_last.copy_from(m_topg);
  }
}

MaxTimestep BedDef::max_timestep_impl(double t) const {

  if (t < m_t_last) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "time %f is less than the previous time %f",
                                  t, m_t_last);
  }

  if (m_update_interval == 0.0) {
    return {};
  }

  // Find the smallest time of the form m_t_last + k * m_update_interval that is greater
  // than t
  double k = std::ceil((t - m_t_last) / m_update_interval);

  double t_next = m_t_last + k * m_update_interval;

  double dt_max = m_update_interval;
  if (t < t_next) {
    dt_max = t_next - t;
  }

  if (dt_max < m_t_eps) {
    dt_max = m_update_interval;
  }

  return MaxTimestep(dt_max, "bed_def");
}

/*!
 * Apply a correction to the bed topography by reading "topg_delta" from `filename`.
 */
void BedDef::apply_topg_offset(const std::string &filename, array::Scalar &bed_topography) {

  auto grid = bed_topography.grid();

  array::Scalar topg_delta(grid, "topg_delta");
  topg_delta.metadata(0).long_name("bed topography correction").units("meters");

  topg_delta.regrid(filename, io::Default::Nil());

  bed_topography.add(1.0, topg_delta);
}

double compute_load(double bed, double ice_thickness, double sea_level,
                    double ice_density, double ocean_density) {

  double
    ice_load    = ice_thickness,
    ocean_depth = std::max(sea_level - bed, 0.0),
    ocean_load  = (ocean_density / ice_density) * ocean_depth;

  // this excludes the load of ice shelves
  return ice_load > ocean_load ? ice_load : 0.0;
}

/*! Compute the load on the bedrock in units of ice-equivalent thickness.
 *
 */
void compute_load(const array::Scalar &bed_elevation,
                  const array::Scalar &ice_thickness,
                  const array::Scalar &sea_level_elevation,
                  array::Scalar &result) {

  Config::ConstPtr config = result.grid()->ctx()->config();

  const double
    ice_density   = config->get_number("constants.ice.density"),
    ocean_density = config->get_number("constants.sea_water.density");

  array::AccessScope list{&bed_elevation, &ice_thickness, &sea_level_elevation, &result};

  for (auto p = result.grid()->points(); p; p.next()) {
    const int i = p.i(), j = p.j();

    result(i, j) = compute_load(bed_elevation(i, j),
                                ice_thickness(i, j),
                                sea_level_elevation(i, j),
                                ice_density, ocean_density);
  }
}

} // end of namespace bed
} // end of namespace pism
