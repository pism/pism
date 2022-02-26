// Copyright (C) 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2022 Constantine Khroulev
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

#include "BedDef.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/util/Time.hh"
#include "pism/util/Vars.hh"
#include "pism/util/IceGrid.hh"
#include "pism/util/ConfigInterface.hh"
#include "pism/util/Context.hh"

namespace pism {
namespace bed {

BedDef::BedDef(IceGrid::ConstPtr grid)
  : Component(grid),
    m_wide_stencil(m_config->get_number("grid.max_stencil_width")),
    m_topg(m_grid, "topg"),
    m_topg_last(m_grid, "topg"),
    m_uplift(m_grid, "dbdt")
{

  m_topg.set_attrs("model_state", "bedrock surface elevation",
                   "m", "m", "bedrock_altitude", 0);

  m_topg_last.set_attrs("model_state", "bedrock surface elevation",
                        "m", "m", "bedrock_altitude", 0);

  m_uplift.set_attrs("model_state", "bedrock uplift rate",
                     "m s-1", "mm year-1", "tendency_of_bedrock_altitude", 0);
}

const array::Scalar& BedDef::bed_elevation() const {
  return m_topg;
}

const array::Scalar& BedDef::uplift() const {
  return m_uplift;
}

void BedDef::define_model_state_impl(const File &output) const {
  m_uplift.define(output);
  m_topg.define(output);
}

void BedDef::write_model_state_impl(const File &output) const {
  m_uplift.write(output);
  m_topg.write(output);
}

DiagnosticList BedDef::diagnostics_impl() const {
  DiagnosticList result;
  result = {
    {"dbdt", Diagnostic::wrap(m_uplift)},
    {"topg", Diagnostic::wrap(m_topg)}
  };

  return result;
}

void BedDef::init(const InputOptions &opts, const array::Scalar &ice_thickness,
                  const array::Scalar &sea_level_elevation) {
  this->init_impl(opts, ice_thickness, sea_level_elevation);
}

//! Initialize using provided bed elevation and uplift.
void BedDef::bootstrap(const array::Scalar &bed_elevation,
                       const array::Scalar &bed_uplift,
                       const array::Scalar &ice_thickness,
                       const array::Scalar &sea_level_elevation) {
  this->bootstrap_impl(bed_elevation, bed_uplift, ice_thickness, sea_level_elevation);
}

void BedDef::bootstrap_impl(const array::Scalar &bed_elevation,
                            const array::Scalar &bed_uplift,
                            const array::Scalar &ice_thickness,
                            const array::Scalar &sea_level_elevation) {
  m_topg.copy_from(bed_elevation);
  m_uplift.copy_from(bed_uplift);

  // suppress a compiler warning:
  (void) ice_thickness;
  (void) sea_level_elevation;
}

void BedDef::update(const array::Scalar &ice_thickness,
                    const array::Scalar &sea_level_elevation,
                    double t, double dt) {
  this->update_impl(ice_thickness, sea_level_elevation, t, dt);
}

//! Initialize from the context (input file and the "variables" database).
void BedDef::init_impl(const InputOptions &opts, const array::Scalar &ice_thickness,
                       const array::Scalar &sea_level_elevation) {
  (void) ice_thickness;
  (void) sea_level_elevation;

  switch (opts.type) {
  case INIT_RESTART:
    // read bed elevation and uplift rate from file
    m_log->message(2,
                   "    reading bed topography and uplift from %s ... \n",
                   opts.filename.c_str());
    // re-starting
    m_topg.read(opts.filename, opts.record); // fails if not found!
    m_uplift.read(opts.filename, opts.record); // fails if not found!
    break;
  case INIT_BOOTSTRAP:
    // bootstrapping
    m_topg.regrid(opts.filename, OPTIONAL,
                  m_config->get_number("bootstrapping.defaults.bed"));
    m_uplift.regrid(opts.filename, OPTIONAL,
                    m_config->get_number("bootstrapping.defaults.uplift"));
    break;
  case INIT_OTHER:
  default:
    {
      // do nothing
    }
  }

  // process -regrid_file and -regrid_vars
  regrid("bed deformation", m_topg);
  // uplift is not a part of the model state, but the user may want to take it from a -regrid_file
  // during bootstrapping
  regrid("bed deformation", m_uplift);

  std::string uplift_file = m_config->get_string("bed_deformation.bed_uplift_file");
  if (not uplift_file.empty()) {
    m_log->message(2,
                   "    reading bed uplift from %s ... \n",
                   uplift_file.c_str());
    m_uplift.regrid(uplift_file, CRITICAL);
  }

  std::string correction_file = m_config->get_string("bed_deformation.bed_topography_delta_file");
  if (not correction_file.empty()) {
    apply_topg_offset(correction_file);
  }

  // this should be the last thing we do here
  m_topg_last.copy_from(m_topg);
}

/*!
 * Apply a correction to the bed topography by reading topg_delta from filename.
 */
void BedDef::apply_topg_offset(const std::string &filename) {
  m_log->message(2, "  Adding a bed topography correction read in from %s...\n",
                 filename.c_str());

  array::Scalar topg_delta(m_grid, "topg_delta");
  topg_delta.set_attrs("internal", "bed topography correction",
                       "meters", "meters", "", 0);

  topg_delta.regrid(filename, CRITICAL);

  m_topg.add(1.0, topg_delta);
}

//! Compute bed uplift (dt is in seconds).
void BedDef::compute_uplift(const array::Scalar &bed, const array::Scalar &bed_last,
                            double dt, array::Scalar &result) {
  bed.add(-1, bed_last, result);
  //! uplift = (topg - topg_last) / dt
  result.scale(1.0 / dt);
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

  IceModelVec::AccessList list{&bed_elevation, &ice_thickness, &sea_level_elevation, &result};

  for (Points p(*result.grid()); p; p.next()) {
    const int i = p.i(), j = p.j();

    result(i, j) = compute_load(bed_elevation(i, j),
                                ice_thickness(i, j),
                                sea_level_elevation(i, j),
                                ice_density, ocean_density);
  }
}

} // end of namespace bed
} // end of namespace pism
