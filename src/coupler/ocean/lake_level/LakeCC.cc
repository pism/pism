/* Copyright (C) 2013, 2014, 2015, 2016 PISM Authors
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
#include <cassert>
#include <vector>

#include "pism/util/Vars.hh"
#include "pism/util/pism_options.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/util/MaxTimestep.hh"
#include "pism/geometry/Geometry.hh"

#include "LakeCC.hh"
#include "LakeLevel_ConnectedComponents.hh"

namespace pism {
namespace ocean {
namespace lake_level {

LakeCC::LakeCC(IceGrid::ConstPtr g)
  : LakeLevel(g),
    m_gc(*m_config) {

  m_log->message(2, "  - Setting up LakeCC Module...\n");

  m_option = "-lakecc";

  m_lake_level_min = 0.;
  m_lake_level_max = 1000.;
  m_lake_level_dh  = 10.;

  m_icefree_thickness = m_config->get_double("geometry.ice_free_thickness_standard");

  const double ice_density        = m_config->get_double("constants.ice.density"),
               freshwater_density = m_config->get_double("constants.fresh_water.density");
  m_drho = ice_density / freshwater_density;

  m_filter_map = true;
  m_n_filter   = 3;
}

LakeCC::~LakeCC() {
  // empty
}

void LakeCC::init_impl(const Geometry &geometry) {

  m_log->message(2, "* Initializing the lakeCC ocean model modifier...\n");

  IceModelVec2S tmp;
  tmp.create(m_grid, "effective_lake_level_elevation", WITHOUT_GHOSTS);
  tmp.set_attrs("diagnostic",
                "lake level elevation, relative to the geoid",
                "meter", "");
  tmp.metadata().set_double("_FillValue", m_fill_value);

  InputOptions opts = process_input_options(m_grid->com, m_config);

  if (opts.type == INIT_RESTART) {

    m_log->message(2, "* Reading lake level forcing from '%s' for re-starting...\n",
                   opts.filename.c_str());

    PIO file(m_grid->com, "guess_mode", opts.filename, PISM_READONLY);
    const unsigned int time_length = file.inq_nrecords(),
                       last_record = time_length > 0 ? time_length - 1 : 0;

    tmp.read(file, last_record);

    file.close();
  } else {
    tmp.set(m_fill_value);
  }

  // Support regridding. This is needed to ensure that initialization using "-i" is
  // equivalent to "-i ... -bootstrap -regrid_file ..."
  {
    regrid("lakecc lake model", tmp,
           REGRID_WITHOUT_REGRID_VARS);
  }

  m_lake_level.copy_from(tmp);

  std::vector<double> default_levels;
  default_levels.push_back(m_lake_level_min);
  default_levels.push_back(m_lake_level_max);
  default_levels.push_back(m_lake_level_dh);
  options::RealList levels(m_option + "_level", "Comma-seperated list of three numbers: z_min, z_max, dz",
                           default_levels);
  m_lake_level_min = levels[0];
  m_lake_level_max = levels[1];
  m_lake_level_dh  = levels[2];

  m_log->message(2, "  LakeCC: lake levels between %g and %gm, with %gm spacing\n",
                 m_lake_level_min, m_lake_level_max, m_lake_level_dh);

  double icefree_thickness = m_icefree_thickness;
  icefree_thickness        = options::Real(m_option + "_icefree_thickness",
                                           "Minimum ice thickness that is regarded as a closed ice cover.", icefree_thickness);
  m_icefree_thickness      = icefree_thickness;
  m_log->message(2, "  LakeCC: icefree thickness: %gm \n", m_icefree_thickness);

  int n_filter = m_n_filter;
  n_filter = options::Integer(m_option + "_n_filter", "Only keep lakes that contain at least one cell that has at least N lake neighbors. (0: no filtering, 1(remove one-cell lakes only)-4(keep lakes that have at least one surrounded cell))", n_filter);
  if (n_filter < 0) {
    n_filter = 0;
  } else if (n_filter > 4) {
    n_filter = 4;
  }
  m_n_filter = n_filter;
  if (m_n_filter == 0) {
    m_filter_map = false;
  }

  if (m_filter_map) {
    m_log->message(2, "  LakeCC: Filtering map (N=%i) \n", m_n_filter);
  }

  m_valid_mask.create(m_grid, "lake_valid_mask", WITHOUT_GHOSTS);
}

void LakeCC::update_impl(const Geometry &geometry, double my_t, double my_dt) {

  const IceModelVec2S &bed = geometry.bed_elevation,
                      &thk = geometry.ice_thickness,
                      &ll  = geometry.lake_level_elevation,
                      &sl  = *m_grid->variables().get_2d_scalar("sea_level");

  do_lake_update(bed, thk, sl, ll);

  if (m_filter_map) {
    do_filter_map();
  }
}

MaxTimestep LakeCC::max_timestep_impl(double t) const {
  return MaxTimestep("lake level forcing");
}

void LakeCC::prepare_mask_validity(const IceModelVec2S &thk, const IceModelVec2S &bed, const IceModelVec2S &old_lake_level, IceModelVec2Int& valid_mask) {
  IsolationCC IsoCC(m_grid, thk, m_icefree_thickness);
  IsoCC.find_isolated_spots(valid_mask);

  IceModelVec::AccessList list{ &valid_mask, &bed, &thk, &old_lake_level };

  if (old_lake_level.state_counter() > 0) {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      //Set valid where a lake already exists
      if ( mask::ocean(m_gc.mask(m_fill_value, bed(i, j), thk(i, j), old_lake_level(i, j))) ) {
        valid_mask(i, j) = 1;
      }
    }
  }

  valid_mask.update_ghosts();
}

void LakeCC::do_lake_update(const IceModelVec2S &bed, const IceModelVec2S &thk, const IceModelVec2S &sea_level, const IceModelVec2S &old_lake_level) {
  IceModelVec2Int pism_mask;
  pism_mask.create(m_grid, "pism_mask", WITHOUT_GHOSTS);
  m_gc.compute_mask(sea_level, bed, thk, pism_mask);
  prepare_mask_validity(thk, bed, old_lake_level, m_valid_mask);

  m_log->message(2, "->LakeCC: Update of Lake Levels! \n");

  ParallelSection ParSec(m_grid->com);
  try {
    // Initialze LakeCC Model
    LakeLevelCC LM(m_grid, m_drho, bed, thk, pism_mask, m_fill_value, m_valid_mask);
    LM.computeLakeLevel(m_lake_level_min, m_lake_level_max, m_lake_level_dh, m_lake_level);
  } catch (...) {
    ParSec.failed();
  }
  ParSec.check();

  m_log->message(2, "          Done!\n");
}

void LakeCC::do_filter_map() {
  ParallelSection ParSec(m_grid->com);
  try {
    FilterLakesCC FL(m_grid, m_fill_value);
    FL.filter_map(m_n_filter, m_lake_level);
  } catch (...) {
    ParSec.failed();
  }
  ParSec.check();
}

// Write diagnostic variables to extra files if requested
DiagnosticList LakeCC::diagnostics_impl() const {

  DiagnosticList result = {
    { "lake_valid_mask",       Diagnostic::wrap(m_valid_mask) },
  };

  return result;
}

const bool LakeCC::expandMargins_impl() const {
  return true;
}

} // end of namespace lake_level
} // end of namespace ocean
} // end of namespace pism

