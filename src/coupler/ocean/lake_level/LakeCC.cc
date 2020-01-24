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
#include "pism/coupler/util/options.hh"

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

  m_icefree_thickness = m_config->get_number("geometry.ice_free_thickness_standard");

  const double ice_density        = m_config->get_number("constants.ice.density"),
               freshwater_density = m_config->get_number("constants.fresh_water.density");
  m_drho = ice_density / freshwater_density;

  m_filter_map = true;
  m_n_filter   = 3;

  m_keep_existing_lakes = false;
  m_check_sl_diagonal = false;
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
                "meter", "meter", "", 0);
  tmp.metadata().set_number("_FillValue", m_fill_value);

  {
    InputOptions opts = process_input_options(m_grid->com, m_config);

    if (opts.type == INIT_RESTART) {

      m_log->message(2, "* Reading lake level forcing from '%s' for re-starting...\n",
                    opts.filename.c_str());

      File file(m_grid->com, opts.filename, PISM_GUESS, PISM_READONLY);
      const unsigned int time_length = file.nrecords(),
                        last_record = time_length > 0 ? time_length - 1 : 0;

      tmp.read(file, last_record);

      file.close();
    } else {
      tmp.set(m_fill_value);
    }
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

  m_keep_existing_lakes = options::Bool(m_option + "_keep_existing_lakes",
                                        "If set lakes are kept even though they are enclosed or covered by ice. This might result in huge sub-glacial lakes.");
  m_check_sl_diagonal = options::Bool(m_option + "_check_sl_diagonal",
                                      "If set invalid SL cells are checked if a diagonal cell is ocean and if so, the mask is updated.");
  m_topg_overlay.create(m_grid, "topg_overlay", WITHOUT_GHOSTS);
  m_topg_overlay.set_attrs("internal",
                           "topography overlay",
                           "meter", "meter", "", 0);
  {
    ForcingOptions opt(*m_grid->ctx(), "lake_level.lakecc.topg_overlay");

    if (not opt.filename.empty()) {
      m_topg_overlay.regrid(opt.filename, OPTIONAL, 0.0);
    } else {
      m_topg_overlay.set(0.0);
    }
  }

  m_bed.create(m_grid, "topg_lakecc", WITHOUT_GHOSTS);
  m_bed.set_attrs("diagnostic",
                "bed topography as seen by LakeCC model",
                "meter", "meter", "", 0);
}

void LakeCC::update_impl(const Geometry &geometry, double my_t, double my_dt) {

  const IceModelVec2S &thk = geometry.ice_thickness,
                      &ll  = geometry.lake_level_elevation,
                      &sl  = *m_grid->variables().get_2d_scalar("sea_level");
  m_topg_overlay.add(1.0, geometry.bed_elevation, m_bed);

  do_lake_update(m_bed, thk, sl, ll);

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

  if (m_keep_existing_lakes) {
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
  }

  valid_mask.update_ghosts();
}

void LakeCC::do_lake_update(const IceModelVec2S &bed, const IceModelVec2S &thk, const IceModelVec2S &sea_level, const IceModelVec2S &old_lake_level) {
  IceModelVec2Int pism_mask;
  pism_mask.create(m_grid, "pism_mask", WITHOUT_GHOSTS);
  m_gc.compute_mask(sea_level, bed, thk, pism_mask);

  if (m_check_sl_diagonal) {
    update_mask_sl_diagonal(sea_level, bed, thk, pism_mask);
  }

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
    { "lakecc_topg",           Diagnostic::wrap(m_bed) },
  };

  return result;
}

bool LakeCC::expandMargins_impl() const {
  return true;
}

void LakeCC::update_mask_sl_diagonal(const IceModelVec2S &sl, const IceModelVec2S &bed, const IceModelVec2S &thk, IceModelVec2Int &mask) {

  IceModelVec2Int mask_wide(m_grid, "mask_wide", WITH_GHOSTS);
  mask_wide.copy_from(mask);

  IceModelVec2S sl_wide(m_grid, "sl_wide", WITH_GHOSTS);
  sl_wide.copy_from(sl);

  std::vector<int>i_diagonals = {1, 1,-1,-1},
                  j_diagonals = {1,-1,-1, 1};

  IceModelVec::AccessList list({ &mask_wide, &mask, &sl_wide, &bed, &thk });

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (sl_wide(i, j) == m_fill_value) {
      for (int k = 0; k < 4; k++) {
        int i_diag = i + i_diagonals[k],
            j_diag = j + j_diagonals[k];
        if (mask::ocean(mask_wide(i_diag, j_diag))) {
          mask(i, j) = m_gc.mask(sl_wide(i_diag, j_diag), bed(i, j), thk(i, j));
          break;
        }
      }
    }
  }

}

} // end of namespace lake_level
} // end of namespace ocean
} // end of namespace pism

