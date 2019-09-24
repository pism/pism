/* Copyright (C) 2018 PISM Authors
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

#include "pism/util/pism_options.hh"
#include "pism/util/Vars.hh"
#include "pism/geometry/Geometry.hh"

#include "Gradual.hh"
#include "LakeLevel_ConnectedComponents.hh"

namespace pism {
namespace ocean {
namespace lake_level {

Gradual::Gradual(IceGrid::ConstPtr grid,
                 std::shared_ptr<LakeLevel> in)
  : LakeLevel(grid, in) {

  m_use_const_fill_rate = false;

  //Set default filling rate for lakes
  m_max_lake_fill_rate = units::convert(m_sys, 1.0, "m/years", "m/seconds");
  m_target_level.create(m_grid, "target_level", WITHOUT_GHOSTS);
  m_target_level.set_attrs("model_state", "target lake level",
                           "m", "target_level");
  m_target_level.metadata().set_double("_FillValue", m_fill_value);

  m_min_level.create(m_grid, "min_lake_level", WITHOUT_GHOSTS);
  m_min_level.set_attrs("model_state", "min lake level",
                        "m", "min_level");
  m_min_level.metadata().set_double("_FillValue", m_fill_value);

  m_max_level.create(m_grid, "max_lake_level", WITHOUT_GHOSTS);
  m_max_level.set_attrs("model_state", "max lake level",
                        "m", "max_level");
  m_max_level.metadata().set_double("_FillValue", m_fill_value);

  m_min_basin.create(m_grid, "min_lake_bed", WITHOUT_GHOSTS);
  m_min_basin.set_attrs("model_state", "min lake bed",
                        "m", "min_bed");
  m_min_basin.metadata().set_double("_FillValue", m_fill_value);

  m_max_sl_basin.create(m_grid, "max_sl_basin", WITHOUT_GHOSTS);
  m_max_sl_basin.set_attrs("model_state", "max sl basin",
                           "m", "max_sl_basin");
  m_max_sl_basin.metadata().set_double("_FillValue", m_fill_value);

  m_expansion_mask.create(m_grid, "expansion_mask", WITHOUT_GHOSTS);
  m_expansion_mask.metadata().set_double("_FillValue", m_fill_value);

  m_lake_area.create(m_grid, "lake_surface_area", WITHOUT_GHOSTS);
  m_lake_area.set_attrs("model_state", "lake_surface_area",
                        "m2", "lake_surface");
  m_lake_area.metadata().set_double("_FillValue", m_fill_value);

  m_lake_mass_input_discharge.create(m_grid, "lake_mass_input_due_to_discharge", WITHOUT_GHOSTS);
  m_lake_mass_input_discharge.set_attrs("model_state", "lake_mass_input_discharge",
                                        "kg s-1", "lake_mass_input_discharge");
  m_lake_mass_input_discharge.metadata().set_double("_FillValue", m_fill_value);
  m_lake_mass_input_discharge.metadata().set_string("glaciological_units", "kg year-1");

  m_lake_mass_input_basal.create(m_grid, "lake_mass_input_due_to_basal melt", WITHOUT_GHOSTS);
  m_lake_mass_input_basal.set_attrs("model_state", "lake_mass_input_basal",
                                    "kg s-1", "lake_mass_input_shelf");
  m_lake_mass_input_basal.metadata().set_double("_FillValue", m_fill_value);
  m_lake_mass_input_basal.metadata().set_string("glaciological_units", "kg year-1");

  m_lake_mass_input_total.create(m_grid, "lake_mass_input_total", WITHOUT_GHOSTS);
  m_lake_mass_input_total.set_attrs("model_state", "lake_mass_input_total",
                                    "kg s-1", "lake_mass_input_total");
  m_lake_mass_input_total.metadata().set_double("_FillValue", m_fill_value);
  m_lake_mass_input_total.metadata().set_string("glaciological_units", "kg year-1");

  m_lake_fill_rate.create(m_grid, "lake_fill_rate", WITHOUT_GHOSTS);
  m_lake_fill_rate.set_attrs("model_state", "lake_fill_rate",
                             "m s-1", "lake_level_rise");
  m_lake_fill_rate.metadata().set_double("_FillValue", m_fill_value);
  m_lake_fill_rate.metadata().set_string("glaciological_units", "m year-1");

  m_alpha_lake = m_config->get_double("constants.ice.density") / m_config->get_double("constants.fresh_water.density");

  m_init_lakes_filled = false;
}

Gradual::~Gradual() {
  //empty
}

void Gradual::init_impl(const Geometry &geometry) {
  m_input_model->init(geometry);

  bool regridded = false;
  bool restarted = false;
  {
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

      restarted = true;
    } else {
      tmp.set(m_fill_value);
    }

    // Support regridding. This is needed to ensure that initialization using "-i" is
    // equivalent to "-i ... -bootstrap -regrid_file ..."
    {
      regrid("lake gradual filling modifier", tmp,
             REGRID_WITHOUT_REGRID_VARS);
    }

    if (tmp.state_counter() == 2) {
      regridded = true;
    }

    m_lake_level.copy_from(tmp);
  }

  bool bootstrap_filled = options::Bool("-lake_gradual_bootstrap_filled",
                                        "When bootstrapping at model start, initialize lakes filled.");

  if (bootstrap_filled and not restarted and not regridded) {
    m_init_lakes_filled = true;
  }

  double max_lake_fill_rate_m_y = units::convert(m_sys, m_max_lake_fill_rate,
                                                 "m/seconds", "m/year");
  max_lake_fill_rate_m_y = options::Real("-lake_gradual_max_fill_rate",
                                         "Maximum rate at which lakes do fill (m/year)",
                                         max_lake_fill_rate_m_y);
  m_max_lake_fill_rate = units::convert(m_sys, max_lake_fill_rate_m_y, "m/year", "m/seconds");

  bool const_fill_rate = m_use_const_fill_rate;
  const_fill_rate = options::Bool("-lake_gradual_const_fill_rate",
                                  "Use constant lake fill rate instead of determining it from ice sheet discharge");
  m_use_const_fill_rate = const_fill_rate;
}

void Gradual::update_impl(const Geometry &geometry, double t, double dt) {
  m_input_model->update(geometry, t, dt);

  m_target_level.copy_from(m_input_model->elevation());

  if (m_init_lakes_filled) {
    m_lake_level.copy_from(m_target_level);
    m_init_lakes_filled = false;
  }

  //Get bed, thk and the 'old' lake and sea level fields from geometry
  //Since sea level was updated just before the lake level is, sea_level_elevation
  //in geometry might not be up-to-date, so we get it from *m_grid->variables()
  const IceModelVec2S &bed    = geometry.bed_elevation,
                      &thk    = geometry.ice_thickness,
                      &old_ll = geometry.lake_level_elevation,
                      &sl     = *m_grid->variables().get_2d_scalar("sea_level"),
                      &old_sl = geometry.sea_level_elevation;

  //Calculate basins that were filled by the ocean
  IceModelVec2S old_sl_basins(m_grid, "sl_basins", WITHOUT_GHOSTS);
  {
    IceModelVec::AccessList list {&old_sl_basins, &old_sl, &bed, &thk};

    GeometryCalculator gc(*m_config);

    ParallelSection ParSec(m_grid->com);
    try {
      for (Points p(*m_grid); p; p.next()) {
        const int i = p.i(), j = p.j();
        if (mask::ocean(gc.mask(old_sl(i, j), bed(i, j), thk(i, j)))) {
          old_sl_basins(i, j) = old_sl(i, j);
        } else {
          old_sl_basins(i, j) = m_fill_value;
        }
      }
    } catch (...) {
      ParSec.failed();
    }
    ParSec.check();

    old_sl_basins.update_ghosts();
  }

  if (not m_use_const_fill_rate) {
    compute_fill_rate(dt, m_target_level, m_lake_area, m_lake_mass_input_discharge, m_lake_mass_input_basal, m_lake_mass_input_total, m_lake_fill_rate);
  }

  updateLakeLevelMinMax(m_lake_level, m_target_level, m_min_level, m_max_level);

  const bool LakeLevelChanged = prepareLakeLevel(m_target_level, bed, m_min_level, old_ll, old_sl_basins, m_expansion_mask, m_min_basin, m_max_sl_basin, m_lake_level);

  if (LakeLevelChanged) {
    //if a new lake basin was added we need to update the min and max lake level
    updateLakeLevelMinMax(m_lake_level, m_target_level, m_min_level, m_max_level);
  }

  gradually_fill(dt, m_max_lake_fill_rate, m_target_level, bed, thk, sl, m_min_level, m_max_level, m_min_basin, m_lake_fill_rate, m_lake_level);
}

void Gradual::compute_fill_rate(const double dt,
                                const IceModelVec2S &lake_level,
                                IceModelVec2S &lake_area,
                                IceModelVec2S &lake_mass_input_discharge,
                                IceModelVec2S &lake_mass_input_basal,
                                IceModelVec2S &lake_mass_input_total,
                                IceModelVec2S &lake_fill_rate) {

  {
    //Initialize Lake accumulator
    LakeAccumulatorCCSerial Lacc(m_grid, m_fill_value);
    Lacc.init(lake_level);

    const IceModelVec2S &cell_area = *m_grid->variables().get_2d_scalar("cell_area"),
                        &discharge = *m_grid->variables().get_2d_scalar("discharge"),
                        &bmb       = *m_grid->variables().get_2d_scalar("effective_BMB");

    IceModelVec2S mass_discharge(m_grid, "mass_discharge", WITHOUT_GHOSTS),
                  mass_basal(m_grid, "mass_basal", WITHOUT_GHOSTS);

    double rho_ice = m_config->get_double("constants.ice.density");

    IceModelVec::AccessList list{ &cell_area, &discharge, &bmb, &mass_discharge, &mass_basal };

    //Update lake extend depending on exp_mask
    ParallelSection ParSec(m_grid->com);
    try {
      for (Points p(*m_grid); p; p.next()) {
        const int i = p.i(), j = p.j();
        const double cell_area_ij = cell_area(i, j);

        //Convert from basal mass balance to basal mass loss [m].
        double bml_ij = -bmb(i, j);
        bml_ij = (bml_ij >= 0.0) ? bml_ij : 0.0;

        double discharge_ij = -discharge(i, j);
        discharge_ij = (discharge_ij >= 0.0) ? discharge_ij : 0.0;

        const double C = cell_area_ij * rho_ice / dt;

        //Convert to mass flux [kg/s]
        mass_discharge(i, j) = C * discharge_ij;
        mass_basal(i, j)     = C * bml_ij;
      }
    } catch (...) {
      ParSec.failed();
    }
    ParSec.check();

    //Lake surface area
    Lacc.accumulate(cell_area, lake_area);

    Lacc.accumulate(mass_discharge, lake_mass_input_discharge);
    Lacc.accumulate(mass_basal, lake_mass_input_basal);
  }

  {
    double rho_fresh_water = m_config->get_double("constants.fresh_water.density");

    IceModelVec::AccessList list{ &lake_area, &lake_mass_input_discharge, &lake_mass_input_basal, &lake_mass_input_total, &lake_fill_rate };

    GeometryCalculator gc(*m_config);

    //Update lake extend depending on exp_mask
    ParallelSection ParSec(m_grid->com);
    try {
      for (Points p(*m_grid); p; p.next()) {
        const int i = p.i(), j = p.j();
        const double lake_area_ij = lake_area(i, j);
        if (gc.islake(lake_area_ij)){
          const double lake_mass_input_total_ij = lake_mass_input_discharge(i, j) + lake_mass_input_basal(i, j);
          lake_mass_input_total(i, j) = lake_mass_input_total_ij;
          lake_fill_rate(i, j)        = lake_mass_input_total_ij / (rho_fresh_water * lake_area_ij);
        } else {
          lake_mass_input_total(i, j) = m_fill_value;
          lake_fill_rate(i, j)        = m_fill_value;
        }
      }
    } catch (...) {
      ParSec.failed();
    }
    ParSec.check();
  }
}

void Gradual::updateLakeLevelMinMax(const IceModelVec2S &lake_level,
                                    const IceModelVec2S &target_level,
                                    IceModelVec2S &min_level,
                                    IceModelVec2S &max_level) {
  //Compute min max level
  ParallelSection ParSec(m_grid->com);
  try {
    // Initialze LakeProperties Model
    LakePropertiesCC LpCC(m_grid, m_fill_value, target_level, lake_level);
    LpCC.getLakeProperties(min_level, max_level);
  } catch (...) {
    ParSec.failed();
  }
  ParSec.check();
}

bool Gradual::prepareLakeLevel(const IceModelVec2S &target_level,
                               const IceModelVec2S &bed,
                               const IceModelVec2S &min_level,
                               const IceModelVec2S &old_ll,
                               const IceModelVec2S &old_sl,
                               IceModelVec2Int &mask,
                               IceModelVec2S &min_basin,
                               IceModelVec2S &max_sl_basin,
                               IceModelVec2S &lake_level) {

  { //Check which lake cells are newly added
    ParallelSection ParSec(m_grid->com);
    try {
      FilterExpansionCC FExCC(m_grid, m_fill_value, bed, old_sl);
      FExCC.filter_ext(lake_level, target_level, mask, min_basin, max_sl_basin);
    } catch (...) {
      ParSec.failed();
    }
    ParSec.check();
  }

  bool MinMaxChanged = false;
  {
    IceModelVec::AccessList list{ &lake_level, &min_level,
                                  &min_basin, &mask,
                                  &old_ll, &max_sl_basin };

    GeometryCalculator gc(*m_config);

    //Update lake extend depending on exp_mask
    ParallelSection ParSec(m_grid->com);
    try {
      for (Points p(*m_grid); p; p.next()) {
        const int i = p.i(), j = p.j();

        const int mask_ij = mask.as_int(i, j);
        const bool new_lake = ( (mask_ij > 0) and not (gc.islake(min_level(i, j))) );
        if ( mask_ij == 1 or new_lake ) {
          if (max_sl_basin(i, j) != m_fill_value) {
            lake_level(i, j) = max_sl_basin(i, j);
          } else {
            if (new_lake) {
              lake_level(i, j) = min_basin(i, j);
            } else {
              //New basin added to existing lake
              lake_level(i, j) = std::min(min_level(i, j), min_basin(i, j));
            }
          }
          MinMaxChanged = true;
        } else if (mask_ij == 2) {
          //Extend existing lake by new cells
          if ( gc.islake(old_ll(i, j)) ) {
            lake_level(i, j) = old_ll(i, j);
          } else {
            lake_level(i, j) = min_level(i, j);
          }
        }
      }
    } catch (...) {
      ParSec.failed();
    }
    ParSec.check();
  }

  lake_level.update_ghosts();

  return (GlobalOr(m_grid->com, MinMaxChanged));
}

void Gradual::gradually_fill(const double dt,
                             const double max_fill_rate,
                             const IceModelVec2S &target_level,
                             const IceModelVec2S &bed,
                             const IceModelVec2S &thk,
                             const IceModelVec2S &sea_level,
                             const IceModelVec2S &min_level,
                             const IceModelVec2S &max_level,
                             const IceModelVec2S &min_bed,
                             const IceModelVec2S &fill_rate,
                             IceModelVec2S &lake_level) {

  double dh_max = max_fill_rate * dt;

  GeometryCalculator gc(*m_config);

  IceModelVec::AccessList list{ &lake_level, &target_level,
                                &min_level, &max_level, &min_bed,
                                &sea_level, &bed, &thk };

  if (not m_use_const_fill_rate) {
    list.add(fill_rate);
  }

  //Update lakes
  ParallelSection ParSec(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      if (gc.islake(lake_level(i, j))) {
        const double min_ij = min_level(i, j),
                     max_ij = max_level(i, j),
                     current_ij = lake_level(i, j),
                     target_level_ij  = target_level(i, j);
        const bool unbalanced_level = (min_ij != max_ij),
                   disappear = not gc.islake(target_level_ij);

        double target_ij = target_level_ij;

        // if lakes vanish ...
        if (disappear) {
          // ... and become ocean, set to sea level, else to floatation level
          if (mask::ocean(gc.mask(sea_level(i, j), bed(i, j), thk(i, j)))) {
            target_ij = sea_level(i, j);
          } else {
            if (mask::ocean(gc.mask(sea_level(i, j), bed(i, j), thk(i, j), current_ij))) {
              // if "real" lake, gradually empty it to floatation level
              target_ij = bed(i, j) + m_alpha_lake * thk(i, j);
            } else {
              //else (if it is below floatation level) immediately get rid of it.
              target_ij = current_ij;
            }
          }
        }

        const bool rising = (current_ij < target_ij);

        if (not m_use_const_fill_rate or not rising) {
          dh_max = max_fill_rate * dt;
        }

        if (rising) {
          //Only if rising, adjust fill rate to mass loss of ice sheet
          if (not (m_use_const_fill_rate or unbalanced_level)) {
            const double fill_rate_ij = fill_rate(i, j);

            if (gc.islake(fill_rate_ij) and (fill_rate_ij <= max_fill_rate) and not disappear) {
              dh_max = fill_rate_ij * dt;
            } else {
              dh_max = max_fill_rate * dt;
            }
          }

          const double new_level = ( disappear ? current_ij : min_ij ) + std::min(dh_max, std::abs(target_ij - min_ij));

          if ( (new_level > current_ij) or disappear ) {
            if (disappear and (new_level >= target_ij)) {
              lake_level(i, j) = m_fill_value;
            } else {
              lake_level(i, j) = new_level;
            }
          }
        } else {
          const double new_level = ( disappear ? current_ij : max_ij ) - std::min(dh_max, std::abs(max_ij - target_ij));

          if ( (new_level < current_ij) or disappear ) {
            if (disappear and (new_level <= target_ij)) {
              lake_level(i, j) = m_fill_value;
            } else {
              lake_level(i, j) = new_level;
            }
          }
        }
      }
    }
  } catch (...) {
    ParSec.failed();
  }
  ParSec.check();
}

// Write diagnostic variables to extra files if requested
DiagnosticList Gradual::diagnostics_impl() const {

  DiagnosticList result = {
    { "lake_gradual_target",       Diagnostic::wrap(m_target_level) },
    { "lake_gradual_min_level",    Diagnostic::wrap(m_min_level) },
    { "lake_gradual_max_level",    Diagnostic::wrap(m_max_level) },
    { "lake_gradual_min_bed",      Diagnostic::wrap(m_min_basin) },
    { "lake_gradual_max_sl_basin", Diagnostic::wrap(m_max_sl_basin) },
    { "lake_expansion_mask",       Diagnostic::wrap(m_expansion_mask) },
    { "lake_area",                 Diagnostic::wrap(m_lake_area) },
    { "lake_mass_input_discharge", Diagnostic::wrap(m_lake_mass_input_discharge) },
    { "lake_mass_input_basal",     Diagnostic::wrap(m_lake_mass_input_basal) },
    { "lake_mass_input_total",     Diagnostic::wrap(m_lake_mass_input_total) },
    { "lake_fill_rate",            Diagnostic::wrap(m_lake_fill_rate) },
  };

  return combine(result, m_input_model->diagnostics());
}

} // end of namespace lake_level
} // end of namespace ocean
} // end of namespace pism
