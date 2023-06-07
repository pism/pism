// Copyright (C) 2004--2021, 2023 Torsten Albrecht and Constantine Khroulev
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

#include "pism/icemodel/IceModel.hh"

#include "pism/util/ConfigInterface.hh"
#include "pism/util/Grid.hh"

#include "pism/frontretreat/FrontRetreat.hh"
#include "pism/frontretreat/calving/CalvingAtThickness.hh"
#include "pism/frontretreat/calving/EigenCalving.hh"
#include "pism/frontretreat/calving/FloatKill.hh"
#include "pism/frontretreat/calving/HayhurstCalving.hh"
#include "pism/frontretreat/calving/vonMisesCalving.hh"

#include "pism/energy/EnergyModel.hh"
#include "pism/coupler/FrontalMelt.hh"
#include "pism/stressbalance/ShallowStressBalance.hh"
#include "pism/hydrology/Hydrology.hh"
#include "pism/frontretreat/util/remove_narrow_tongues.hh"
#include "pism/frontretreat/PrescribedRetreat.hh"
#include "pism/util/ScalarForcing.hh"
#include "pism/util/label_components.hh"

namespace pism {

void IceModel::identify_open_ocean(const array::CellType &cell_type, array::Scalar &result) {

  auto &tmp_p0 = *m_work2d_proc0;

  array::AccessScope list{ &cell_type, &result };

  auto grid = cell_type.grid();

  // assume that ice-free ocean points at the edge of the domain belong to the "global
  // ocean"
  for (auto p = grid->points(); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (cell_type.ice_free_water(i, j)) {
      result(i, j) = 1.0;

      if (grid::domain_edge(*grid, i, j)) {
        result(i, j) = 2.0;
      }
    } else {
      result(i, j) = 0.0;
    }
  }

  label_components(result, tmp_p0, true, 2);

  // now `result` contains ones in "ice free ocean" cells that are not connected to the edge
  // of the domain and zeros elsewhere

  // create a mask that contains ones at "ice free ocean" locations connected to the edge
  // of the domain and zeros elsewhere:
  for (auto p = grid->points(); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (cell_type.ice_free_water(i, j) and result(i, j) < 0.5) {
      result(i, j) = 1;
    } else {
      result(i, j) = 0;
    }
  }

  result.update_ghosts();
}

void IceModel::front_retreat_step() {

  bool retreat_rate_based_calving = m_eigen_calving or m_vonmises_calving or m_hayhurst_calving;
  bool calving_is_active =
      retreat_rate_based_calving or m_float_kill_calving or m_thickness_threshold_calving;
  bool frontal_melt_only_open_ocean = m_config->get_flag("frontal_melt.open_ocean_margins_only");

  auto &open_ocean_mask = *m_work2d[3];
  if (calving_is_active or (m_frontal_melt and frontal_melt_only_open_ocean)) {
    identify_open_ocean(m_geometry.cell_type, open_ocean_mask);
  } else {
    open_ocean_mask.set(1.0);
  }

  // compute retreat rates due to eigencalving, von Mises calving, Hayhurst calving,
  // and frontal melt.
  // We do this first to make sure that all three mechanisms use the same ice geometry.
  {
    if (m_eigen_calving) {
      m_eigen_calving->update(m_geometry.cell_type,
                              m_stress_balance->shallow()->velocity());
    }

    if (m_hayhurst_calving) {
      m_hayhurst_calving->update(m_geometry.cell_type,
                                 m_geometry.ice_thickness,
                                 m_geometry.sea_level_elevation,
                                 m_geometry.bed_elevation);
    }

    if (m_vonmises_calving) {
      // FIXME: consider computing vertically-averaged hardness here and providing that
      // instead of using ice thickness and enthalpy.
      m_vonmises_calving->update(m_geometry.cell_type,
                                 m_geometry.ice_thickness,
                                 m_stress_balance->shallow()->velocity(),
                                 m_energy_model->enthalpy());
    }

    if (m_frontal_melt) {
      array::Scalar &flux_magnitude = *m_work2d[0];

      compute_magnitude(m_subglacial_hydrology->flux(), flux_magnitude);

      FrontalMeltInputs inputs;

      inputs.geometry = &m_geometry;
      inputs.subglacial_water_flux = &flux_magnitude;

      m_frontal_melt->update(inputs, m_time->current(), m_dt);
    }
  }

  array::Scalar
    &old_H    = *m_work2d[0],
    &old_Href = *m_work2d[1];

  // frontal melt
  if (m_frontal_melt) {
    assert(m_front_retreat);

    old_H.copy_from(m_geometry.ice_thickness);
    old_Href.copy_from(m_geometry.ice_area_specific_volume);

    array::Scalar &retreat_rate = *m_work2d[2];
    retreat_rate.copy_from(m_frontal_melt->retreat_rate());

    if (frontal_melt_only_open_ocean) {
      array::AccessScope list{ &retreat_rate, &open_ocean_mask };

      for (Points p(*m_grid); p; p.next()) {
        const int i = p.i(), j = p.j();

        if (open_ocean_mask(i, j) < 0.5) {
          retreat_rate(i, j) = 0.0;
        }
      }
    }

    // apply the frontal melt rate
    m_front_retreat->update_geometry(m_dt, m_geometry, m_ice_thickness_bc_mask,
                                     retreat_rate,
                                     m_geometry.ice_area_specific_volume,
                                     m_geometry.ice_thickness);
    bool add_values = false;
    compute_geometry_change(m_geometry.ice_thickness,
                            m_geometry.ice_area_specific_volume,
                            old_H, old_Href,
                            add_values,
                            m_thickness_change.frontal_melt);
  } else {
    m_thickness_change.frontal_melt.set(0.0);
  }

  // calving
  if (calving_is_active) {

    old_H.copy_from(m_geometry.ice_thickness);
    old_Href.copy_from(m_geometry.ice_area_specific_volume);

    // retreat-rate-based calving parameterizations:
    if (retreat_rate_based_calving) {
      assert(m_front_retreat);

      array::Scalar &retreat_rate = *m_work2d[2];
      retreat_rate.set(0.0);

      if (m_eigen_calving) {
        retreat_rate.add(1.0, m_eigen_calving->calving_rate());
      }

      if (m_hayhurst_calving) {
        retreat_rate.add(1.0, m_hayhurst_calving->calving_rate());
      }

      if (m_vonmises_calving) {
        retreat_rate.add(1.0, m_vonmises_calving->calving_rate());
      }

      if (m_calving_rate_factor) {
        double T = m_time->current() + 0.5 * m_dt;
        retreat_rate.scale(m_calving_rate_factor->value(T));
      }

      // modify the retreat rate to avoid calving into "holes" in ice shelves that are not
      // connected to the open ocean
      {
        array::AccessScope list{ &open_ocean_mask, &retreat_rate };

        for (Points p(*m_grid); p; p.next()) {
          const int i = p.i(), j = p.j();

          if (open_ocean_mask(i, j) < 0.5) {
            retreat_rate(i, j) = 0.0;
          }
        }
      }

      m_front_retreat->update_geometry(m_dt, m_geometry, m_ice_thickness_bc_mask,
                                       retreat_rate,
                                       m_geometry.ice_area_specific_volume,
                                       m_geometry.ice_thickness);

      auto thickness_threshold = m_config->get_number("stress_balance.ice_free_thickness_standard");

      m_geometry.ensure_consistency(thickness_threshold);

      if (m_eigen_calving or m_vonmises_calving or m_hayhurst_calving) {
        remove_narrow_tongues(m_geometry, m_geometry.ice_thickness);

        m_geometry.ensure_consistency(thickness_threshold);
      }
    }

    // calving using local geometry (usually calving one grid cell per time step)
    {
      auto &modified_cell_type = *m_work2d[2];

      // create a modified cell type mask to avoid calving into holes in ice shelves that
      // are not connected to the open ocean
      {
        const auto &cell_type = m_geometry.cell_type;

        array::AccessScope list{ &modified_cell_type, &cell_type, &open_ocean_mask };

        for (Points p(*m_grid); p; p.next()) {
          const int i = p.i(), j = p.j();

          if (cell_type.ice_free_water(i, j) and open_ocean_mask(i, j) < 0.5) {
            // This modification will ensure that cells next to *these* ice free ocean
            // cells will not be considered "marginal" and so thickness threshold and
            // float-kill parameterizations will not apply.
            modified_cell_type(i, j) = cell_type::UNKNOWN;
          } else {
            modified_cell_type(i, j) = cell_type(i, j);
          }
        }
      }

      if (m_float_kill_calving) {
        m_float_kill_calving->update(modified_cell_type, m_geometry.ice_thickness);
      }

      if (m_thickness_threshold_calving) {
        m_thickness_threshold_calving->update(m_time->current(), m_dt, modified_cell_type,
                                              m_geometry.ice_thickness);
      }
    }

    bool add_values = false;
    compute_geometry_change(m_geometry.ice_thickness,
                            m_geometry.ice_area_specific_volume,
                            old_H, old_Href,
                            add_values,
                            m_thickness_change.calving);
  } else {
    m_thickness_change.calving.set(0.0);
  }

  // prescribed retreat

  if (m_prescribed_retreat) {
    old_H.copy_from(m_geometry.ice_thickness);
    old_Href.copy_from(m_geometry.ice_area_specific_volume);

    m_prescribed_retreat->update(m_time->current(), m_dt,
                                 m_geometry.ice_thickness,
                                 m_geometry.ice_area_specific_volume);

    bool add_values = false;
    compute_geometry_change(m_geometry.ice_thickness,
                            m_geometry.ice_area_specific_volume,
                            old_H, old_Href,
                            add_values,
                            m_thickness_change.forced_retreat);

  } else {
    m_thickness_change.forced_retreat.set(0.0);
  }

  // Changes above may create icebergs; here we remove them and account for additional
  // mass losses.
  {
    old_H.copy_from(m_geometry.ice_thickness);
    old_Href.copy_from(m_geometry.ice_area_specific_volume);

    enforce_consistency_of_geometry(REMOVE_ICEBERGS);

    bool add_values = true;
    compute_geometry_change(m_geometry.ice_thickness,
                            m_geometry.ice_area_specific_volume,
                            old_H, old_Href,
                            add_values,
                            m_thickness_change.calving);
  }
}

/**
 * Compute the change in ice geometry from "old" to "current".
 *
 * Units: ice equivalent meters.
 *
 * @param thickness current ice thickness
 * @param Href current "reference ice thickness"
 * @param thickness_old old ice thickness
 * @param Href_old old "reference ice thickness"
 * @param[in] add_values if `true`, add computed values to `output`, otherwise
 *            overwrite them
 * @param[in,out] output computed change
 */
void IceModel::compute_geometry_change(const array::Scalar &thickness,
                                       const array::Scalar &Href,
                                       const array::Scalar &thickness_old,
                                       const array::Scalar &Href_old,
                                       bool add_values,
                                       array::Scalar &output) {

  array::AccessScope list{&thickness, &thickness_old,
      &Href, &Href_old, &output};

  for (auto p = m_grid->points(); p; p.next()) {
    const int i = p.i(), j = p.j();

    const double
      H_old  = thickness_old(i, j) + Href_old(i, j),
      H_new  = thickness(i, j) + Href(i, j),
      change = H_new - H_old;

    if (add_values) {
      output(i, j) += change;
    } else {
      output(i, j) = change;
    }
  }
}

} // end of namespace pism
