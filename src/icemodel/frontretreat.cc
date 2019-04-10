// Copyright (C) 2004--2019 Torsten Albrecht and Constantine Khroulev
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

#include "IceModel.hh"

#include "pism/util/IceGrid.hh"
#include "pism/util/Mask.hh"
#include "pism/util/ConfigInterface.hh"
#include "pism/util/pism_utilities.hh"

#include "pism/frontretreat/FrontRetreat.hh"
#include "pism/frontretreat/util/IcebergRemover.hh"
#include "pism/frontretreat/calving/CalvingAtThickness.hh"
#include "pism/frontretreat/calving/EigenCalving.hh"
#include "pism/frontretreat/calving/FloatKill.hh"
#include "pism/frontretreat/calving/OceanKill.hh"
#include "pism/frontretreat/calving/vonMisesCalving.hh"

#include "pism/energy/EnergyModel.hh"
#include "pism/coupler/FrontalMelt.hh"
#include "pism/stressbalance/ShallowStressBalance.hh"
#include "pism/hydrology/Hydrology.hh"
#include "pism/frontretreat/util/remove_narrow_tongues.hh"

namespace pism {

void IceModel::front_retreat_step() {
  IceModelVec2S
    &old_H    = m_work2d[0],
    &old_Href = m_work2d[1];

  {
    old_H.copy_from(m_geometry.ice_thickness);
    old_Href.copy_from(m_geometry.ice_area_specific_volume);
    m_thickness_change.calving.set(0.0);
  }

  {
    // mechanisms that use a retreat rate
    if (m_eigen_calving or m_vonmises_calving or m_frontal_melt) {
      // at least one of front retreat mechanisms is active

      if (m_frontal_melt) {

        IceModelVec2S &flux_magnitude = m_work2d[3];

        flux_magnitude.set_to_magnitude(m_subglacial_hydrology->flux());

        FrontalMeltInputs inputs;

        inputs.geometry = &m_geometry;
        inputs.subglacial_water_flux = &flux_magnitude;

        m_frontal_melt->update(inputs, m_time->current(), m_dt);

        // apply frontal melt rate
        // retreat_rate.add(1.0, m_frontal_melt->retreat_rate());

        compute_geometry_change(m_geometry.ice_thickness,
                                m_geometry.ice_area_specific_volume,
                                old_H, old_Href,
                                REPLACE_CHANGES,
                                m_thickness_change.frontal_melt);

      } else {
        m_thickness_change.frontal_melt.set(0.0);
      }

      IceModelVec2S &retreat_rate = m_work2d[2];
      retreat_rate.set(0.0);

      if (m_eigen_calving) {
        m_eigen_calving->update(m_geometry.cell_type,
                                m_stress_balance->shallow()->velocity());
        retreat_rate.add(1.0, m_eigen_calving->calving_rate());
      }

      if (m_vonmises_calving) {
        // FIXME: consider computing vertically-averaged hardness here and providing that
        // instead of using ice thickness and enthalpy.
        m_vonmises_calving->update(m_geometry.cell_type,
                                   m_geometry.ice_thickness,
                                   m_stress_balance->shallow()->velocity(),
                                   m_energy_model->enthalpy());
        retreat_rate.add(1.0, m_vonmises_calving->calving_rate());
      }


      assert(m_front_retreat);

      m_front_retreat->update_geometry(m_dt, m_geometry, m_ssa_dirichlet_bc_mask,
                                       retreat_rate,
                                       m_geometry.ice_area_specific_volume,
                                       m_geometry.ice_thickness);

      auto thickness_threshold = m_config->get_double("stress_balance.ice_free_thickness_standard");

      m_geometry.ensure_consistency(thickness_threshold);

      if (m_eigen_calving or m_vonmises_calving) {
        remove_narrow_tongues(m_geometry.cell_type, m_geometry.ice_thickness);

        m_geometry.ensure_consistency(thickness_threshold);
      }
    }

    // calving mechanisms that remove ice at a given location instead of applying a retreat
    // rate
    {
      if (m_ocean_kill_calving) {
        m_ocean_kill_calving->update(m_geometry.cell_type, m_geometry.ice_thickness);
      }

      if (m_float_kill_calving) {
        m_float_kill_calving->update(m_geometry.cell_type, m_geometry.ice_thickness);
      }

      if (m_thickness_threshold_calving) {
        m_thickness_threshold_calving->update(m_geometry.cell_type, m_geometry.ice_thickness);
      }
    }
  }

  enforce_consistency_of_geometry(REMOVE_ICEBERGS);

  // clean up partially-filled cells that are not next to ice
  {
    IceModelVec::AccessList list{&m_geometry.ice_area_specific_volume,
                                 &m_geometry.cell_type};

    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      if (m_geometry.ice_area_specific_volume(i, j) > 0.0 and
          not m_geometry.cell_type.next_to_ice(i, j)) {
        m_geometry.ice_area_specific_volume(i, j) = 0.0;
      }
    }
  }

  compute_geometry_change(m_geometry.ice_thickness,
                          m_geometry.ice_area_specific_volume,
                          old_H, old_Href,
                          ADD_CHANGES,
                          m_thickness_change.calving);
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
 * @param[in,out] output computed change
 */
void IceModel::compute_geometry_change(const IceModelVec2S &thickness,
                                       const IceModelVec2S &Href,
                                       const IceModelVec2S &thickness_old,
                                       const IceModelVec2S &Href_old,
                                       GeometryChangeFlag flag,
                                       IceModelVec2S &output) {

  IceModelVec::AccessList list{&thickness, &thickness_old,
      &Href, &Href_old, &output};

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    const double
      H_old  = thickness_old(i, j) + Href_old(i, j),
      H_new  = thickness(i, j) + Href(i, j),
      change = H_new - H_old;

    if (flag == ADD_CHANGES) {
      output(i, j) += change;
    } else {
      output(i, j) = change;
    }
  }
}

} // end of namespace pism
