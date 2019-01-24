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

#include "pism/frontretreat/FrontalMelt.hh"
#include "pism/frontretreat/IcebergRemover.hh"
#include "pism/frontretreat/calving/CalvingAtThickness.hh"
#include "pism/frontretreat/calving/EigenCalving.hh"
#include "pism/frontretreat/calving/FloatKill.hh"
#include "pism/frontretreat/calving/OceanKill.hh"
#include "pism/frontretreat/calving/vonMisesCalving.hh"

#include "pism/energy/EnergyModel.hh"
#include "pism/coupler/FrontalMeltModel.hh"
#include "pism/stressbalance/ShallowStressBalance.hh"

namespace pism {

void IceModel::do_calving() {

  {
    IceModelVec2S &flux_magnitude = m_work2d[0];

    flux_magnitude.set_to_magnitude(m_subglacial_hydrology->flux());

    FrontalMeltInputs inputs;
    inputs.geometry              = &m_geometry;
    inputs.subglacial_water_flux = &flux_magnitude;

    m_frontalmelt->update(inputs, m_time->current(), m_dt);
  }

  CalvingInputs inputs;

  inputs.geometry = &m_geometry;
  inputs.bc_mask  = &m_ssa_dirichlet_bc_mask;

  inputs.ice_velocity         = &m_stress_balance->shallow()->velocity();
  inputs.ice_enthalpy         = &m_energy_model->enthalpy();
  inputs.frontal_melt_rate    = &m_frontalmelt->frontal_melt_rate();

  // eigen-calving should go first: it uses the ice velocity field,
  // which is defined at grid points that were icy at the *beginning*
  // of a time-step.
  if (m_eigen_calving) {
    m_eigen_calving->update(m_dt,
                            inputs,
                            m_geometry.cell_type,
                            m_geometry.ice_area_specific_volume,
                            m_geometry.ice_thickness);
  }

  if (m_vonmises_calving) {
    m_vonmises_calving->update(m_dt,
                               inputs,
                               m_geometry.cell_type,
                               m_geometry.ice_area_specific_volume,
                               m_geometry.ice_thickness);
  }

  if (m_frontal_melt) {
    m_frontal_melt->update(m_dt,
                           inputs,
                           m_geometry.cell_type,
                           m_geometry.ice_area_specific_volume,
                           m_geometry.ice_thickness);
  }

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

/**
 * Compute the ice discharge into the ocean during the current time step.
 *
 * Units: ice equivalent meters.
 *
 * @param thickness current ice thickness
 * @param Href current "reference ice thickness"
 * @param thickness_old old ice thickness
 * @param Href_old old "reference ice thickness"
 * @param[in,out] output computed discharge during the current time step
 */
void IceModel::accumulate_discharge(const IceModelVec2S &thickness,
                                    const IceModelVec2S &Href,
                                    const IceModelVec2S &thickness_old,
                                    const IceModelVec2S &Href_old,
                                    IceModelVec2S &output) {

  IceModelVec::AccessList list{&thickness, &thickness_old,
      &Href, &Href_old, &output};

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    const double
      H_old       = thickness_old(i, j) + Href_old(i, j),
      H_new       = thickness(i, j) + Href(i, j);
    output(i, j) += H_new - H_old;
  }
}

} // end of namespace pism
