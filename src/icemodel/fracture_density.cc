// Copyright (C) 2011-2020, 2023 Torsten Albrecht and Constantine Khroulev
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 2 of the License, or (at your option) any later
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

#include "pism/energy/EnergyModel.hh"
#include "pism/stressbalance/ShallowStressBalance.hh"
#include "pism/fracturedensity/FractureDensity.hh"

namespace pism {

void IceModel::update_fracture_density() {
  // generate the BC mask for the fracture density model
  //
  // This mask contains ones at the in-flow boundary according to the SSA Dirichlet BC
  // mask (if it is used) and at grounded grid points if
  // fracture_density.include_grounded_ice is not set.
  auto &bc_mask = *m_work2d[0];
  {
    bool do_fracground = m_config->get_flag("fracture_density.include_grounded_ice");
    const bool dirichlet_bc = m_config->get_flag("stress_balance.ssa.dirichlet_bc");

    bc_mask.set(0.0);

    array::AccessScope list{&bc_mask, &m_geometry.cell_type};

    if (dirichlet_bc) {
      list.add(m_velocity_bc_mask);
      list.add(m_velocity_bc_values);
    }

    for (auto p = m_grid->points(); p; p.next()) {
      const int i = p.i(), j = p.j();

      if (m_geometry.cell_type.land(i, j) and not do_fracground) {
        bc_mask(i, j) = 1.0;
      }

      if (dirichlet_bc) {
        if (m_velocity_bc_mask(i, j) > 0.5 and
            (m_velocity_bc_values(i, j).u != 0.0 or
             m_velocity_bc_values(i, j).v != 0.0)) {
          bc_mask(i, j) = 1.0;
        }
      }
    } // end of the loop over grid points
  }

  // compute the vertically-averaged ice hardness
  auto &hardness = *m_work2d[1];
  {
    rheology::averaged_hardness_vec(*m_stress_balance->shallow()->flow_law(),
                                    m_geometry.ice_thickness,
                                    m_energy_model->enthalpy(),
                                    hardness);
  }

  // This model has the same time-step restriction as the mass transport code so we don't
  // check if this time step is short enough.
  m_fracture->update(m_dt, m_geometry,
                     m_stress_balance->shallow()->velocity(),
                     hardness, bc_mask);
}

} // end of namespace pism
