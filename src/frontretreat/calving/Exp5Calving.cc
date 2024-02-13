/* Copyright (C) 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2021, 2022, 2024 PISM Authors
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

#include "Exp5Calving.hh"

#include "pism/util/Mask.hh"
#include "pism/util/IceGrid.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/coupler/util/options.hh"

#include "pism/util/array/CellType.hh"
#include "pism/geometry/Geometry.hh"

#include "pism/util/array/Vector.hh"


namespace pism {

//! @brief Calving and iceberg removal code as in https://github.com/JRowanJordan/CalvingMIP/wiki/Experiment-5
namespace calving {


Exp5Calving::Exp5Calving(IceGrid::ConstPtr grid)
  : Component(grid),
    m_calving_rate(grid, "exp5_calving_rate")
{
  m_calving_rate.metadata().set_name("exp5_calving_rate");
  m_calving_rate.set_attrs("diagnostic",
                           "horizontal calving rate due to Exp5 calving",
                           "m s-1", "m year-1", "", 0);
}

void Exp5Calving::init() {

  m_log->message(2, "* Initializing the 'EXP5 calving' mechanism...\n");

  m_calving_threshold = m_config->get_number("calving.exp5_calving.threshold");
    
  m_log->message(2, "  Exp5 thickness threshold: %3.3f meters.\n", calving_threshold);
}

/**
 * Updates calving rate according to the
 * calving rule Exp5 removing ice at the shelf front 
 * Cr = max(0,1+(Hc-H)/Hc) * Iv;.
 *
 * @param[in] pism_mask ice cover mask
 * @param[in] ice_velocity ice velocity
 * @param[in] ice_thickness ice thickness
 *
 * @return 0 on success
 */

void Exp5Calving::update(const array::CellType1 &cell_type,
                         const array::Vector1 &ice_velocity,
                         const array::Scalar &ice_thickness) {


  m_log->message(3,
                   "    Update Exp5 calving rate.\n");

  m_calving_rate.set(0.0);

  array::AccessScope list{&cell_type, &ice_thickness, &m_calving_rate, &ice_velocity};

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (cell_type.floating_ice(i, j) && cell_type.next_to_ice_free_ocean(i, j)) {

      double terminal_velocity = ice_velocity(i, j).magnitude();
      double Hcr               = 1.0 + (m_calving_threshold - ice_thickness(i, j)) / m_calving_threshold;
      if (Hcr < 0.0) {
        m_calving_rate(i, j) = 0.0;
      } else {
        m_calving_rate(i, j) = Hcr * terminal_velocity;
      }
    }
  }
  m_calving_rate.update_ghosts();

  const Direction dirs[] = {North, East, South, West};

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (cell_type.ice_free_ocean(i, j) and cell_type.next_to_ice(i, j) ) {

      auto R = m_calving_rate.star(i, j);
      auto M = cell_type.star(i, j);

      int N = 0;
      double R_sum = 0.0;
      for (int n = 0; n < 4; ++n) {
        Direction direction = dirs[n];
        if (mask::icy(M[direction])) {
          R_sum += R[direction];
          N++;
        }
      }

      if (N > 0) {
        m_calving_rate(i, j) = R_sum / N;
      }

    m_calving_rate(i, j) *= convert(m_sys, 1.0, "m year-1", "m second-1");
    //m_log->message(2,
    //               "!!!!!!!!!!!!!!!!!! Exp5 cr=%f at %d,%d.\n",m_calving_rate(i, j),i,j);
    }
  }

    //        u_abs = fabs(ice_velocity(i, j).u),
    //        v_abs = fabs(ice_velocity(i, j).v);

}


DiagnosticList Exp5Calving::diagnostics_impl() const {
  return {{"exp5_calving_rate", Diagnostic::wrap(m_calving_rate)}};
}

const array::Scalar &Exp5Calving::calving_rate() const {
  return m_calving_rate;
}



} // end of namespace calving
} // end of namespace pism
