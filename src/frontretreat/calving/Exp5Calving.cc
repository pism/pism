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
#include "pism/coupler/util/options.hh"

#include "pism/util/array/CellType.hh"
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

  if (m_calving_threshold <= 0.0) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "'calving.exp5_calving.threshold' has to be positive");
  }
    
  m_log->message(2, "  Exp5 thickness threshold: %3.3f meters.\n", m_calving_threshold);

  m_calving_along_flow = m_config->get_flag("calving.exp5_calving.calve_along_flow_direction");

  if (m_calving_along_flow) {
    m_log->message(2, "  Exp5 calving along terminal ice flow.\n");
  }

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


  m_log->message(3, "    Update Exp5 calving rate.\n");

  array::AccessScope list{&cell_type, &ice_thickness, &m_calving_rate, &ice_velocity};

  // a shortcut for readability:
  const auto &Hc = m_calving_threshold;
  double C = convert(m_sys, 1.0, "m year-1", "m second-1");


  if (m_calving_along_flow) {

    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      m_calving_rate(i, j) = 0.0;

      if (cell_type.ice_free_ocean(i, j) and cell_type.next_to_floating_ice(i, j)) {

        //m_log->message(3, "  Exp5 at %d,%d: %3.1f, %3.3f, %3.3f\n",i,j,ice_thickness(i,j),ice_velocity(i, j).u,ice_velocity(i, j).v);
        if (ice_velocity(i-1, j).u > 0.0 and ice_thickness(i-1,j) > 0.0) {
          m_calving_rate(i, j) += std::max(0.0, 1.0 + (Hc - ice_thickness(i-1,j)) / Hc) * ice_velocity(i-1, j).u;
          //m_log->message(3, "  Exp5 n: %3.1f, %3.3f, %3.3f at %d,%d.\n", ice_thickness(i-1,j),ice_velocity(i-1, j).u/C,m_calving_rate(i, j)/C,i,j);
        }
        if (ice_velocity(i+1, j).u < 0.0 and ice_thickness(i+1,j) > 0.0) {
          m_calving_rate(i, j) -= std::max(0.0, 1.0 + (Hc - ice_thickness(i+1,j)) / Hc) * ice_velocity(i+1, j).u;
          //m_log->message(3, "  Exp5 s: %3.1f, %3.3f, %3.3f at %d,%d.\n", ice_thickness(i+1,j),ice_velocity(i+1, j).u/C,m_calving_rate(i, j)/C,i,j);
        }
        if (ice_velocity(i, j-1).v > 0.0 and ice_thickness(i,j-1) > 0.0) {
          m_calving_rate(i, j) += std::max(0.0, 1.0 + (Hc - ice_thickness(i,j-1)) / Hc) * ice_velocity(i, j-1).v;
          //m_log->message(3, "  Exp5 w: %3.1f, %3.3f, %3.3f at %d,%d.\n", ice_thickness(i,j-1),ice_velocity(i, j-1).v/C,m_calving_rate(i, j)/C,i,j);
        }
        if (ice_velocity(i, j+1).v < 0.0 and ice_thickness(i,j+1) > 0.0) {
          m_calving_rate(i, j) -= std::max(0.0, 1.0 + (Hc - ice_thickness(i,j+1)) / Hc) * ice_velocity(i, j+1).v;
          //m_log->message(3, "  Exp5 e: %3.1f, %3.3f, %3.3f at %d,%d.\n", ice_thickness(i,j+1),ice_velocity(i, j+1).v/C,m_calving_rate(i, j)/C,i,j);
        }
        
        //m_calving_rate(i, j) /= C;
      }
    }
  } 
  else { //only magnitudes

    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      if (cell_type.floating_ice(i, j) && cell_type.next_to_ice_free_ocean(i, j)) {
        double Iv            = ice_velocity(i, j).magnitude();
        double H             = ice_thickness(i, j);
        m_calving_rate(i, j) = std::max(0.0, 1.0 + (Hc - H) / Hc) * Iv;
        //m_calving_rate(i, j) /= C;
      } else {
        m_calving_rate(i, j) = 0.0;
      }
    }
    m_calving_rate.update_ghosts();


    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      if (cell_type.ice_free_ocean(i, j) and cell_type.next_to_ice(i, j)) {

        auto R = m_calving_rate.star(i, j);
        auto M = cell_type.star_int(i, j);

        int N        = 0;
        double R_sum = 0.0;
        for (auto direction : { North, East, South, West }) {
          if (mask::icy(M[direction])) {
            R_sum += R[direction];
            N++;
          }
        }

        if (N > 0) {
          m_calving_rate(i, j) = R_sum / N;
        }
      }
    } 
  }
}


DiagnosticList Exp5Calving::diagnostics_impl() const {
  return {{"exp5_calving_rate", Diagnostic::wrap(m_calving_rate)}};
}

const array::Scalar &Exp5Calving::calving_rate() const {
  return m_calving_rate;
}



} // end of namespace calving
} // end of namespace pism
