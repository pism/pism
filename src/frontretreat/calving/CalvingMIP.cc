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

#include "CalvingMIP.hh"

#include "pism/util/Mask.hh"
#include "pism/util/IceGrid.hh"
#include "pism/coupler/util/options.hh"

#include "pism/util/array/CellType.hh"
#include "pism/util/array/Vector.hh"

#include "pism/util/Time.hh"

#include <cmath>                // pow, tan, atan


namespace pism {

//! @brief Calving and iceberg removal code as in https://github.com/JRowanJordan/CalvingMIP/wiki
namespace calving {


CalvingMIP::CalvingMIP(IceGrid::ConstPtr grid)
  : Component(grid),
    m_calving_rate(grid, "calvingmip_calving_rate"),
    m_calving_rate_tmp(grid, "calvingmip_calving_rate_temporary"),
    m_cell_type(grid, "cell_type_where_calving_rate_set")
{
  m_calving_rate.metadata().set_name("calvingmip_calving_rate");
  m_calving_rate.set_attrs("diagnostic",
                           "horizontal calving rate due to CalvingMIP calving",
                           "m s-1", "m year-1", "", 0);

  m_cell_type.set_attrs("internal", "cell type mask", "", "", "", 0);
}

void CalvingMIP::init() {

  m_log->message(2, "* Initializing the 'CalvingMIP calving' mechanism...\n");

  m_experiment = m_config->get_number("calving.calvingmip_calving.experiment");

  m_retreat_and_advance = false;
  if (m_experiment == 2 or m_experiment == 4)
    m_retreat_and_advance = true;

  if (m_experiment == 5) {
    m_calving_threshold = m_config->get_number("calving.calvingmip_calving.threshold");
    if (m_calving_threshold <= 0.0) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION, "'calving.calvingmip_calving.threshold' has to be positive");
    }
    m_log->message(2, "  CalvingMIP thickness threshold: %3.3f meters.\n", m_calving_threshold);
  }

  m_calving_along_flow = m_config->get_flag("calving.calvingmip_calving.calve_along_flow_direction");
  if (m_calving_along_flow) {
    m_log->message(2, "  CalvingMIP calving along terminal ice flow.\n");
  }

}


/**
 * Updates calving rate according to the
 * CalvingMIP calving rules removing ice at the shelf front
 * Cr = - Iv;.
 *
 * @param[in] pism_mask ice cover mask
 * @param[in] ice_velocity ice velocity
 * @param[in] ice_thickness ice thickness
 *
 * @return 0 on success
 */

void CalvingMIP::update(const array::CellType1 &cell_type,
                         const array::Vector1 &ice_velocity,
                         const array::Scalar &ice_thickness) {


  bool exp5 = false;
  if (m_experiment == 5)
    exp5 = true;

  m_cell_type.set(0);

  // a shortcut for readability:
  const auto &Hc = m_calving_threshold;
  double C = convert(m_sys, 1.0, "m year-1", "m second-1");
  // lower calving rate threshold
  double Cmin = 1.0e-7;
  // lower threshold for terminal velocity
  double vcr = 1e-20;
  double Wv = 0.0,
         thistime = m_grid->ctx()->time()->current(),
         starttime = m_grid->ctx()->time()->start(),
         thisyear = convert(m_sys, thistime-starttime, "seconds","year");

  if (m_retreat_and_advance) {
    if (m_experiment == 2)
      Wv = -300.0 * sin(2.0*M_PI*thisyear/1000.0);
    else if (m_experiment == 4)
      Wv = -750.0 * sin(2.0*M_PI*thisyear/1000.0);
    m_log->message(2, "    Update CalvingMIP Exp%d calving rate at year %f with Wv=%f m/yr\n",m_experiment,thisyear,Wv);
  } else {
    m_log->message(3, "    Update CalvingMIP Exp%d calving rate.\n",m_experiment);
  }

  array::AccessScope list{&cell_type, &m_cell_type, &ice_thickness, &m_calving_rate, &m_calving_rate_tmp, &ice_velocity};


  if (m_calving_along_flow) { // Here, the marginal ice velocity vectors are considered in calving rate

    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      double old_calv_rate = m_calving_rate(i, j);
      m_calving_rate(i, j) = 0.0;

      // This is relevant in case of a propagating ice shelf front, for which no terminal velocities exist yet
      if (cell_type.floating_ice(i, j) and cell_type.next_to_ice_free_ocean(i, j)) {
        if (ice_velocity(i, j).magnitude() < vcr) {
          m_calving_rate(i, j) = old_calv_rate;
          m_cell_type(i, j) = 1;
        }
      }

      if (cell_type.ice_free_ocean(i, j) and cell_type.next_to_floating_ice(i, j)) {

        double vw  = ice_velocity(i-1, j).u,
              vwm = ice_velocity(i-1, j).magnitude(),
              hw  = ice_thickness(i-1,j),
              ve  = ice_velocity(i+1, j).u,
              vem = ice_velocity(i+1, j).magnitude(),
              he  = ice_thickness(i+1,j),
              vs  = ice_velocity(i, j-1).v,
              vsm = ice_velocity(i, j-1).magnitude(),
              hs  = ice_thickness(i,j-1),
              vn  = ice_velocity(i, j+1).v,
              vnm = ice_velocity(i, j+1).magnitude(),
              hn  = ice_thickness(i,j+1);


        bool divide_by_null = false;

        // west
        if (vw > 0.0 and hw > 0.0) {
          if (exp5) {
            m_calving_rate(i, j) += std::max(0.0, 1.0 + (Hc - hw) / Hc) * vw;
          } else {
            m_calving_rate(i, j) += vw * (1.0 - (vwm > 0.0 ? Wv*C/vwm : 0.0));
          }
          if (vwm==0.0) divide_by_null=true;
        }
        // east
        if (ve < 0.0 and he > 0.0) {
          if (exp5) {
            m_calving_rate(i, j) -= std::max(0.0, 1.0 + (Hc - he) / Hc) * ve;
          } else {
            m_calving_rate(i, j) -= ve * (1.0 - (vem > 0.0 ? Wv*C/vem : 0.0));
          }
          if (vem==0.0) divide_by_null=true;
        }
        // south
        if (vs > 0.0 and hs > 0.0) {
          if (exp5) {
            m_calving_rate(i, j) += std::max(0.0, 1.0 + (Hc - hs) / Hc) * vs;
          } else {
            m_calving_rate(i, j) += vs * (1.0 - (vsm > 0.0 ? Wv*C/vsm : 0.0));
          }
          if (vsm==0.0) divide_by_null=true;
        }
        // north
        if (vn < 0.0 and hn > 0.0) {
          if (exp5) {
            m_calving_rate(i, j) -= std::max(0.0, 1.0 + (Hc - hn) / Hc) * vn;
          } else {
            m_calving_rate(i, j) -= vn * (1.0 - (vnm > 0.0 ? Wv*C/vnm : 0.0));
          }
          if (vnm==0.0) divide_by_null=true;
        }

        m_calving_rate(i, j) = std::max(0.0,m_calving_rate(i, j));
        m_cell_type(i, j) = 1;

        if (divide_by_null) {
          m_calving_rate(i, j) = old_calv_rate;
        }
      } // end of next to front
    } //end of i,j loop

    m_calving_rate.update_ghosts();
    m_cell_type.update_ghosts();
    m_calving_rate_tmp.copy_from(m_calving_rate);


    // This part shall fill up ocean cells next to calving front with mean calving rates,
    // for cases, when calving front propagates forward with once cell per adaptive timestep (CFL)
    for (Points p(*m_grid); p; p.next()) {
        const int i = p.i(), j = p.j();

        // where calving rate was defined in previous step
        bool next_to_calving_front = (m_cell_type(i+1, j)==1 or m_cell_type(i-1, j)==1 or m_cell_type(i, j+1)==1 or m_cell_type(i, j-1)==1);

        // this yields a band of three cell widths around calving front
        if (m_calving_rate(i, j) < Cmin and next_to_calving_front) {
          auto R = m_calving_rate.star(i, j);

          int N        = 0;
          double R_sum = 0.0;
          for (auto direction : { North, East, South, West }) {
            if (R[direction] > Cmin) {
              R_sum += R[direction];
              N++;
            }
          }

          if (N > 0) {
            m_calving_rate_tmp(i, j) = R_sum / N;
          }
        }
    }
    m_calving_rate.copy_from(m_calving_rate_tmp);
  } // end of calving_along_directions
  else { //only velocity magnitudes are considered in calving rate

    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      if (cell_type.floating_ice(i, j) && cell_type.next_to_ice_free_ocean(i, j)) {
        double Iv            = ice_velocity(i, j).magnitude();
        double H             = ice_thickness(i, j);
        if (exp5) {
          m_calving_rate(i, j) = std::max(0.0, 1.0 + (Hc - H) / Hc) * Iv;
        } else {
          m_calving_rate(i, j) = Iv - Wv;
        }
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


DiagnosticList CalvingMIP::diagnostics_impl() const {
  return {{"calvingmip_calving_rate", Diagnostic::wrap(m_calving_rate)}};
}

const array::Scalar &CalvingMIP::calving_rate() const {
  return m_calving_rate;
}



} // end of namespace calving
} // end of namespace pism
