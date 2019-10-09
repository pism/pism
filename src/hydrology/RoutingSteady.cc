/* Copyright (C) 2019 PISM Authors
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

#include "RoutingSteady.hh"

#include "pism/geometry/Geometry.hh"
#include "pism/util/pism_utilities.hh"

namespace pism {
namespace hydrology {


RoutingSteady::RoutingSteady(IceGrid::ConstPtr grid)
  : Routing(grid) {
  m_Wtill.set(0.0);
}

RoutingSteady::~RoutingSteady() {
  // empty
}

void RoutingSteady::restart_impl(const PIO& input_file, int record) {
  m_W.read(input_file, record);

  regrid("Hydrology", m_W);
}

void RoutingSteady::bootstrap_impl(const PIO& input_file, const IceModelVec2S& ice_thickness) {
  (void) ice_thickness;

  double bwat_default = m_config->get_double("bootstrapping.defaults.bwat");
  m_W.regrid(input_file, OPTIONAL, bwat_default);

  regrid("Hydrology", m_W);

}

void RoutingSteady::init_impl(const IceModelVec2S& W_till,
                              const IceModelVec2S& W,
                              const IceModelVec2S& P) {
  (void) W_till;
  (void) P;

  m_W.copy_from(W);
}

void RoutingSteady::update_impl(double t, double dt, const Inputs& inputs) {
  (void) t;

  ice_bottom_surface(*inputs.geometry, m_bottom_surface);

  m_Qstag_average.set(0.0);

  // this model does not have till storage
  m_Wtill.set(0.0);
  m_Wtillnew.set(0.0);

  // add all the water in the beginning (this updates ghosts)
  m_W.add(dt, m_surface_input_rate);
  m_input_change.add(dt, m_surface_input_rate);

  const unsigned int n_iterations = m_config->get_double("hydrology.routing_steady.n_iterations");
  unsigned int step_counter = 0;
  for (; step_counter < n_iterations; ++step_counter) {

    // updates ghosts of m_Wstag
    water_thickness_staggered(m_W,
                              inputs.geometry->cell_type,
                              m_Wstag);

    // ghosts of m_Vstag are not updated
    compute_velocity(m_Wstag,
                     subglacial_water_pressure(),
                     m_bottom_surface,
                     inputs.no_model_mask,
                     m_Vstag);
    double hdt = max_timestep_W_cfl();

    // to get Q, W needs valid ghosts (ghosts of m_Vstag are not used)
    // updates ghosts of m_Qstag
    advective_fluxes(m_Vstag, m_W, m_Qstag);

    m_Qstag_average.add(hdt, m_Qstag);

    // update Wnew from W and Q
    // uses ghosts of m_Qstag
    {
      IceModelVec::AccessList list{&m_W, &m_Qstag, &m_Wnew, &m_flow_change};

      for (Points p(*m_grid); p; p.next()) {
        const int i = p.i(), j = p.j();

        auto q = m_Qstag.star(i, j);
        const double divQ = (q.e - q.w) / m_dx + (q.n - q.s) / m_dy;

        double dW = hdt * (- divQ);

        m_Wnew(i, j) = m_W(i, j) + dW;

        m_flow_change(i, j) += dW;
      }

      // remove water in ice-free areas and account for changes
      enforce_bounds(inputs.geometry->cell_type,
                     inputs.no_model_mask,
                     0.0,        // do not limit maximum thickness
                     m_Wnew,
                     m_grounded_margin_change,
                     m_grounding_line_change,
                     m_conservation_error_change,
                     m_no_model_mask_change);

      // transfer new into old (updates ghosts of m_W)
      m_W.copy_from(m_Wnew);
    }
  } // end of the loop

  staggered_to_regular(inputs.geometry->cell_type, m_Qstag_average,
                       m_config->get_boolean("hydrology.routing.include_floating_ice"),
                       m_Q);
  m_Q.scale(1.0 / dt);
}


//! Get the advection velocity V at the center of cell edges.
/*!
  Computes the advection velocity @f$\mathbf{V}@f$ on the staggered
  (edge-centered) grid.  If V = (u, v) in components then we have
  <code> result(i, j, 0) = u(i+1/2, j) </code> and
  <code> result(i, j, 1) = v(i, j+1/2) </code>

  The advection velocity is given by the formula

  @f[ \mathbf{V} = - \nabla (P + \rho_w g (b + W)) / \left| \nabla (P + \rho_w g (b + W)) \right| @f]

  where @f$ \mathbf{V} @f$ is the water velocity, @f$ P @f$ is the water
  pressure, and @f$ b @f$ is the bedrock elevation.

  If the corresponding staggered grid value of the water thickness is zero then that
  component of V is set to zero. This does not change the flux value (which would be zero
  anyway) but it does provide the correct max velocity in the CFL calculation. We assume
  bed has valid ghosts.
*/
void RoutingSteady::compute_velocity(const IceModelVec2Stag &W,
                                     const IceModelVec2S &pressure,
                                     const IceModelVec2S &bed,
                                     const IceModelVec2Int *no_model_mask,
                                     IceModelVec2Stag &result) const {
  IceModelVec2S &P = m_R;
  P.copy_from(pressure);  // yes, it updates ghosts

  IceModelVec::AccessList list{&P, &W, &bed, &result};

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (W(i, j, 0) > 0.0) {
      double
        W_x = (W(i + 1, j, 0) - W(i - 1, j, 0)) / (2.0 * m_dx),
        W_y = (W(i, j + 1, 0) - W(i, j - 1, 0)) / (2.0 * m_dy),
        P_x = (P(i + 1, j) - P(i, j)) / m_dx,
        b_x = (bed(i + 1, j) - bed(i, j)) / m_dx,
        P_y = (+ P(i + 1, j + 1) + P(i, j + 1)
               - P(i + 1, j - 1) - P(i, j - 1)) / (4.0 * m_dy),
        b_y = (+ bed(i + 1, j + 1) + bed(i, j + 1)
               - bed(i + 1, j - 1) - bed(i, j - 1)) / (4.0 * m_dy),
        u   = - (P_x + m_rg * (b_x + W_x)),
        v   = - (P_y + m_rg * (b_y + W_y)),
        S   = Vector2(u, v).magnitude();

      result(i, j, 0) = S > 0.0 ? u / S : 0.0;
    } else {
      result(i, j, 0) = 0.0;
    }

    if (W(i, j, 1) > 0.0) {
      double
        W_x = (W(i + 1, j, 1) - W(i - 1, j, 1)) / (2.0 * m_dx),
        W_y = (W(i, j + 1, 1) - W(i, j - 1, 1)) / (2.0 * m_dy),
        P_x = (+ P(i + 1, j + 1) + P(i + 1, j)
               - P(i - 1, j + 1) - P(i - 1, j)) / (4.0 * m_dx),
        b_x = (+ bed(i + 1, j + 1) + bed(i + 1, j)
               - bed(i - 1, j + 1) - bed(i - 1, j)) / (4.0 * m_dx),
        P_y = (P(i, j + 1) - P(i, j)) / m_dy,
        b_y = (bed(i, j + 1) - bed(i, j)) / m_dy,
        u   = - (P_x + m_rg * (b_x + W_x)),
        v   = - (P_y + m_rg * (b_y + W_y)),
        S   = Vector2(u, v).magnitude();

      result(i, j, 1) = S > 0.0 ? v / S : 0.0;
    } else {
      result(i, j, 1) = 0.0;
    }
  }

  if (no_model_mask) {
    list.add(*no_model_mask);

    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      auto M = no_model_mask->int_star(i, j);

      if (M.ij or M.e) {
        result(i, j, 0) = 0.0;
      }

      if (M.ij or M.n) {
        result(i, j, 1) = 0.0;
      }
    }
  }
}

void RoutingSteady::define_model_state_impl(const PIO& output) const {
  m_W.define(output);
}

void RoutingSteady::write_model_state_impl(const PIO& output) const {
  m_W.write(output);
}


} // end of namespace hydrology
} // end of namespace pism
