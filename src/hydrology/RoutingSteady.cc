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

  ice_bottom_surface(*inputs.geometry, m_bottom_surface);

  m_Qstag_average.set(0.0);

  // make sure W has valid ghosts before starting hydrology steps
  m_W.update_ghosts();

  unsigned int step_counter = 0;
  for (; ht < t_final; ht += hdt) {
    step_counter++;

#if (PISM_DEBUG==1)
    double huge_number = 1e6;
    check_bounds(m_W, huge_number);

    check_bounds(m_Wtill, m_config->get_double("hydrology.tillwat_max"));
#endif

    // updates ghosts of m_Wstag
    water_thickness_staggered(m_W,
                              inputs.geometry->cell_type,
                              m_Wstag);

    double maxKW = 0.0;
    // updates ghosts of m_Kstag
    m_grid->ctx()->profiling().begin("routing_conductivity");
    compute_conductivity(m_Wstag,
                         subglacial_water_pressure(),
                         m_bottom_surface,
                         m_Kstag, maxKW);
    m_grid->ctx()->profiling().end("routing_conductivity");

    // ghosts of m_Vstag are not updated
    m_grid->ctx()->profiling().begin("routing_velocity");
    compute_velocity(m_Wstag,
                     subglacial_water_pressure(),
                     m_bottom_surface,
                     m_Kstag,
                     inputs.no_model_mask,
                     m_Vstag);
    m_grid->ctx()->profiling().end("routing_velocity");

    // to get Q, W needs valid ghosts (ghosts of m_Vstag are not used)
    // updates ghosts of m_Qstag
    m_grid->ctx()->profiling().begin("routing_flux");
    advective_fluxes(m_Vstag, m_W, m_Qstag);
    m_grid->ctx()->profiling().end("routing_flux");

    m_Qstag_average.add(hdt, m_Qstag);

    {
      const double
        dt_cfl    = max_timestep_W_cfl(),
        dt_diff_w = max_timestep_W_diff(maxKW);

      hdt = std::min(t_final - ht, dt_max);
      hdt = std::min(hdt, dt_cfl);
      hdt = std::min(hdt, dt_diff_w);
    }

    m_log->message(3, "  hydrology step %05d, dt = %f s\n", step_counter, hdt);

    // update Wtillnew from Wtill and input_rate
    {
      m_grid->ctx()->profiling().begin("routing_Wtill");
      update_Wtill(hdt,
                   m_Wtill,
                   m_surface_input_rate,
                   m_basal_melt_rate,
                   m_Wtillnew);
      // remove water in ice-free areas and account for changes
      enforce_bounds(inputs.geometry->cell_type,
                     inputs.no_model_mask,
                     0.0,        // do not limit maximum thickness
                     m_Wtillnew,
                     m_grounded_margin_change,
                     m_grounding_line_change,
                     m_conservation_error_change,
                     m_no_model_mask_change);
      m_grid->ctx()->profiling().end("routing_Wtill");
    }

    // update Wnew from W, Wtill, Wtillnew, Wstag, Q, input_rate
    // uses ghosts of m_W, m_Wstag, m_Qstag, m_Kstag
    {
      m_grid->ctx()->profiling().begin("routing_W");
      update_W(hdt,
               m_surface_input_rate,
               m_basal_melt_rate,
               m_W, m_Wstag,
               m_Wtill, m_Wtillnew,
               m_Kstag, m_Qstag,
               m_Wnew);
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
      m_grid->ctx()->profiling().end("routing_W");
    }

    // m_Wtill has no ghosts
    m_Wtill.copy_from(m_Wtillnew);
  } // end of the time-stepping loop

  staggered_to_regular(inputs.geometry->cell_type, m_Qstag_average,
                       m_config->get_boolean("hydrology.routing.include_floating_ice"),
                       m_Q);
  m_Q.scale(1.0 / dt);

  m_log->message(2,
                 "  took %d hydrology sub-steps with average dt = %.6f years (%.3f s or %.3f hours)\n",
                 step_counter,
                 units::convert(m_sys, dt / step_counter, "seconds", "years"),
                 dt / step_counter,
                 (dt / step_counter) / 3600.0);
}

void RoutingSteady::define_model_state_impl(const PIO& output) const {
  m_W.define(output);
}

void RoutingSteady::write_model_state_impl(const PIO& output) const {
  m_W.write(output);
}


} // end of namespace hydrology
} // end of namespace pism
