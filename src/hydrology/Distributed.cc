// Copyright (C) 2012-2018 PISM Authors
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

#include <algorithm>            // std::min, std::max

#include "Distributed.hh"
#include "pism/util/Mask.hh"
#include "pism/util/Vars.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/io/PIO.hh"
#include "pism/util/pism_options.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/util/IceModelVec2CellType.hh"

namespace pism {
namespace hydrology {

Distributed::Distributed(IceGrid::ConstPtr g)
  : Routing(g) {

  // additional variables beyond hydrology::Routing
  m_P.create(m_grid, "bwp", WITH_GHOSTS, 1);
  m_P.set_attrs("model_state",
                "pressure of transportable water in subglacial layer",
                "Pa", "");
  m_P.metadata().set_double("valid_min", 0.0);

  m_Pnew.create(m_grid, "Pnew_internal", WITHOUT_GHOSTS);
  m_Pnew.set_attrs("internal",
                   "new transportable subglacial water pressure during update",
                   "Pa", "");
  m_Pnew.metadata().set_double("valid_min", 0.0);
}

Distributed::~Distributed() {
  // empty
}

void Distributed::initialization_message() const {
  m_log->message(2,
                 "* Initializing the distributed, linked-cavities subglacial hydrology model...\n");
}

void Distributed::restart_impl(const PIO &input_file, int record) {
  Routing::restart_impl(input_file, record);

  m_P.read(input_file, record);

  regrid("Hydrology", m_P);
}

void Distributed::bootstrap_impl(const PIO &input_file,
                                 const IceModelVec2S &ice_thickness) {
  Routing::bootstrap_impl(input_file, ice_thickness);

  double bwp_default = m_config->get_double("bootstrapping.defaults.bwp");
  m_P.regrid(input_file, OPTIONAL, bwp_default);

  regrid("Hydrology", m_P);

  bool init_P_from_steady = m_config->get_boolean("hydrology.distributed.init_p_from_steady");

  if (init_P_from_steady) { // if so, just overwrite -i or -bootstrap value of P=bwp
    m_log->message(2,
                   "  initializing P from P(W) formula which applies in steady state\n");

    compute_overburden_pressure(ice_thickness, m_Pover);

    IceModelVec2S sliding_speed;
    sliding_speed.create(m_grid, "velbase_mag", WITHOUT_GHOSTS);
    sliding_speed.set_attrs("internal", "basal sliding speed", "m s-1", "");

    std::string filename = m_config->get_string("hydrology.distributed.sliding_speed_file");

    if (filename.empty()) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                    "hydrology.distributed.sliding_speed_file is not set");
    }

    sliding_speed.regrid(filename, CRITICAL);

    P_from_W_steady(m_W, m_Pover, sliding_speed,
                    m_P);
  }
}

void Distributed::initialize_impl(const IceModelVec2S &W_till,
                              const IceModelVec2S &W,
                              const IceModelVec2S &P) {
  Routing::initialize_impl(W_till, W, P);

  m_P.copy_from(P);
}

void Distributed::define_model_state_impl(const PIO &output) const {
  Routing::define_model_state_impl(output);
  m_P.define(output);
}

void Distributed::write_model_state_impl(const PIO &output) const {
  Routing::write_model_state_impl(output);
  m_P.write(output);
}

std::map<std::string, TSDiagnostic::Ptr> Distributed::ts_diagnostics_impl() const {
  std::map<std::string, TSDiagnostic::Ptr> result = {
    // FIXME: add mass-conservation time-series diagnostics
  };
  return result;
}

//! Copies the P state variable which is the modeled water pressure.
const IceModelVec2S& Distributed::subglacial_water_pressure() const {
  return m_P;
}

//! Check bounds on P and fail with message if not satisfied. Optionally, enforces the
//! upper bound instead of checking it.
/*!
  The bounds are \f$0 \le P \le P_o\f$ where \f$P_o\f$ is the overburden pressure.
*/
void Distributed::check_P_bounds(IceModelVec2S &P,
                                 const IceModelVec2S &P_o,
                                 bool enforce_upper) {

  IceModelVec::AccessList list{&P, &P_o};

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      if (P(i,j) < 0.0) {
        throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                      "negative subglacial water pressure\n"
                                      "P(%d, %d) = %.6f Pa",
                                      i, j, P(i, j));
      }

      if (enforce_upper) {
        P(i,j) = std::min(P(i,j), P_o(i,j));
      } else if (P(i,j) > P_o(i,j) + 0.001) {
        throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                      "subglacial water pressure P = %.16f Pa exceeds\n"
                                      "overburden pressure Po = %.16f Pa at (i,j)=(%d,%d)",
                                      P(i, j), P_o(i, j), i, j);
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

}


//! Compute functional relationship P(W) which applies only in steady state.
/*!
  In steady state in this model, water pressure is determined by a balance of
  cavitation (opening) caused by sliding and creep closure.

  This will be used in initialization when P is otherwise unknown, and
  in verification and/or reporting.  It is not used during time-dependent
  model runs.  To be more complete, \f$P = P(W,P_o,|v_b|)\f$.
*/
void Distributed::P_from_W_steady(const IceModelVec2S &W,
                                  const IceModelVec2S &P_overburden,
                                  const IceModelVec2S &sliding_speed,
                                  IceModelVec2S &result) {

  const double
    ice_softness                   = m_config->get_double("flow_law.isothermal_Glen.ice_softness"),
    creep_closure_coefficient      = m_config->get_double("hydrology.creep_closure_coefficient"),
    cavitation_opening_coefficient = m_config->get_double("hydrology.cavitation_opening_coefficient"),
    Glen_exponent                  = m_config->get_double("stress_balance.sia.Glen_exponent"),
    Wr                             = m_config->get_double("hydrology.roughness_scale");

  const double CC = cavitation_opening_coefficient / (creep_closure_coefficient * ice_softness);

  IceModelVec::AccessList list{&W, &P_overburden, &sliding_speed, &result};

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    double sb = pow(CC * sliding_speed(i, j), 1.0 / Glen_exponent);
    if (W(i, j) == 0.0) {
      // see P(W) formula in steady state; note P(W) is continuous (in steady
      // state); these facts imply:
      if (sb > 0.0) {
        // no water + cavitation = underpressure
        result(i, j) = 0.0;
      } else {
        // no water + no cavitation = creep repressurizes = overburden
        result(i, j) = P_overburden(i, j);
      }
    } else {
      double Wratio = std::max(0.0, Wr - W(i, j)) / W(i, j);
      // in cases where steady state is actually possible this will come out positive, but
      // otherwise we should get underpressure P=0, and that is what it yields
      result(i, j) = std::max(0.0, P_overburden(i, j) - sb * pow(Wratio, 1.0 / Glen_exponent));
    }
  } // end of the loop over grid points
}

double Distributed::max_timestep_P_diff(double phi0, double dt_diff_w) const {
  return 2.0 * phi0 * dt_diff_w;
}

void Distributed::update_P(double dt,
                           const IceModelVec2CellType &cell_type,
                           const IceModelVec2S &sliding_speed,
                           const IceModelVec2S &total_input,
                           const IceModelVec2S &P_overburden,
                           const IceModelVec2S &Wtill,
                           const IceModelVec2S &Wtill_new,
                           const IceModelVec2S &P,
                           const IceModelVec2S &W,
                           const IceModelVec2Stag &Ws,
                           const IceModelVec2Stag &K,
                           const IceModelVec2Stag &Q,
                           IceModelVec2S &P_new) const {

  const double
    n    = m_config->get_double("stress_balance.sia.Glen_exponent"),
    A    = m_config->get_double("flow_law.isothermal_Glen.ice_softness"),
    c1   = m_config->get_double("hydrology.cavitation_opening_coefficient"),
    c2   = m_config->get_double("hydrology.creep_closure_coefficient"),
    Wr   = m_config->get_double("hydrology.roughness_scale"),
    phi0 = m_config->get_double("hydrology.regularizing_porosity");

  // update Pnew from time step
  const double
    CC  = (m_rg * dt) / phi0,
    wux = 1.0 / (m_dx * m_dx),
    wuy = 1.0 / (m_dy * m_dy);

  IceModelVec::AccessList list{&P, &W, &Wtill, &Wtill_new, &sliding_speed, &Ws,
      &K, &Q, &total_input, &cell_type, &P_overburden, &P_new};

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    auto w = W.star(i, j);
    double P_o = P_overburden(i, j);

    if (cell_type.ice_free_land(i, j)) {
      P_new(i, j) = 0.0;
    } else if (cell_type.ocean(i, j)) {
      P_new(i, j) = P_o;
    } else if (w.ij <= 0.0) {
      P_new(i, j) = P_o;
    } else {
      auto q = Q.star(i, j);
      auto k = K.star(i, j);
      auto ws = Ws.star(i, j);

      double
        Open  = c1 * sliding_speed(i, j) * std::max(0.0, Wr - w.ij),
        Close = c2 * A * pow(P_o - P(i, j), n) * w.ij;

      // compute the flux divergence the same way as in update_W()
      const double divadflux = (q.e - q.w) / m_dx + (q.n - q.s) / m_dy;
      const double
        De = m_rg * k.e * ws.e,
        Dw = m_rg * k.w * ws.w,
        Dn = m_rg * k.n * ws.n,
        Ds = m_rg * k.s * ws.s;

      double diffW = (wux * (De * (w.e - w.ij) - Dw * (w.ij - w.w)) +
                      wuy * (Dn * (w.n - w.ij) - Ds * (w.ij - w.s)));

      double divflux = -divadflux + diffW;

      // pressure update equation
      double Wtill_change = Wtill_new(i, j) - Wtill(i, j);
      double ZZ = Close - Open + total_input(i, j) - Wtill_change / dt;

      P_new(i, j) = P(i, j) + CC * (divflux + ZZ);

      // projection to enforce  0 <= P <= P_o
      P_new(i, j) = clip(P_new(i, j), 0.0, P_o);
    }
  } // end of the loop over grid points
}


//! Update the model state variables W,P by running the subglacial hydrology model.
/*!
  Runs the hydrology model from time t to time t + dt.  Here [t,dt]
  is generally on the order of months to years.  This hydrology model will take its
  own shorter time steps, perhaps hours to weeks.
*/
void Distributed::update_impl(double t, double dt, const Inputs& inputs) {

  double
    ht  = t,
    hdt = 0.0;

  const double
    t_final = t + dt,
    dt_max  = m_config->get_double("hydrology.maximum_time_step", "seconds"),
    phi0    = m_config->get_double("hydrology.regularizing_porosity");

  // make sure W,P have valid ghosts before starting hydrology steps
  m_W.update_ghosts();
  m_P.update_ghosts();

#if (PISM_DEBUG==1)
  double tillwat_max = m_config->get_double("hydrology.tillwat_max");
#endif

  unsigned int step_counter = 0;
  for (; ht < t_final; ht += hdt) {
    step_counter++;

#if (PISM_DEBUG==1)
    double huge_number = 1e6;
    check_bounds(m_W, huge_number);
    check_bounds(m_Wtill, tillwat_max);
#endif

    // note that ice dynamics can change overburden pressure, so we can only check P
    // bounds if thk has not changed; if thk could have just changed, such as in the first
    // time through the current loop, we enforce them
    bool enforce_upper = (step_counter == 1);
    check_P_bounds(m_P, m_Pover, enforce_upper);

    water_thickness_staggered(m_W,
                              *inputs.cell_type,
                              m_Wstag);

    double maxKW = 0.0;
    compute_conductivity(m_Wstag,
                         subglacial_water_pressure(),
                         *inputs.bed_elevation,
                         m_K, maxKW);

    compute_velocity(m_Wstag,
                     subglacial_water_pressure(),
                     *inputs.bed_elevation,
                     m_K,
                     inputs.no_model_mask,
                     m_V);

    // to get Q, W needs valid ghosts
    advective_fluxes(m_V, m_W, m_Q);

    {
      const double
        dt_cfl    = max_timestep_W_cfl(),
        dt_diff_w = max_timestep_W_diff(maxKW),
        dt_diff_p = max_timestep_P_diff(phi0, dt_diff_w);

      hdt = std::min(t_final - ht, dt_max);
      hdt = std::min(hdt, dt_cfl);
      hdt = std::min(hdt, dt_diff_w);
      hdt = std::min(hdt, dt_diff_p);
    }

    m_log->message(3, "  hydrology step %05d, dt = %f s\n", step_counter, hdt);

    // update Wtillnew from Wtill and input_rate
    update_Wtill(hdt,
                 m_Wtill,
                 m_input_rate,
                 m_Wtillnew);
    // remove water in ice-free areas and account for changes
    enforce_bounds(*inputs.cell_type,
                   inputs.no_model_mask,
                   0.0,        // do not limit maximum thickness
                   m_Wtillnew,
                   m_grounded_margin_change,
                   m_grounding_line_change,
                   m_conservation_error_change,
                   m_no_model_mask_change);

    update_P(hdt,
             *inputs.cell_type,
             *inputs.ice_sliding_speed,
             m_input_rate,
             m_Pover,
             m_Wtill, m_Wtillnew,
             subglacial_water_pressure(),
             m_W, m_Wstag,
             m_K, m_Q,
             m_Pnew);

    // update Wnew from W, Wtill, Wtillnew, Wstag, Q, input_rate
    update_W(hdt,
             m_input_rate,
             m_W, m_Wstag,
             m_Wtill, m_Wtillnew,
             m_K, m_Q,
             m_Wnew);
    // remove water in ice-free areas and account for changes
    enforce_bounds(*inputs.cell_type,
                   inputs.no_model_mask,
                   0.0, // do  not limit maximum thickness
                   m_Wnew,
                   m_grounded_margin_change,
                   m_grounding_line_change,
                   m_conservation_error_change,
                   m_no_model_mask_change);

    // transfer new into old
    m_W.copy_from(m_Wnew);
    m_Wtill.copy_from(m_Wtillnew);
    m_P.copy_from(m_Pnew);
  } // end of the time-stepping loop

  m_log->message(2,
                 "  took %d hydrology sub-steps with average dt = %.6f years (%.6f s)\n",
                 step_counter,
                 units::convert(m_sys, dt/step_counter, "seconds", "years"),
                 dt/step_counter);
}

} // end of namespace hydrology
} // end of namespace pism
