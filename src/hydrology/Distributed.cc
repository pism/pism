// Copyright (C) 2012-2019, 2021, 2022, 2023 PISM Authors
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

#include <algorithm> // std::min, std::max

#include "pism/geometry/Geometry.hh"
#include "pism/hydrology/Distributed.hh"
#include "pism/util/array/CellType.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/io/File.hh"
#include "pism/util/pism_utilities.hh"

namespace pism {
namespace hydrology {

Distributed::Distributed(std::shared_ptr<const Grid> g)
    : Routing(g), m_P(m_grid, "bwp"), m_Pnew(m_grid, "Pnew_internal") {

  // additional variables beyond hydrology::Routing
  m_P.metadata(0)
      .long_name("pressure of transportable water in subglacial layer")
      .units("Pa");
  m_P.metadata()["valid_min"] = { 0.0 };

  m_Pnew.metadata(0)
      .long_name("new transportable subglacial water pressure during update")
      .units("Pa");
  m_Pnew.metadata()["valid_min"] = { 0.0 };
}

void Distributed::initialization_message() const {
  m_log->message(2,
                 "* Initializing the distributed, linked-cavities subglacial hydrology model...\n");
}

void Distributed::restart_impl(const File &input_file, int record) {
  Routing::restart_impl(input_file, record);

  m_P.read(input_file, record);

  regrid("Hydrology", m_P);
}

void Distributed::bootstrap_impl(const File &input_file, const array::Scalar &ice_thickness) {
  Routing::bootstrap_impl(input_file, ice_thickness);

  double bwp_default = m_config->get_number("bootstrapping.defaults.bwp");
  m_P.regrid(input_file, io::Default(bwp_default));

  regrid("Hydrology", m_P);

  bool init_P_from_steady = m_config->get_flag("hydrology.distributed.init_p_from_steady");

  if (init_P_from_steady) { // if so, just overwrite -i or -bootstrap value of P=bwp
    m_log->message(2, "  initializing P from P(W) formula which applies in steady state\n");

    compute_overburden_pressure(ice_thickness, m_Pover);

    array::Scalar sliding_speed(m_grid, "velbase_mag");
    sliding_speed.metadata(0).long_name("basal sliding speed").units("m s^-1");

    std::string filename = m_config->get_string("hydrology.distributed.sliding_speed_file");

    if (filename.empty()) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                    "hydrology.distributed.sliding_speed_file is not set");
    }

    sliding_speed.regrid(filename, io::Default::Nil());

    P_from_W_steady(m_W, m_Pover, sliding_speed,
                    m_P);
  }
}

void Distributed::init_impl(const array::Scalar &W_till,
                              const array::Scalar &W,
                              const array::Scalar &P) {
  Routing::init_impl(W_till, W, P);

  m_P.copy_from(P);
}

void Distributed::define_model_state_impl(const File &output) const {
  Routing::define_model_state_impl(output);
  m_P.define(output, io::PISM_DOUBLE);
}

void Distributed::write_model_state_impl(const File &output) const {
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
const array::Scalar& Distributed::subglacial_water_pressure() const {
  return m_P;
}

//! Check bounds on P and fail with message if not satisfied. Optionally, enforces the
//! upper bound instead of checking it.
/*!
  The bounds are \f$0 \le P \le P_o\f$ where \f$P_o\f$ is the overburden pressure.
*/
void Distributed::check_P_bounds(array::Scalar &P,
                                 const array::Scalar &P_o,
                                 bool enforce_upper) {

  array::AccessScope list{&P, &P_o};

  ParallelSection loop(m_grid->com);
  try {
    for (auto p = m_grid->points(); p; p.next()) {
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
void Distributed::P_from_W_steady(const array::Scalar &W,
                                  const array::Scalar &P_overburden,
                                  const array::Scalar &sliding_speed,
                                  array::Scalar &result) {

  const double
    ice_softness                   = m_config->get_number("flow_law.isothermal_Glen.ice_softness"),
    creep_closure_coefficient      = m_config->get_number("hydrology.creep_closure_coefficient"),
    cavitation_opening_coefficient = m_config->get_number("hydrology.cavitation_opening_coefficient"),
    Glen_exponent                  = m_config->get_number("stress_balance.sia.Glen_exponent"),
    Wr                             = m_config->get_number("hydrology.roughness_scale");

  const double CC = cavitation_opening_coefficient / (creep_closure_coefficient * ice_softness);

  array::AccessScope list{&W, &P_overburden, &sliding_speed, &result};

  for (auto p = m_grid->points(); p; p.next()) {
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
                           const array::CellType &cell_type,
                           const array::Scalar &sliding_speed,
                           const array::Scalar &surface_input_rate,
                           const array::Scalar &basal_melt_rate,
                           const array::Scalar &P_overburden,
                           const array::Scalar &Wtill,
                           const array::Scalar &Wtill_new,
                           const array::Scalar &P,
                           const array::Scalar1 &W,
                           const array::Staggered1 &Ws,
                           const array::Staggered1 &K,
                           const array::Staggered1 &Q,
                           array::Scalar &P_new) const {

  const double
    n    = m_config->get_number("stress_balance.sia.Glen_exponent"),
    A    = m_config->get_number("flow_law.isothermal_Glen.ice_softness"),
    c1   = m_config->get_number("hydrology.cavitation_opening_coefficient"),
    c2   = m_config->get_number("hydrology.creep_closure_coefficient"),
    Wr   = m_config->get_number("hydrology.roughness_scale"),
    phi0 = m_config->get_number("hydrology.regularizing_porosity");

  // update Pnew from time step
  const double
    CC  = (m_rg * dt) / phi0,
    wux = 1.0 / (m_dx * m_dx),
    wuy = 1.0 / (m_dy * m_dy);

  array::AccessScope list{&P, &W, &Wtill, &Wtill_new, &sliding_speed, &Ws,
                               &K, &Q, &surface_input_rate, &basal_melt_rate,
                               &cell_type, &P_overburden, &P_new};

  for (auto p = m_grid->points(); p; p.next()) {
    const int i = p.i(), j = p.j();

    auto w = W.star(i, j);
    double P_o = P_overburden(i, j);

    if (cell_type.ice_free_land(i, j)) {
      P_new(i, j) = 0.0;
    } else if (cell_type.ocean(i, j)) {
      P_new(i, j) = P_o;
    } else if (w.c <= 0.0) {
      P_new(i, j) = P_o;
    } else {
      auto q = Q.star(i, j);
      auto k = K.star(i, j);
      auto ws = Ws.star(i, j);

      double
        Open  = c1 * sliding_speed(i, j) * std::max(0.0, Wr - w.c),
        Close = c2 * A * pow(P_o - P(i, j), n) * w.c;

      // compute the flux divergence the same way as in update_W()
      const double divadflux = (q.e - q.w) / m_dx + (q.n - q.s) / m_dy;
      const double
        De = m_rg * k.e * ws.e,
        Dw = m_rg * k.w * ws.w,
        Dn = m_rg * k.n * ws.n,
        Ds = m_rg * k.s * ws.s;

      double diffW = (wux * (De * (w.e - w.c) - Dw * (w.c - w.w)) +
                      wuy * (Dn * (w.n - w.c) - Ds * (w.c - w.s)));

      double divflux = -divadflux + diffW;

      // pressure update equation
      double Wtill_change = Wtill_new(i, j) - Wtill(i, j);
      double total_input = surface_input_rate(i, j) + basal_melt_rate(i, j);
      double ZZ = Close - Open + total_input - Wtill_change / dt;

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

  ice_bottom_surface(*inputs.geometry, m_bottom_surface);

  double
    ht  = t,
    hdt = 0.0;

  const double
    t_final     = t + dt,
    dt_max      = m_config->get_number("hydrology.maximum_time_step", "seconds"),
    phi0        = m_config->get_number("hydrology.regularizing_porosity"),
    tillwat_max = m_config->get_number("hydrology.tillwat_max");

  m_Qstag_average.set(0.0);

  // make sure W,P have valid ghosts before starting hydrology steps
  m_W.update_ghosts();
  m_P.update_ghosts();

  unsigned int step_counter = 0;
  for (; ht < t_final; ht += hdt) {
    step_counter++;

#if (Pism_DEBUG==1)
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
                              inputs.geometry->cell_type,
                              m_Wstag);

    double maxKW = 0.0;
    compute_conductivity(m_Wstag,
                         subglacial_water_pressure(),
                         m_bottom_surface,
                         m_Kstag, maxKW);

    compute_velocity(m_Wstag,
                     subglacial_water_pressure(),
                     m_bottom_surface,
                     m_Kstag,
                     inputs.no_model_mask,
                     m_Vstag);

    // to get Q, W needs valid ghosts
    advective_fluxes(m_Vstag, m_W, m_Qstag);

    m_Qstag_average.add(hdt, m_Qstag);

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
                 m_surface_input_rate,
                 m_basal_melt_rate,
                 m_Wtillnew);
    // remove water in ice-free areas and account for changes
    enforce_bounds(inputs.geometry->cell_type,
                   inputs.no_model_mask,
                   0.0,           // do not limit maximum thickness
                   tillwat_max,   // till water thickness under the ocean
                   m_Wtillnew,
                   m_grounded_margin_change,
                   m_grounding_line_change,
                   m_conservation_error_change,
                   m_no_model_mask_change);

    update_P(hdt,
             inputs.geometry->cell_type,
             *inputs.ice_sliding_speed,
             m_surface_input_rate,
             m_basal_melt_rate,
             m_Pover,
             m_Wtill, m_Wtillnew,
             subglacial_water_pressure(),
             m_W, m_Wstag,
             m_Kstag, m_Qstag,
             m_Pnew);

    // update Wnew from W, Wtill, Wtillnew, Wstag, Q, input_rate
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
                   0.0, // do not limit maximum thickness
                   0.0, // transportable water layer thickness under the ocean
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

  staggered_to_regular(inputs.geometry->cell_type, m_Qstag_average,
                       m_config->get_flag("hydrology.routing.include_floating_ice"),
                       m_Q);
  m_Q.scale(1.0 / dt);

  m_log->message(2,
                 "  took %d hydrology sub-steps with average dt = %.6f years (%.6f s)\n",
                 step_counter,
                 units::convert(m_sys, dt/step_counter, "seconds", "years"),
                 dt/step_counter);
}

} // end of namespace hydrology
} // end of namespace pism
