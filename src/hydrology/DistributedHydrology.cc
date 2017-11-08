// Copyright (C) 2012-2017 PISM Authors
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

#include "Distributed.hh"
#include "pism/util/Mask.hh"
#include "pism/util/Vars.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/io/PIO.hh"
#include "pism/util/pism_options.hh"
#include "hydrology_diagnostics.hh"
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

void Distributed::init() {
  m_log->message(2,
             "* Initializing the distributed, linked-cavities subglacial hydrology model...\n");

  {
    m_stripwidth = units::convert(m_sys, m_stripwidth, "m", "km");
    options::Real hydrology_null_strip("-hydrology_null_strip",
                                       "set the width, in km, of the strip around the edge "
                                       "of the computational domain in which hydrology is inactivated",
                                       m_stripwidth);
    m_stripwidth = units::convert(m_sys, hydrology_null_strip, "km", "m");
  }

  bool init_P_from_steady = options::Bool("-init_P_from_steady",
                                          "initialize P from formula P(W) which applies in steady state");

  options::String
    hydrology_velbase_mag_file("-hydrology_velbase_mag_file",
                               "Specifies a file to get velbase_mag from,"
                               " for 'distributed' hydrology model");

  Hydrology::init();

  Routing::init_bwat();

  init_bwp();

  if (init_P_from_steady) { // if so, just overwrite -i or -bootstrap value of P=bwp
    m_log->message(2,
               "  option -init_P_from_steady seen ...\n"
               "  initializing P from P(W) formula which applies in steady state\n");
    const IceModelVec2S &ice_thickness = *m_grid->variables().get_2d_scalar("land_ice_thickness");

    compute_overburden_pressure(ice_thickness, m_Pover);

    IceModelVec2S velbase_mag;  // FIXME: pass this is as an argument

    P_from_W_steady(m_W, m_Pover, velbase_mag,
                    m_P);
  }
}


void Distributed::init_bwp() {

  // initialize water layer thickness from the context if present, otherwise from -i
  // otherwise with constant value

  InputOptions opts = process_input_options(m_grid->com);

  // initialize P: present or -i file or -bootstrap file or set to constant;
  //   then overwrite by regrid; then overwrite by -init_P_from_steady
  const double bwp_default = m_config->get_double("bootstrapping.defaults.bwp");

  switch (opts.type) {
  case INIT_RESTART:
  case INIT_BOOTSTRAP:
    // regridding is equivalent to reading in if grids match, but this way we can start
    // from a file that does not have 'bwp', setting it to bwp_default
    m_P.regrid(opts.filename, OPTIONAL, bwp_default);
    break;
  case INIT_OTHER:
  default:
    m_P.set(bwp_default);
  }

  regrid("hydrology::Distributed", m_P); //  we could be asked to regrid from file
}


void Distributed::define_model_state_impl(const PIO &output) const {
  Routing::define_model_state_impl(output);
  m_P.define(output);
}

void Distributed::write_model_state_impl(const PIO &output) const {
  Routing::write_model_state_impl(output);
  m_P.write(output);
}

std::map<std::string, Diagnostic::Ptr> Distributed::diagnostics_impl() const {
  std::map<std::string, Diagnostic::Ptr> result = {
    {"effbwp",           Diagnostic::Ptr(new EffectiveBasalWaterPressure(this))},
    {"wallmelt",         Diagnostic::Ptr(new WallMelt(this))},
  };
  return combine(result, Routing::diagnostics_impl());
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

//! Check bounds on P and fail with message if not satisfied.  Optionally, enforces the upper bound instead of checking it.
/*!
The bounds are \f$0 \le P \le P_o\f$ where \f$P_o\f$ is the overburden pressure.
 */
void Distributed::check_P_bounds(bool enforce_upper) {

  IceModelVec::AccessList list{&m_P, &m_Pover};

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      if (m_P(i,j) < 0.0) {
        throw RuntimeError::formatted(PISM_ERROR_LOCATION, "disallowed negative subglacial water pressure\n"
                                      "P = %.6f Pa at (i,j)=(%d,%d)",
                                      m_P(i, j), i, j);
      }

      if (enforce_upper) {
        m_P(i,j) = std::min(m_P(i,j), m_Pover(i,j));
      } else if (m_P(i,j) > m_Pover(i,j) + 0.001) {
        throw RuntimeError::formatted(PISM_ERROR_LOCATION, "subglacial water pressure P = %.16f Pa exceeds\n"
                                      "overburden pressure Po = %.16f Pa at (i,j)=(%d,%d)",
                                      m_P(i, j), m_Pover(i, j), i, j);
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
  }
}

//! Update the model state variables W,P by running the subglacial hydrology model.
/*!
Runs the hydrology model from time icet to time icet + icedt.  Here [icet,icedt]
is generally on the order of months to years.  This hydrology model will take its
own shorter time steps, perhaps hours to weeks.
 */
void Distributed::update_impl(double icet, double icedt, const Inputs& inputs) {

  // if asked for the identical time interval versus last time, then
  //   do nothing; otherwise assume that [my_t,my_t+my_dt] is the time
  //   interval on which we are solving
  if ((fabs(icet - m_t) < 1e-12) && (fabs(icedt - m_dt) < 1e-12)) {
    return;
  }
  // update Component times: t = current time, t+dt = target time
  m_t = icet;
  m_dt = icedt;

  // make sure W,P have valid ghosts before starting hydrology steps
  m_W.update_ghosts();
  m_P.update_ghosts();

  const double
    nglen = m_config->get_double("stress_balance.sia.Glen_exponent"),
    Aglen = m_config->get_double("flow_law.isothermal_Glen.ice_softness"),
    c1    = m_config->get_double("hydrology.cavitation_opening_coefficient"),
    c2    = m_config->get_double("hydrology.creep_closure_coefficient"),
    Wr    = m_config->get_double("hydrology.roughness_scale"),
    phi0  = m_config->get_double("hydrology.regularizing_porosity");

  double
    ht    = m_t,
    hdt   = 0.0,                  // hydrology model time and time step
    maxKW = 0.0,
    maxV  = 0.0,
    maxD  = 0.0;

  double
    PtoCFLratio = 0.0,          // for reporting ratio of dtCFL to dtDIFFP
    cumratio    = 0.0;

  const IceModelVec2CellType &cell_type     = *inputs.cell_type;
  const IceModelVec2S        &bed_elevation = *inputs.bed_elevation;
  const IceModelVec2S        &ice_thickness = *inputs.ice_thickness;
  const IceModelVec2S        &sliding_speed = *inputs.ice_sliding_speed;

  get_input_rate(*inputs.basal_melt_rate,
                 *inputs.surface_input_rate,
                 cell_type,
                 m_total_input);

  compute_overburden_pressure(ice_thickness, m_Pover);

  const double
    t_final = m_t + m_dt,
    dt_max  = m_config->get_double("hydrology.maximum_time_step");

  unsigned int step_counter = 0;
  while (ht < m_t + m_dt) {
    step_counter++;

#if (PISM_DEBUG==1)
    check_water_thickness_nonnegative(m_W);
    check_Wtil_bounds();
#endif

    // note that ice dynamics can change overburden pressure, so we can only check P
    //   bounds if thk has not changed; if thk could have just changed, such as in the
    //   first time through the current loop, we enforce them
    check_P_bounds((step_counter == 1));

    water_thickness_staggered(m_W,
                              cell_type,
                              m_Wstag);

    compute_conductivity(m_Wstag,
                         subglacial_water_pressure(),
                         bed_elevation,
                         m_K, maxKW);

    compute_velocity(m_Wstag,
                     subglacial_water_pressure(),
                     bed_elevation,
                     m_K,
                     m_V);

    // to get Qstag, W needs valid ghosts
    advective_fluxes(m_V, m_W, m_Q);

    {
      const double
        dt_cfl    = max_timestep_cfl(),
        dt_diff_w = max_timestep_diffusivity(maxKW),
        dt_diff_p = 2.0 * phi0 * dt_diff_w;

      // dt = min([te-t dtmax dtCFL dtDIFFW dtDIFFP]);

      hdt = std::min(t_final - ht, dt_max);
      hdt = std::min(hdt, dt_cfl);
      hdt = std::min(hdt, dt_diff_w);
      hdt = std::min(hdt, dt_diff_p);

      if (dt_diff_p > 0.0) {
        PtoCFLratio = std::max(1.0, dt_cfl / dt_diff_p);
      } else {
        PtoCFLratio = 1.0;
      }
    }


    // {
    //   units::Converter years(m_sys, "seconds", "years");

    //   using units::convert;
    //   m_log->message(4,
    //                  "   [%.5e  %.7f  %.6f  %.9f  -->  dt = %.9f (a)  at  t = %.6f (a)]\n",
    //                  convert(m_sys, maxV, "m second-1", "m year-1"),
    //                  years(dtCFL), years(dtDIFFW), years(dtDIFFP), years(hdt), years(ht));
    // }

    cumratio += PtoCFLratio;

    // update Wtilnew from Wtil
    raw_update_Wtil(hdt);
    // correct water thickness and account for the changes

    // update Pnew from time step
    const double
      CC  = (m_rg * hdt) / phi0,
      wux = 1.0 / (m_dx * m_dx),
      wuy = 1.0 / (m_dy * m_dy);
    double diffW;

    IceModelVec::AccessList list{&m_P, &m_W, &m_Wtil, &m_Wtilnew, &sliding_speed, &m_Wstag,
        &m_K, &m_Q, &m_total_input, &cell_type, &m_Pover, &m_Pnew};

    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      if (cell_type.ice_free_land(i,j)) {
        m_Pnew(i,j) = 0.0;
      } else if (cell_type.ocean(i,j)) {
        m_Pnew(i,j) = m_Pover(i,j);
      } else if (m_W(i,j) <= 0.0) {
        m_Pnew(i,j) = m_Pover(i,j);
      } else {
        // opening and closure terms in pressure equation
        double Open = std::max(0.0,c1 * sliding_speed(i,j) * (Wr - m_W(i,j)));
        double Close = c2 * Aglen * pow(m_Pover(i,j) - m_P(i,j),nglen) * m_W(i,j);

        // compute the flux divergence the same way as in raw_update_W()
        const double divadflux =
          (m_Q(i,j,0) - m_Q(i-1,j  ,0)) / m_dx +
          (m_Q(i,j,1) - m_Q(i,  j-1,1)) / m_dy;
        const double
          De = m_rg * m_K(i,  j,0) * m_Wstag(i,  j,0),
          Dw = m_rg * m_K(i-1,j,0) * m_Wstag(i-1,j,0),
          Dn = m_rg * m_K(i,j  ,1) * m_Wstag(i,j  ,1),
          Ds = m_rg * m_K(i,j-1,1) * m_Wstag(i,j-1,1);
        diffW =   wux * (De * (m_W(i+1,j) - m_W(i,j)) - Dw * (m_W(i,j) - m_W(i-1,j)))
          + wuy * (Dn * (m_W(i,j+1) - m_W(i,j)) - Ds * (m_W(i,j) - m_W(i,j-1)));
        double divflux = - divadflux + diffW;

        // pressure update equation
        double ZZ = Close - Open + m_total_input(i,j) - (m_Wtilnew(i,j) - m_Wtil(i,j)) / hdt;
        m_Pnew(i,j) = m_P(i,j) + CC * (divflux + ZZ);
        // projection to enforce  0 <= P <= P_o
        m_Pnew(i,j) = std::min(std::max(0.0, m_Pnew(i,j)), m_Pover(i,j));
      }
    }

    // update Wnew from W, Wtil, Wtilnew, Wstag, Qstag, total_input
    raw_update_W(hdt);
    // correct water thickness and account for the changes

    // transfer new into old
    m_Wnew.update_ghosts(m_W);
    m_Wtil.copy_from(m_Wtilnew);
    m_Pnew.update_ghosts(m_P);

    ht += hdt;
  } // end of hydrology model time-stepping loop

  m_log->message(2,
             "  'distributed' hydrology took %d hydrology sub-steps"
             " with average dt = %.6f years\n",
             step_counter, units::convert(m_sys, m_dt/step_counter, "seconds", "years"));
  m_log->message(3,
             "  (hydrology info: dt = %.2f s,  av %.2f steps per CFL,  max |V| = %.2e m s-1,"
             "  max D = %.2e m^2 s-1)\n",
             m_dt/step_counter, cumratio/step_counter, maxV, maxD);
}

} // end of namespace hydrology
} // end of namespace pism
