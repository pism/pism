// Copyright (C) 2012-2016 PISM Authors
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

#include "PISMHydrology.hh"
#include "base/stressbalance/PISMStressBalance.hh"
#include "base/util/Mask.hh"
#include "base/util/PISMVars.hh"
#include "base/util/error_handling.hh"
#include "base/util/io/PIO.hh"
#include "base/util/pism_options.hh"
#include "hydrology_diagnostics.hh"
#include "base/util/pism_utilities.hh"
#include "base/util/IceModelVec2CellType.hh"

namespace pism {
namespace hydrology {

Distributed::Distributed(IceGrid::ConstPtr g, stressbalance::StressBalance *sb)
  : Routing(g) {
  m_stressbalance = sb;
  m_hold_velbase_mag = false;

  // additional variables beyond hydrology::Routing::allocate()
  m_P.create(m_grid, "bwp", WITH_GHOSTS, 1);
  m_P.set_attrs("model_state",
              "pressure of transportable water in subglacial layer",
              "Pa", "");
  m_P.metadata().set_double("valid_min", 0.0);
  m_velbase_mag.create(m_grid, "velbase_mag", WITHOUT_GHOSTS);
  m_velbase_mag.set_attrs("internal",
                        "ice sliding speed seen by subglacial hydrology",
                        "m s-1", "");
  m_velbase_mag.metadata().set_double("valid_min", 0.0);
  m_Pnew.create(m_grid, "Pnew_internal", WITHOUT_GHOSTS);
  m_Pnew.set_attrs("internal",
                 "new transportable subglacial water pressure during update",
                 "Pa", "");
  m_Pnew.metadata().set_double("valid_min", 0.0);
  m_psi.create(m_grid, "hydraulic_potential", WITH_GHOSTS, 1);
  m_psi.set_attrs("internal",
                "hydraulic potential of water in subglacial layer",
                "Pa", "");
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

  m_ice_free_land_loss_cumulative      = 0.0;
  m_ocean_loss_cumulative              = 0.0;
  m_negative_thickness_gain_cumulative = 0.0;
  m_null_strip_loss_cumulative         = 0.0;

  if (init_P_from_steady) { // if so, just overwrite -i or -bootstrap value of P=bwp
    m_log->message(2,
               "  option -init_P_from_steady seen ...\n"
               "  initializing P from P(W) formula which applies in steady state\n");
    P_from_W_steady(m_P);
  }

  if (hydrology_velbase_mag_file.is_set()) {
    m_log->message(2,
               "  reading velbase_mag for 'distributed' hydrology from '%s'.\n",
               hydrology_velbase_mag_file->c_str());
    m_velbase_mag.regrid(hydrology_velbase_mag_file, CRITICAL_FILL_MISSING, 0.0);
    m_hold_velbase_mag = true;
  }
}


void Distributed::init_bwp() {

  // initialize water layer thickness from the context if present, otherwise from -i otherwise with
  // constant value

  InputOptions opts = process_input_options(m_grid->com);

  // initialize P: present or -i file or -bootstrap file or set to constant;
  //   then overwrite by regrid; then overwrite by -init_P_from_steady
  const double bwp_default = m_config->get_double("bootstrapping.defaults.bwp");

  switch (opts.type) {
  case INIT_RESTART:
  case INIT_BOOTSTRAP:
    // regridding is equivalent to reading in if grids match, but this way we can start from a file
    // that does not have 'bwp', setting it to bwp_default
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
    {"bwprel",           Diagnostic::Ptr(new Hydrology_bwprel(this))},
    {"effbwp",           Diagnostic::Ptr(new Hydrology_effbwp(this))},
    {"hydrobmelt",       Diagnostic::Ptr(new Hydrology_hydrobmelt(this))},
    {"hydroinput",       Diagnostic::Ptr(new Hydrology_hydroinput(this))},
    {"wallmelt",         Diagnostic::Ptr(new Hydrology_wallmelt(this))},
    {"bwatvel",          Diagnostic::Ptr(new Routing_bwatvel(this))},
    {"hydrovelbase_mag", Diagnostic::Ptr(new Distributed_hydrovelbase_mag(this))}
  };
  return result;
}

std::map<std::string, TSDiagnostic::Ptr> Distributed::ts_diagnostics_impl() const {
  std::map<std::string, TSDiagnostic::Ptr> result = {
    // add mass-conservation time-series diagnostics
    {"hydro_ice_free_land_loss_cumulative",      TSDiagnostic::Ptr(new MCHydrology_ice_free_land_loss_cumulative(this))},
    {"hydro_ice_free_land_loss",                 TSDiagnostic::Ptr(new MCHydrology_ice_free_land_loss(this))},
    {"hydro_ocean_loss_cumulative",              TSDiagnostic::Ptr(new MCHydrology_ocean_loss_cumulative(this))},
    {"hydro_ocean_loss",                         TSDiagnostic::Ptr(new MCHydrology_ocean_loss(this))},
    {"hydro_negative_thickness_gain_cumulative", TSDiagnostic::Ptr(new MCHydrology_negative_thickness_gain_cumulative(this))},
    {"hydro_negative_thickness_gain",            TSDiagnostic::Ptr(new MCHydrology_negative_thickness_gain(this))},
    {"hydro_null_strip_loss_cumulative",         TSDiagnostic::Ptr(new MCHydrology_null_strip_loss_cumulative(this))},
    {"hydro_null_strip_loss",                    TSDiagnostic::Ptr(new MCHydrology_null_strip_loss(this))}
  };
  return result;
}



//! Copies the P state variable which is the modeled water pressure.
void Distributed::subglacial_water_pressure(IceModelVec2S &result) const {
  result.copy_from(m_P);
}


//! Check bounds on P and fail with message if not satisfied.  Optionally, enforces the upper bound instead of checking it.
/*!
The bounds are \f$0 \le P \le P_o\f$ where \f$P_o\f$ is the overburden pressure.
 */
void Distributed::check_P_bounds(bool enforce_upper) {

  overburden_pressure(m_Pover);

  IceModelVec::AccessList list;
  list.add(m_P);
  list.add(m_Pover);

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
model runs.  To be more complete, \f$P=P(W,P_o,|v_b|)\f$.
 */
void Distributed::P_from_W_steady(IceModelVec2S &result) {
  double CC = m_config->get_double("hydrology.cavitation_opening_coefficient") /
                    (m_config->get_double("hydrology.creep_closure_coefficient") * m_config->get_double("flow_law.isothermal_Glen.ice_softness")),
    powglen = 1.0 / m_config->get_double("stress_balance.sia.Glen_exponent"), // choice is SIA; see #285
    Wr = m_config->get_double("hydrology.roughness_scale");

  overburden_pressure(m_Pover);

  IceModelVec::AccessList list;
  list.add(m_W);
  list.add(m_Pover);
  list.add(m_velbase_mag);
  list.add(result);

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    double sb = pow(CC * m_velbase_mag(i, j), powglen);
    if (m_W(i, j) == 0.0) {
      // see P(W) formula in steady state; note P(W) is continuous (in steady
      // state); these facts imply:
      if (sb > 0.0) {
        result(i, j) = 0.0;        // no water + cavitation = underpressure
      } else {
        result(i, j) = m_Pover(i, j); // no water + no cavitation = creep repressurizes = overburden
      }
    } else {
      double Wratio = std::max(0.0, Wr - m_W(i, j)) / m_W(i, j);
      // in cases where steady state is actually possible this will
      //   come out positive, but otherwise we should get underpressure P=0,
      //   and that is what it yields
      result(i, j) = std::max(0.0, m_Pover(i, j) - sb * pow(Wratio, powglen));
    }
  }
}


//! Update the the sliding speed |v_b| from ice quantities.
/*!
Calls a StressBalance method to get the vector basal velocity of the ice,
and then computes the magnitude of that.
 */
void Distributed::update_velbase_mag(IceModelVec2S &result) {
  // velbase_mag = |v_b|
  result.set_to_magnitude(m_stressbalance->advective_velocity());
}


//! Computes the adaptive time step for this (W,P) state space model.
void Distributed::adaptive_for_WandP_evolution(double t_current, double t_end, double maxKW,
                                                        double &dt_result,
                                                        double &maxV_result, double &maxD_result,
                                                        double &PtoCFLratio) {
  double dtCFL, dtDIFFW, dtDIFFP;

  adaptive_for_W_evolution(t_current,t_end, maxKW,
                           dt_result,maxV_result,maxD_result,dtCFL,dtDIFFW);

  const double phi0 = m_config->get_double("hydrology.regularizing_porosity");
  dtDIFFP = 2.0 * phi0 * dtDIFFW;

  // dt = min([te-t dtmax dtCFL dtDIFFW dtDIFFP]);
  dt_result = std::min(dt_result, dtDIFFP);

  if (dtDIFFP > 0.0) {
    PtoCFLratio = std::max(1.0, dtCFL / dtDIFFP);
  } else {
    PtoCFLratio = 1.0;
  }

  using units::convert;
  m_log->message(4,
                 "   [%.5e  %.7f  %.6f  %.9f  -->  dt = %.9f (a)  at  t = %.6f (a)]\n",
                 convert(m_sys, maxV_result, "m second-1", "m year-1"),
                 convert(m_sys, dtCFL,       "seconds",  "years"),
                 convert(m_sys, dtDIFFW,     "seconds",  "years"),
                 convert(m_sys, dtDIFFP,     "seconds",  "years"),
                 convert(m_sys, dt_result,   "seconds",  "years"),
                 convert(m_sys, t_current,   "seconds",  "years"));
}


//! Update the model state variables W,P by running the subglacial hydrology model.
/*!
Runs the hydrology model from time icet to time icet + icedt.  Here [icet,icedt]
is generally on the order of months to years.  This hydrology model will take its
own shorter time steps, perhaps hours to weeks.
 */
void Distributed::update_impl(double icet, double icedt) {

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

  // from current ice geometry/velocity variables, initialize Po and velbase_mag
  if (!m_hold_velbase_mag) {
    update_velbase_mag(m_velbase_mag);
  }

  const double
            rg    = m_config->get_double("constants.fresh_water.density") * m_config->get_double("constants.standard_gravity"),
            nglen = m_config->get_double("stress_balance.sia.Glen_exponent"), // choice is SIA; see #285
            Aglen = m_config->get_double("flow_law.isothermal_Glen.ice_softness"),
            c1    = m_config->get_double("hydrology.cavitation_opening_coefficient"),
            c2    = m_config->get_double("hydrology.creep_closure_coefficient"),
            Wr    = m_config->get_double("hydrology.roughness_scale"),
            phi0  = m_config->get_double("hydrology.regularizing_porosity");

  double ht = m_t, hdt = 0, // hydrology model time and time step
            maxKW = 0, maxV = 0, maxD = 0;
  double icefreelost = 0, oceanlost = 0, negativegain = 0, nullstriplost = 0,
            delta_icefree = 0, delta_ocean = 0, delta_neggain = 0, delta_nullstrip = 0;

  double PtoCFLratio = 0,  // for reporting ratio of dtCFL to dtDIFFP
            cumratio = 0.0;
  int hydrocount = 0; // count hydrology time steps

  while (ht < m_t + m_dt) {
    hydrocount++;

#if (PISM_DEBUG==1)
    check_water_thickness_nonnegative(m_W);
    check_Wtil_bounds();
#endif

    // note that ice dynamics can change overburden pressure, so we can only check P
    //   bounds if thk has not changed; if thk could have just changed, such as in the
    //   first time through the current loop, we enforce them
    check_P_bounds((hydrocount == 1));

    water_thickness_staggered(m_Wstag);
    m_Wstag.update_ghosts();

    conductivity_staggered(m_K,maxKW);
    m_K.update_ghosts();

    velocity_staggered(m_V);

    // to get Qstag, W needs valid ghosts
    advective_fluxes(m_Q);
    m_Q.update_ghosts();

    adaptive_for_WandP_evolution(ht, m_t+m_dt, maxKW, hdt, maxV, maxD, PtoCFLratio);
    cumratio += PtoCFLratio;

    if ((m_inputtobed != NULL) || (hydrocount==1)) {
      get_input_rate(ht,hdt,m_total_input);
    }

    // update Wtilnew from Wtil
    raw_update_Wtil(hdt);
    boundary_mass_changes(m_Wtilnew, delta_icefree, delta_ocean,
                          delta_neggain, delta_nullstrip);
    icefreelost  += delta_icefree;
    oceanlost    += delta_ocean;
    negativegain += delta_neggain;
    nullstriplost+= delta_nullstrip;

    // update Pnew from time step
    const double  CC = (rg * hdt) / phi0,
                     wux  = 1.0 / (m_dx * m_dx),
                     wuy  = 1.0 / (m_dy * m_dy);
    double Open, Close, divflux, ZZ, diffW;
    overburden_pressure(m_Pover);

    const IceModelVec2CellType &mask = *m_grid->variables().get_2d_cell_type("mask");

    IceModelVec::AccessList list;
    list.add(m_P);
    list.add(m_W);
    list.add(m_Wtil);
    list.add(m_Wtilnew);
    list.add(m_velbase_mag);
    list.add(m_Wstag);
    list.add(m_K);
    list.add(m_Q);
    list.add(m_total_input);
    list.add(mask);
    list.add(m_Pover);
    list.add(m_Pnew);

    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      if (mask.ice_free_land(i,j)) {
        m_Pnew(i,j) = 0.0;
      } else if (mask.ocean(i,j)) {
        m_Pnew(i,j) = m_Pover(i,j);
      } else if (m_W(i,j) <= 0.0) {
        m_Pnew(i,j) = m_Pover(i,j);
      } else {
        // opening and closure terms in pressure equation
        Open = std::max(0.0,c1 * m_velbase_mag(i,j) * (Wr - m_W(i,j)));
        Close = c2 * Aglen * pow(m_Pover(i,j) - m_P(i,j),nglen) * m_W(i,j);

        // compute the flux divergence the same way as in raw_update_W()
        const double divadflux =
          (m_Q(i,j,0) - m_Q(i-1,j  ,0)) / m_dx +
          (m_Q(i,j,1) - m_Q(i,  j-1,1)) / m_dy;
        const double
          De = rg * m_K(i,  j,0) * m_Wstag(i,  j,0),
          Dw = rg * m_K(i-1,j,0) * m_Wstag(i-1,j,0),
          Dn = rg * m_K(i,j  ,1) * m_Wstag(i,j  ,1),
          Ds = rg * m_K(i,j-1,1) * m_Wstag(i,j-1,1);
        diffW =   wux * (De * (m_W(i+1,j) - m_W(i,j)) - Dw * (m_W(i,j) - m_W(i-1,j)))
          + wuy * (Dn * (m_W(i,j+1) - m_W(i,j)) - Ds * (m_W(i,j) - m_W(i,j-1)));
        divflux = - divadflux + diffW;

        // pressure update equation
        ZZ = Close - Open + m_total_input(i,j) - (m_Wtilnew(i,j) - m_Wtil(i,j)) / hdt;
        m_Pnew(i,j) = m_P(i,j) + CC * (divflux + ZZ);
        // projection to enforce  0 <= P <= P_o
        m_Pnew(i,j) = std::min(std::max(0.0, m_Pnew(i,j)), m_Pover(i,j));
      }
    }

    // update Wnew from W, Wtil, Wtilnew, Wstag, Qstag, total_input
    raw_update_W(hdt);
    boundary_mass_changes(m_Wnew, delta_icefree, delta_ocean,
                          delta_neggain, delta_nullstrip);
    icefreelost  += delta_icefree;
    oceanlost    += delta_ocean;
    negativegain += delta_neggain;
    nullstriplost+= delta_nullstrip;

    // transfer new into old
    m_Wnew.update_ghosts(m_W);
    m_Wtil.copy_from(m_Wtilnew);
    m_Pnew.update_ghosts(m_P);

    ht += hdt;
  } // end of hydrology model time-stepping loop

  m_log->message(2,
             "  'distributed' hydrology took %d hydrology sub-steps"
             " with average dt = %.6f years\n",
             hydrocount, units::convert(m_sys, m_dt/hydrocount, "seconds", "years"));
  m_log->message(3,
             "  (hydrology info: dt = %.2f s,  av %.2f steps per CFL,  max |V| = %.2e m s-1,"
             "  max D = %.2e m^2 s-1)\n",
             m_dt/hydrocount, cumratio/hydrocount, maxV, maxD);

  m_ice_free_land_loss_cumulative      += icefreelost;
  m_ocean_loss_cumulative              += oceanlost;
  m_negative_thickness_gain_cumulative += negativegain;
  m_null_strip_loss_cumulative         += nullstriplost;
}


Distributed_hydrovelbase_mag::Distributed_hydrovelbase_mag(const Distributed *m)
  : Diag<Distributed>(m) {
  m_vars = {SpatialVariableMetadata(m_sys, "hydrovelbase_mag")};
  set_attrs("the version of velbase_mag seen by the 'distributed' hydrology model",
            "", "m s-1", "m year-1", 0);
}


IceModelVec::Ptr Distributed_hydrovelbase_mag::compute_impl() {
  IceModelVec2S::Ptr result(new IceModelVec2S);
  result->create(m_grid, "hydrovelbase_mag", WITHOUT_GHOSTS);
  result->metadata(0) = m_vars[0];
  result->write_in_glaciological_units = true;

  // the value reported diagnostically is merely the last value filled
  result->copy_from(model->m_velbase_mag);

  return result;
}


} // end of namespace hydrology
} // end of namespace pism
