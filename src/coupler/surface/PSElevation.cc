// Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016 Andy Aschwanden and Constantine Khroulev
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

#include <gsl/gsl_math.h>

#include "PSElevation.hh"

#include "base/util/iceModelVec.hh"
#include "base/util/io/PIO.hh"
#include "base/util/PISMVars.hh"
#include "base/util/IceGrid.hh"
#include "base/util/PISMConfigInterface.hh"
#include "base/util/error_handling.hh"
#include "base/util/pism_options.hh"
#include "base/util/io/io_helpers.hh"
#include "base/util/MaxTimestep.hh"
#include "base/util/pism_utilities.hh"

namespace pism {
namespace surface {

///// Elevation-dependent temperature and surface mass balance.
Elevation::Elevation(IceGrid::ConstPtr g)
  : SurfaceModel(g),
    m_climatic_mass_balance(m_sys, "climatic_mass_balance"),
    m_ice_surface_temp(m_sys, "ice_surface_temp")
{
  // empty
}

void Elevation::init_impl() {
  bool m_limits_set = false;

  using units::convert;

  m_t = m_dt = GSL_NAN;  // every re-init restarts the clock

  m_log->message(2,
             "* Initializing the constant-in-time surface processes model Elevation. Setting...\n");

  // options
  {
    // ice surface temperature
    {
      // set defaults:
      m_T_min   = convert(m_sys, -5.0, "Celsius", "Kelvin");
      m_T_max   = convert(m_sys, 0.0, "Celsius", "Kelvin");
      m_z_T_min = 1325.0;
      m_z_T_max = 1350.0;

      options::RealList IST("-ice_surface_temp", "ice surface temperature parameterization");
      if (IST.is_set()) {
        if (IST->size() != 4) {
          throw RuntimeError("option -ice_surface_temp requires an argument"
                             " (comma-separated list of 4 numbers)");
        }
        m_T_min   = convert(m_sys, IST[0], "Celsius", "Kelvin");
        m_T_max   = convert(m_sys, IST[1], "Celsius", "Kelvin");
        m_z_T_min = IST[2];
        m_z_T_max = IST[3];
      }
    }

    // climatic mass balance
    {
      // set defaults:
      m_M_min   = convert(m_sys, -3.0, "m year-1", "m s-1");
      m_M_max   = convert(m_sys, 4.0, "m year-1", "m s-1");
      m_z_M_min = 1100.0;
      m_z_ELA   = 1450.0;
      m_z_M_max = 1700.0;

      options::RealList CMB("-climatic_mass_balance",
                            "climatic mass balance parameterization");
      if (CMB.is_set()) {
        if (CMB->size() != 5) {
          throw RuntimeError("-climatic_mass_balance requires an argument"
                             " (comma-separated list of 5 numbers)");
        }
        m_M_min   = convert(m_sys, CMB[0], "m year-1", "m s-1");
        m_M_max   = convert(m_sys, CMB[1], "m year-1", "m s-1");
        m_z_M_min = CMB[2];
        m_z_ELA   = CMB[3];
        m_z_M_max = CMB[4];
      }
    }

    // limits of the climatic mass balance
    {
      options::RealList m_limits("-climatic_mass_balance_limits",
                                 "lower and upper limits of the climatic mass balance");
      m_limits_set = m_limits.is_set();
      if (m_limits.is_set()) {
        if (m_limits->size() != 2) {
          throw RuntimeError("-climatic_mass_balance_limits requires an argument"
                             " (a comma-separated list of 2 numbers)");
        }
        m_M_limit_min = convert(m_sys, m_limits[0], "m year-1", "m s-1");
        m_M_limit_max = convert(m_sys, m_limits[1], "m year-1", "m s-1");
      } else {
        m_M_limit_min = m_M_min;
        m_M_limit_max = m_M_max;
      }
    }
  }

  m_log->message(3,
             "     temperature at %.0f m a.s.l. = %.2f deg C\n"
             "     temperature at %.0f m a.s.l. = %.2f deg C\n"
             "     mass balance below %.0f m a.s.l. = %.2f m year-1\n"
             "     mass balance at  %.0f m a.s.l. = %.2f m year-1\n"
             "     mass balance at  %.0f m a.s.l. = %.2f m year-1\n"
             "     mass balance above %.0f m a.s.l. = %.2f m year-1\n"
             "     equilibrium line altitude z_ELA = %.2f m a.s.l.\n",
             m_z_T_min, m_T_min, m_z_T_max, m_T_max, m_z_M_min,
             convert(m_sys, m_M_limit_min, "m s-1", "m year-1"),
             m_z_M_min, m_M_min, m_z_M_max,
             convert(m_sys, m_M_max, "m s-1", "m year-1"),
             m_z_M_max,
             convert(m_sys, m_M_limit_max, "m s-1", "m year-1"), m_z_ELA);

  // SpatialVariableMetadatas storing temperature and surface mass balance metadata

  m_climatic_mass_balance.set_string("pism_intent", "diagnostic");
  m_climatic_mass_balance.set_string("long_name",
                                   "surface mass balance (accumulation/ablation) rate");
  m_climatic_mass_balance.set_string("standard_name",
                                   "land_ice_surface_specific_mass_balance_flux");
  m_climatic_mass_balance.set_string("units", "kg m-2 s-1");
  m_climatic_mass_balance.set_string("glaciological_units", "kg m-2 year-1");

  m_ice_surface_temp.set_string("pism_intent", "diagnostic");
  m_ice_surface_temp.set_string("long_name",
                              "ice temperature at the ice surface");
  m_ice_surface_temp.set_string("units", "K");

  // parameterizing the ice surface temperature 'ice_surface_temp'
  m_log->message(2, "    - parameterizing the ice surface temperature 'ice_surface_temp' ... \n");
  m_log->message(2,
             "      ice temperature at the ice surface (T = ice_surface_temp) is piecewise-linear function\n"
             "        of surface altitude (usurf):\n"
             "                 /  %2.2f K                            for            usurf < %.f m\n"
             "            T = |   %5.2f K + %5.3f * (usurf - %.f m) for   %.f m < usurf < %.f m\n"
             "                 \\  %5.2f K                            for   %.f m < usurf\n",
             m_T_min, m_z_T_min,
             m_T_min, (m_T_max-m_T_min)/(m_z_T_max-m_z_T_min), m_z_T_min, m_z_T_min, m_z_T_max,
             m_T_max, m_z_T_max);

  // parameterizing the ice surface mass balance 'climatic_mass_balance'
  m_log->message(2,
             "    - parameterizing the ice surface mass balance 'climatic_mass_balance' ... \n");

  if (m_limits_set) {
    m_log->message(2,
               "    - option '-climatic_mass_balance_limits' seen, limiting upper and lower bounds ... \n");
  }

  m_log->message(2,
             "      surface mass balance (M = climatic_mass_balance) is piecewise-linear function\n"
             "        of surface altitue (usurf):\n"
             "                  /  %5.2f m year-1                       for          usurf < %3.f m\n"
             "             M = |    %5.3f 1/a * (usurf-%.0f m)     for %3.f m < usurf < %3.f m\n"
             "                  \\   %5.3f 1/a * (usurf-%.0f m)     for %3.f m < usurf < %3.f m\n"
             "                   \\ %5.2f m year-1                       for %3.f m < usurf\n",
             convert(m_sys, m_M_limit_min, "m s-1", "m year-1"), m_z_M_min,
             convert(m_sys, -m_M_min, "m s-1", "m year-1")/(m_z_ELA - m_z_M_min), m_z_ELA, m_z_M_min, m_z_ELA,
             convert(m_sys, m_M_max, "m s-1", "m year-1")/(m_z_M_max - m_z_ELA), m_z_ELA, m_z_ELA, m_z_M_max,
             convert(m_sys, m_M_limit_max, "m s-1", "m year-1"), m_z_M_max);
}

MaxTimestep Elevation::max_timestep_impl(double t) {
  (void) t;
  return MaxTimestep();
}

void Elevation::attach_atmosphere_model_impl(atmosphere::AtmosphereModel *input)
{
  delete input;
}

void Elevation::get_diagnostics_impl(std::map<std::string, Diagnostic::Ptr> &/*dict*/,
                                  std::map<std::string, TSDiagnostic::Ptr> &/*ts_dict*/) {
  // empty
}

void Elevation::update_impl(double my_t, double my_dt)
{
  m_t = my_t;
  m_dt = my_dt;
}


void Elevation::ice_surface_mass_flux_impl(IceModelVec2S &result) {
  double dabdz = -m_M_min/(m_z_ELA - m_z_M_min);
  double dacdz = m_M_max/(m_z_M_max - m_z_ELA);

  // get access to ice upper surface elevation
  const IceModelVec2S *usurf = m_grid->variables().get_2d_scalar("surface_altitude");

  IceModelVec::AccessList list;
  list.add(result);
  list.add(*usurf);
  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      double z = (*usurf)(i, j);
      if (z < m_z_M_min) {
        result(i, j) = m_M_limit_min;
      }
      else if ((z >= m_z_M_min) && (z < m_z_ELA)) {
        result(i, j) = dabdz * (z - m_z_ELA);
      }
      else if ((z >= m_z_ELA) && (z <= m_z_M_max)) {
        result(i, j) = dacdz * (z - m_z_ELA);
      }
      else if (z > m_z_M_max) {
        result(i, j) = m_M_limit_max;
      }
      else {
        throw RuntimeError("Elevation::ice_surface_mass_flux: HOW DID I GET HERE?");
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  // convert from m second-1 ice equivalent to kg m-2 s-1:
  result.scale(m_config->get_double("ice_density"));
}

void Elevation::ice_surface_temperature_impl(IceModelVec2S &result) {

  const IceModelVec2S *usurf = m_grid->variables().get_2d_scalar("surface_altitude");

  IceModelVec::AccessList list;
  list.add(result);
  list.add(*usurf);
  double dTdz = (m_T_max - m_T_min)/(m_z_T_max - m_z_T_min);
  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      double z = (*usurf)(i, j);
      if (z <= m_z_T_min) {
        result(i, j) = m_T_min;
      }
      else if ((z > m_z_T_min) && (z < m_z_T_max)) {
        result(i, j) = m_T_min + dTdz * (z - m_z_T_min);
      }
      else if (z >= m_z_T_max) {
        result(i, j) = m_T_max;
      }
      else {
        throw RuntimeError("Elevation::ice_surface_temperature: HOW DID I GET HERE?");
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();
}

void Elevation::add_vars_to_output_impl(const std::string &keyword, std::set<std::string> &result) {
  if (keyword == "medium" || keyword == "big" || keyword == "2dbig") {
    result.insert("ice_surface_temp");
    result.insert("climatic_mass_balance");
  }
}

void Elevation::define_variables_impl(const std::set<std::string> &vars, const PIO &nc, IO_Type nctype) {
  std::string order = m_config->get_string("output_variable_order");
  SurfaceModel::define_variables_impl(vars, nc, nctype);

  if (set_contains(vars, "ice_surface_temp")) {
    io::define_spatial_variable(m_ice_surface_temp, *m_grid, nc, nctype, order, true);
  }

  if (set_contains(vars, "climatic_mass_balance")) {
    io::define_spatial_variable(m_climatic_mass_balance, *m_grid, nc, nctype, order, true);
  }
}

void Elevation::write_variables_impl(const std::set<std::string> &vars, const PIO &nc) {
  if (set_contains(vars, "ice_surface_temp")) {
    IceModelVec2S tmp;
    tmp.create(m_grid, "ice_surface_temp", WITHOUT_GHOSTS);
    tmp.metadata() = m_ice_surface_temp;

    ice_surface_temperature(tmp);

    tmp.write(nc);
  }

  if (set_contains(vars, "climatic_mass_balance")) {
    IceModelVec2S tmp;
    tmp.create(m_grid, "climatic_mass_balance", WITHOUT_GHOSTS);
    tmp.metadata() = m_climatic_mass_balance;

    ice_surface_mass_flux(tmp);
    tmp.write_in_glaciological_units = true;
    tmp.write(nc);
  }
}

} // end of namespace surface
} // end of namespace pism
