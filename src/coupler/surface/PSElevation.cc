// Copyright (C) 2011, 2012, 2013, 2014, 2015 Andy Aschwanden and Constantine Khroulev
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

#include "PSElevation.hh"
#include "PIO.hh"
#include "PISMVars.hh"
#include "IceGrid.hh"
#include "PISMConfig.hh"
#include "error_handling.hh"
#include "pism_options.hh"

namespace pism {

///// Elevation-dependent temperature and surface mass balance.
PSElevation::PSElevation(const IceGrid &g)
  : SurfaceModel(g),
    climatic_mass_balance(g.config.get_unit_system(), "climatic_mass_balance", m_grid),
    ice_surface_temp(g.config.get_unit_system(), "ice_surface_temp", m_grid)
{
  // empty
}

void PSElevation::init() {
  bool m_limits_set = false;

  m_t = m_dt = GSL_NAN;  // every re-init restarts the clock

  verbPrintf(2, m_grid.com,
             "* Initializing the constant-in-time surface processes model PSElevation. Setting...\n");

  // options
  {
    // ice surface temperature
    {
      // set defaults:
      T_min   = m_grid.convert(-5.0, "Celsius", "Kelvin");
      T_max   = m_grid.convert(0.0, "Celsius", "Kelvin");
      z_T_min = 1325.0;
      z_T_max = 1350.0;

      options::RealList IST("-ice_surface_temp", "ice surface temperature parameterization");
      if (IST.is_set()) {
        if (IST->size() != 4) {
          throw RuntimeError("option -ice_surface_temp requires an argument"
                             " (comma-separated list of 4 numbers)");
        }
        T_min   = m_grid.convert(IST[0], "Celsius", "Kelvin");
        T_max   = m_grid.convert(IST[1], "Celsius", "Kelvin");
        z_T_min = IST[2];
        z_T_max = IST[3];
      }
    }

    // climatic mass balance
    {
      // set defaults:
      m_min   = m_grid.convert(-3.0, "m year-1", "m s-1");
      m_max   = m_grid.convert(4.0, "m year-1", "m s-1");
      z_m_min = 1100.0;
      z_ELA   = 1450.0;
      z_m_max = 1700.0;

      options::RealList CMB("-climatic_mass_balance",
                            "climatic mass balance parameterization");
      if (CMB.is_set()) {
        if (CMB->size() != 5) {
          throw RuntimeError("-climatic_mass_balance requires an argument"
                             " (comma-separated list of 5 numbers)");
        }
        m_min   = m_grid.convert(CMB[0], "m year-1", "m s-1");
        m_max   = m_grid.convert(CMB[1], "m year-1", "m s-1");
        z_m_min = CMB[2];
        z_ELA   = CMB[3];
        z_m_max = CMB[4];
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
        m_limit_min = m_grid.convert(m_limits[0], "m year-1", "m s-1");
        m_limit_max = m_grid.convert(m_limits[1], "m year-1", "m s-1");
      } else {
        m_limit_min = m_min;
        m_limit_max = m_max;
      }
    }
  }

  verbPrintf(3, m_grid.com,
             "     temperature at %.0f m a.s.l. = %.2f deg C\n"
             "     temperature at %.0f m a.s.l. = %.2f deg C\n"
             "     mass balance below %.0f m a.s.l. = %.2f m/year\n"
             "     mass balance at  %.0f m a.s.l. = %.2f m/year\n"
             "     mass balance at  %.0f m a.s.l. = %.2f m/year\n"
             "     mass balance above %.0f m a.s.l. = %.2f m/year\n"
             "     equilibrium line altitude z_ELA = %.2f m a.s.l.\n",
             z_T_min, T_min, z_T_max, T_max, z_m_min,
             m_grid.convert(m_limit_min, "m s-1", "m year-1"),
             z_m_min, m_min, z_m_max,
             m_grid.convert(m_max, "m s-1", "m year-1"),
             z_m_max,
             m_grid.convert(m_limit_max, "m s-1", "m year-1"), z_ELA);

  // get access to ice upper surface elevation
  usurf = m_grid.variables().get_2d_scalar("surface_altitude");

  // NCSpatialVariables storing temperature and surface mass balance metadata

  climatic_mass_balance.set_string("pism_intent", "diagnostic");
  climatic_mass_balance.set_string("long_name",
                                   "surface mass balance (accumulation/ablation) rate");
  climatic_mass_balance.set_string("standard_name",
                                   "land_ice_surface_specific_mass_balance_flux");
  climatic_mass_balance.set_units("kg m-2 s-1");
  climatic_mass_balance.set_glaciological_units("kg m-2 year-1");

  ice_surface_temp.set_string("pism_intent", "diagnostic");
  ice_surface_temp.set_string("long_name",
                              "ice temperature at the ice surface");
  ice_surface_temp.set_units("K");

  // parameterizing the ice surface temperature 'ice_surface_temp'
  verbPrintf(2, m_grid.com, "    - parameterizing the ice surface temperature 'ice_surface_temp' ... \n");
  verbPrintf(2, m_grid.com,
             "      ice temperature at the ice surface (T = ice_surface_temp) is piecewise-linear function\n"
             "        of surface altitude (usurf):\n"
             "                 /  %2.2f K                            for            usurf < %.f m\n"
             "            T = |   %5.2f K + %5.3f * (usurf - %.f m) for   %.f m < usurf < %.f m\n"
             "                 \\  %5.2f K                            for   %.f m < usurf\n",
             T_min, z_T_min,
             T_min, (T_max-T_min)/(z_T_max-z_T_min), z_T_min, z_T_min, z_T_max,
             T_max, z_T_max);

  // parameterizing the ice surface mass balance 'climatic_mass_balance'
  verbPrintf(2, m_grid.com, "    - parameterizing the ice surface mass balance 'climatic_mass_balance' ... \n");

  if (m_limits_set) {
    verbPrintf(2, m_grid.com, "    - option '-climatic_mass_balance_limits' seen, limiting upper and lower bounds ... \n");
  }

  verbPrintf(2, m_grid.com,
             "      surface mass balance (M = climatic_mass_balance) is piecewise-linear function\n"
             "        of surface altitue (usurf):\n"
             "                  /  %5.2f m/year                       for          usurf < %3.f m\n"
             "             M = |    %5.3f 1/a * (usurf-%.0f m)     for %3.f m < usurf < %3.f m\n"
             "                  \\   %5.3f 1/a * (usurf-%.0f m)     for %3.f m < usurf < %3.f m\n"
             "                   \\ %5.2f m/year                       for %3.f m < usurf\n",
             m_grid.convert(m_limit_min, "m s-1", "m year-1"), z_m_min,
             m_grid.convert(-m_min, "m s-1", "m year-1")/(z_ELA - z_m_min), z_ELA, z_m_min, z_ELA,
             m_grid.convert(m_max, "m s-1", "m year-1")/(z_m_max - z_ELA), z_ELA, z_ELA, z_m_max,
             m_grid.convert(m_limit_max, "m s-1", "m year-1"), z_m_max);
}

void PSElevation::attach_atmosphere_model(AtmosphereModel *input)
{
  delete input;
}

void PSElevation::get_diagnostics(std::map<std::string, Diagnostic*> &/*dict*/,
                                  std::map<std::string, TSDiagnostic*> &/*ts_dict*/) {
  // empty
}

void PSElevation::update(double my_t, double my_dt)
{
  m_t = my_t;
  m_dt = my_dt;
}


void PSElevation::ice_surface_mass_flux(IceModelVec2S &result) {
  double dabdz = -m_min/(z_ELA - z_m_min);
  double dacdz = m_max/(z_m_max - z_ELA);

  IceModelVec::AccessList list;
  list.add(result);
  list.add(*usurf);
  for (Points p(m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    double z = (*usurf)(i, j);
    if (z < z_m_min) {
      result(i, j) = m_limit_min;
    }
    else if ((z >= z_m_min) && (z < z_ELA)) {
      result(i, j) = dabdz * (z - z_ELA);
    }
    else if ((z >= z_ELA) && (z <= z_m_max)) {
      result(i, j) = dacdz * (z - z_ELA);
    }
    else if (z > z_m_max) {
      result(i, j) = m_limit_max;
    }
    else {
      throw RuntimeError("PSElevation::ice_surface_mass_flux: HOW DID I GET HERE?");
    }
  }

  // convert from m/s ice equivalent to kg m-2 s-1:
  result.scale(m_config.get("ice_density"));
}

void PSElevation::ice_surface_temperature(IceModelVec2S &result) {
  IceModelVec::AccessList list;
  list.add(result);
  list.add(*usurf);
  double dTdz = (T_max - T_min)/(z_T_max - z_T_min);
  for (Points p(m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    double z = (*usurf)(i, j);
    if (z <= z_T_min) {
      result(i, j) = T_min;
    }
    else if ((z > z_T_min) && (z < z_T_max)) {
      result(i, j) = T_min + dTdz * (z - z_T_min);
    }
    else if (z >= z_T_max) {
      result(i, j) = T_max;
    }
    else {
      throw RuntimeError("PSElevation::ice_surface_temperature: HOW DID I GET HERE?");
    }
  }
}

void PSElevation::add_vars_to_output(const std::string &keyword, std::set<std::string> &result) {
  if (keyword == "medium" || keyword == "big") {
    result.insert("ice_surface_temp");
    result.insert("climatic_mass_balance");
  }
}

void PSElevation::define_variables(const std::set<std::string> &vars, const PIO &nc, IO_Type nctype) {
  SurfaceModel::define_variables(vars, nc, nctype);

  if (set_contains(vars, "ice_surface_temp")) {
    ice_surface_temp.define(nc, nctype, true);
  }

  if (set_contains(vars, "climatic_mass_balance")) {
    climatic_mass_balance.define(nc, nctype, true);
  }
}

void PSElevation::write_variables(const std::set<std::string> &vars, const PIO &nc) {
  if (set_contains(vars, "ice_surface_temp")) {
    IceModelVec2S tmp;
    tmp.create(m_grid, "ice_surface_temp", WITHOUT_GHOSTS);
    tmp.metadata() = ice_surface_temp;

    ice_surface_temperature(tmp);

    tmp.write(nc);
  }

  if (set_contains(vars, "climatic_mass_balance")) {
    IceModelVec2S tmp;
    tmp.create(m_grid, "climatic_mass_balance", WITHOUT_GHOSTS);
    tmp.metadata() = climatic_mass_balance;

    ice_surface_mass_flux(tmp);
    tmp.write_in_glaciological_units = true;
    tmp.write(nc);
  }
}

} // end of namespace pism
