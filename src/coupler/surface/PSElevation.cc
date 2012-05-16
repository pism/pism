// Copyright (C) 2011, 2012 Andy Aschwanden and Constantine Khroulev
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 2 of the License, or (at your option) any later
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


///// Elevation-dependent temperature and surface mass balance.

PetscErrorCode PSElevation::init(PISMVars &vars) {
  PetscErrorCode ierr;
  PetscBool T_is_set, m_is_set, m_limits_set;

  ierr = verbPrintf(2, grid.com,
                    "* Initializing the constant-in-time surface processes model PSElevation. Setting...\n"); CHKERRQ(ierr);

  ierr = PetscOptionsBegin(grid.com, "", "Elevation-dependent surface model options", ""); CHKERRQ(ierr);
  {
    PetscInt T_param_number = 4;
    PetscReal T_array[4] = {-5, 0, 1325, 1350};

    ierr = PetscOptionsGetRealArray(PETSC_NULL, "-ice_surface_temp", T_array, &T_param_number, &T_is_set);
    CHKERRQ(ierr);

    T_min = convert(T_array[0], "Celsius", "Kelvin");
    T_max = convert(T_array[1], "Celsius", "Kelvin");
    z_T_min = T_array[2];
    z_T_max = T_array[3];

    PetscInt m_param_number = 5;
    PetscReal m_array[5] = {-3, 4, 1100, 1450, 1700};

    ierr = PetscOptionsGetRealArray(PETSC_NULL, "-climatic_mass_balance", m_array, &m_param_number, &m_is_set);
    CHKERRQ(ierr);

    m_min = convert(m_array[0], "m year-1", "m s-1");
    m_max = convert(m_array[1], "m year-1", "m s-1");
    z_m_min = m_array[2];
    z_ELA = m_array[3];
    z_m_max = m_array[4];

    PetscInt Nlimitsparam= 2;
    PetscReal limitsarray[2] = {0, 0};

    ierr = PetscOptionsGetRealArray(PETSC_NULL,
                                    "-climatic_mass_balance_limits",
                                    limitsarray, &Nlimitsparam, &m_limits_set); CHKERRQ(ierr);

    if (m_limits_set)
      {
        m_limit_min = convert(limitsarray[0], "m year-1", "m s-1");
        m_limit_max = convert(limitsarray[1], "m year-1", "m s-1");
      }
    else
      {
        m_limit_min = m_min;
        m_limit_max = m_max;
      }

  }
  ierr = PetscOptionsEnd(); CHKERRQ(ierr);

  ierr = verbPrintf(3, grid.com,
                    "     temperature at %.0f m a.s.l. = %.2f deg C\n"
                    "     temperature at %.0f m a.s.l. = %.2f deg C\n"
                    "     mass balance below %.0f m a.s.l. = %.2f m/a\n"
                    "     mass balance at  %.0f m a.s.l. = %.2f m/a\n"
                    "     mass balance at  %.0f m a.s.l. = %.2f m/a\n"
                    "     mass balance above %.0f m a.s.l. = %.2f m/a\n"
                    "     equilibrium line altitude z_ELA = %.2f m a.s.l.\n",
                    z_T_min, T_min, z_T_max, T_max, z_m_min,
                    convert(m_limit_min, "m s-1", "m year-1"),
                    z_m_min, m_min, z_m_max,
                    convert(m_max, "m s-1", "m year-1"),
                    z_m_max,
                    convert(m_limit_max, "m s-1", "m year-1"), z_ELA); CHKERRQ(ierr);

  // get access to ice upper surface elevation
  usurf = dynamic_cast<IceModelVec2S*>(vars.get("surface_altitude"));
  if (!usurf) SETERRQ(grid.com, 12, "ERROR: 'usurf' is not available or is wrong type in dictionary");


  // NCSpatialVariables storing temperature and surface mass balance metadata

  climatic_mass_balance.init_2d("climatic_mass_balance", grid);
  climatic_mass_balance.set_string("pism_intent", "diagnostic");
  climatic_mass_balance.set_string("long_name",
                  "ice-equivalent surface mass balance (accumulation/ablation) rate");
  climatic_mass_balance.set_string("standard_name",
                  "land_ice_surface_specific_mass_balance");
  ierr = climatic_mass_balance.set_units("m s-1"); CHKERRQ(ierr);
  ierr = climatic_mass_balance.set_glaciological_units("m year-1"); CHKERRQ(ierr);

  ice_surface_temp.init_2d("ice_surface_temp", grid);
  ice_surface_temp.set_string("pism_intent", "diagnostic");
  ice_surface_temp.set_string("long_name",
                  "ice temperature at the ice surface");
  ierr = ice_surface_temp.set_units("K"); CHKERRQ(ierr);

  // parameterizing the ice surface temperature 'ice_surface_temp'
  ierr = verbPrintf(2, grid.com, "    - parameterizing the ice surface temperature 'ice_surface_temp' ... \n"); CHKERRQ(ierr);
  ierr = verbPrintf(2, grid.com,
                    "      ice temperature at the ice surface (T = ice_surface_temp) is piecewise-linear function\n"
                    "        of surface altitude (usurf):\n"
                    "                 /  %2.2f K                            for            usurf < %.f m\n"
                    "            T = |   %5.2f K + %5.3f * (usurf - %.f m) for   %.f m < usurf < %.f m\n"
                    "                 \\  %5.2f K                            for   %.f m < usurf\n",
                    T_min, z_T_min,
                    T_min, (T_max-T_min)/(z_T_max-z_T_min), z_T_min, z_T_min, z_T_max,
                    T_max, z_T_max); CHKERRQ(ierr);

  // parameterizing the ice surface mass balance 'climatic_mass_balance'
  ierr = verbPrintf(2, grid.com, "    - parameterizing the ice surface mass balance 'climatic_mass_balance' ... \n"); CHKERRQ(ierr);

  if (m_limits_set)
    {
      ierr = verbPrintf(2, grid.com, "    - option '-climatic_mass_balance_limits' seen, limiting upper and lower bounds ... \n"); CHKERRQ(ierr);
    }

  ierr = verbPrintf(2, grid.com,
                    "      surface mass balance (M = climatic_mass_balance) is piecewise-linear function\n"
                    "        of surface altitue (usurf):\n"
                    "                  /  %5.2f m/a                       for          usurf < %3.f m\n"
                    "             M = |    %5.3f 1/a * (usurf-%.0f m)     for %3.f m < usurf < %3.f m\n"
                    "                  \\   %5.3f 1/a * (usurf-%.0f m)     for %3.f m < usurf < %3.f m\n"
                    "                   \\ %5.2f m/a                       for %3.f m < usurf\n",
                    convert(m_limit_min, "m s-1", "m year-1"), z_m_min,
                    convert(-m_min, "m s-1", "m year-1")/(z_ELA - z_m_min), z_ELA, z_m_min, z_ELA,
                    convert(m_max, "m s-1", "m year-1")/(z_m_max - z_ELA), z_ELA, z_ELA, z_m_max,
                    convert(m_limit_max, "m s-1", "m year-1"), z_m_max); CHKERRQ(ierr);


  return 0;
}

PetscErrorCode PSElevation::ice_surface_mass_flux(IceModelVec2S &result) {
  PetscErrorCode ierr;
  PetscReal dabdz = -m_min/(z_ELA - z_m_min);
  PetscReal dacdz = m_max/(z_m_max - z_ELA);
  string history  = "elevation-dependent surface mass balance\n";

  ierr = result.begin_access(); CHKERRQ(ierr);
  ierr = usurf->begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      PetscReal z = (*usurf)(i, j);
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
        SETERRQ(grid.com, 1, "HOW DID I GET HERE? ... ending...\n");
      }
      ierr = verbPrintf(5, grid.com, "!!!!! z=%.2f, climatic_mass_balance=%.2f\n", z, result(i, j)*secpera); CHKERRQ(ierr);
    }
  }
  ierr = usurf->end_access(); CHKERRQ(ierr);
  ierr = result.end_access(); CHKERRQ(ierr);

  ierr = result.set_attr("history", history); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PSElevation::ice_surface_temperature(IceModelVec2S &result) {
  PetscErrorCode ierr;
  string history  = "elevation-dependent ice surface temperature \n";

  ierr = result.begin_access(); CHKERRQ(ierr);
  ierr = usurf->begin_access(); CHKERRQ(ierr);
  PetscReal dTdz = (T_max - T_min)/(z_T_max - z_T_min);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      PetscReal z = (*usurf)(i, j);
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
        SETERRQ(grid.com, 1, "HOW DID I GET HERE? ... ending...\n");
      }
      ierr = verbPrintf(5, grid.com,
                        "!!!!! z=%.2f, T_min=%.2f, dTdz=%.2f, dz=%.2f, T=%.2f\n",
                        z, T_min, dTdz, z - z_T_min, result(i, j)); CHKERRQ(ierr);
    }
  }
  ierr = usurf->end_access(); CHKERRQ(ierr);
  ierr = result.end_access(); CHKERRQ(ierr);

  ierr = result.set_attr("history", history); CHKERRQ(ierr);

  return 0;
}

void PSElevation::add_vars_to_output(string keyword, map<string, NCSpatialVariable> &result) {
  if (keyword == "medium" || keyword == "big") {
    result["ice_surface_temp"] = ice_surface_temp;
    result["climatic_mass_balance"] = climatic_mass_balance;
  }
}

PetscErrorCode PSElevation::define_variables(set<string> vars, const PIO &nc, PISM_IO_Type nctype) {
  PetscErrorCode ierr;

  ierr = PISMSurfaceModel::define_variables(vars, nc, nctype); CHKERRQ(ierr);

  if (set_contains(vars, "ice_surface_temp")) {
    ierr = ice_surface_temp.define(nc, nctype, true); CHKERRQ(ierr);
  }

  if (set_contains(vars, "climatic_mass_balance")) {
    ierr = climatic_mass_balance.define(nc, nctype, true); CHKERRQ(ierr);
  }

  return 0;
}

PetscErrorCode PSElevation::write_variables(set<string> vars, string filename) {
  PetscErrorCode ierr;

  if (set_contains(vars, "ice_surface_temp")) {
    IceModelVec2S tmp;
    ierr = tmp.create(grid, "ice_surface_temp", false); CHKERRQ(ierr);
    ierr = tmp.set_metadata(ice_surface_temp, 0); CHKERRQ(ierr);

    ierr = ice_surface_temperature(tmp); CHKERRQ(ierr);

    ierr = tmp.write(filename.c_str()); CHKERRQ(ierr);
  }

  if (set_contains(vars, "climatic_mass_balance")) {
    IceModelVec2S tmp;
    ierr = tmp.create(grid, "climatic_mass_balance", false); CHKERRQ(ierr);
    ierr = tmp.set_metadata(climatic_mass_balance, 0); CHKERRQ(ierr);

    ierr = ice_surface_mass_flux(tmp); CHKERRQ(ierr);
    tmp.write_in_glaciological_units = true;
    ierr = tmp.write(filename.c_str()); CHKERRQ(ierr);
  }

  return 0;
}
