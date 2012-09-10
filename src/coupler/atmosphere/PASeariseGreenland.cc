// Copyright (C) 2008-2012 Ed Bueler, Constantine Khroulev, Ricarda Winkelmann,
// Gudfinna Adalgeirsdottir and Andy Aschwanden
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

// Implementation of the atmosphere model using constant-in-time precipitation
// and a cosine yearly cycle for near-surface air temperatures.

// This includes the SeaRISE Greenland parameterization.

#include "PASeariseGreenland.hh"
#include "PISMVars.hh"
#include "IceGrid.hh"
#include "pism_options.hh"
#include "PISMTime.hh"

///// PA_SeaRISE_Greenland

PetscErrorCode PA_SeaRISE_Greenland::init(PISMVars &vars) {
  PetscErrorCode ierr;

  ierr = verbPrintf(2, grid.com,
		    "* Initializing SeaRISE-Greenland atmosphere model based on the Fausto et al (2009)\n"
		    "  air temperature parameterization and using stored time-independent precipitation...\n");
  CHKERRQ(ierr);

  reference =
    "R. S. Fausto, A. P. Ahlstrom, D. V. As, C. E. Boggild, and S. J. Johnsen, 2009. "
    "A new present-day temperature parameterization for Greenland. J. Glaciol. 55 (189), 95-105.";

  ierr = PAYearlyCycle::init(vars); CHKERRQ(ierr);

  // initialize pointers to fields the parameterization depends on:
  surfelev = dynamic_cast<IceModelVec2S*>(vars.get("surface_altitude"));
  if (!surfelev) SETERRQ(grid.com, 1, "ERROR: surface_altitude is not available");

  lat = dynamic_cast<IceModelVec2S*>(vars.get("latitude"));
  if (!lat) SETERRQ(grid.com, 1, "ERROR: latitude is not available");

  lon = dynamic_cast<IceModelVec2S*>(vars.get("longitude"));
  if (!lon) SETERRQ(grid.com, 1, "ERROR: longitude is not available");

  ierr = PISMOptionsIsSet("-paleo_precip", paleo_precipitation_correction); CHKERRQ(ierr);

  if (paleo_precipitation_correction) {
    bool delta_T_set;
    string delta_T_file;

    ierr = PISMOptionsString("-paleo_precip",
                             "Specifies the air temperature offsets file to use with -paleo_precip",
                             delta_T_file, delta_T_set); CHKERRQ(ierr);

    ierr = verbPrintf(2, grid.com, 
                      "  reading delta_T data from forcing file %s for -paleo_precip actions ...\n",
                      delta_T_file.c_str());  CHKERRQ(ierr);

    delta_T = new Timeseries(grid.com, grid.rank, "delta_T",
                             grid.config.get_string("time_dimension_name"));
    ierr = delta_T->set_units("Kelvin", ""); CHKERRQ(ierr);
    ierr = delta_T->set_dimension_units(grid.time->units(), ""); CHKERRQ(ierr);
    ierr = delta_T->set_attr("long_name", "near-surface air temperature offsets");
    CHKERRQ(ierr);
    ierr = delta_T->read(delta_T_file, grid.time->use_reference_date()); CHKERRQ(ierr);
  }

  return 0;
}

PetscErrorCode PA_SeaRISE_Greenland::mean_precipitation(IceModelVec2S &result) {
  PetscErrorCode ierr;

  ierr = PAYearlyCycle::mean_precipitation(result); CHKERRQ(ierr);

  if ((delta_T != NULL) && paleo_precipitation_correction) {
    PetscReal precipexpfactor = config.get("precip_exponential_factor_for_temperature");
    ierr = result.scale(exp( precipexpfactor * (*delta_T)(t + 0.5 * dt) )); CHKERRQ(ierr);
  }

  return 0;
}

//! \brief Updates mean annual and mean July near-surface air temperatures.
//! Note that the precipitation rate is time-independent and does not need
//! to be updated.
PetscErrorCode PA_SeaRISE_Greenland::update(PetscReal my_t, PetscReal my_dt) {
  PetscErrorCode ierr;

  if (lat->has_attr("missing_at_bootstrap")) {
    ierr = PetscPrintf(grid.com, "PISM ERROR: latitude variable was missing at bootstrap;\n"
      "  SeaRISE-Greenland atmosphere model depends on latitude and would return nonsense!!\n");
      CHKERRQ(ierr);
    PISMEnd();
  }
  if (lon->has_attr("missing_at_bootstrap")) {
    ierr = PetscPrintf(grid.com, "PISM ERROR: longitude variable was missing at bootstrap;\n"
      "  SeaRISE-Greenland atmosphere model depends on longitude and would return nonsense!!\n");
      CHKERRQ(ierr);
    PISMEnd();
  }

  if ((fabs(my_t - t) < 1e-12) &&
      (fabs(my_dt - dt) < 1e-12))
    return 0;

  t  = my_t;
  dt = my_dt;

  const PetscReal 
    d_ma     = config.get("snow_temp_fausto_d_ma"),      // K
    gamma_ma = config.get("snow_temp_fausto_gamma_ma"),  // K m-1
    c_ma     = config.get("snow_temp_fausto_c_ma"),      // K (degN)-1
    kappa_ma = config.get("snow_temp_fausto_kappa_ma"),  // K (degW)-1
    d_mj     = config.get("snow_temp_fausto_d_mj"),      // SAME UNITS as for _ma ...
    gamma_mj = config.get("snow_temp_fausto_gamma_mj"),
    c_mj     = config.get("snow_temp_fausto_c_mj"),
    kappa_mj = config.get("snow_temp_fausto_kappa_mj");

  PetscScalar **lat_degN, **lon_degE, **h;
  ierr = surfelev->get_array(h);   CHKERRQ(ierr);
  ierr = lat->get_array(lat_degN); CHKERRQ(ierr);
  ierr = lon->get_array(lon_degE); CHKERRQ(ierr);
  ierr = air_temp_mean_annual.begin_access();  CHKERRQ(ierr);
  ierr = air_temp_mean_july.begin_access();  CHKERRQ(ierr);

  for (PetscInt i = grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j<grid.ys+grid.ym; ++j) {
      air_temp_mean_annual(i,j) = d_ma + gamma_ma * h[i][j] + c_ma * lat_degN[i][j] +
        kappa_ma * (-lon_degE[i][j]);
      air_temp_mean_july(i,j) = d_mj + gamma_mj * h[i][j] + c_mj * lat_degN[i][j] +
        kappa_mj * (-lon_degE[i][j]);
    }
  }

  ierr = surfelev->end_access();   CHKERRQ(ierr);
  ierr = lat->end_access(); CHKERRQ(ierr);
  ierr = lon->end_access(); CHKERRQ(ierr);
  ierr = air_temp_mean_annual.end_access();  CHKERRQ(ierr);
  ierr = air_temp_mean_july.end_access();  CHKERRQ(ierr);

  return 0;
}
