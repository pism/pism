// Copyright (C) 2010, 2011 Constantine Khroulev
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

#include "PSDirectForcing.hh"

//! \brief Computes year modulo bc_period if bc_period is active.
PetscReal PSDirectForcing::my_mod(PetscReal input) {
  if (bc_period < 0.01) return input;

  // compute time since the reference year:
  PetscReal delta = input - bc_reference_year;

  // compute delta mod bc_period:
  return delta - floor(delta / bc_period) * bc_period;
}


PetscErrorCode PSDirectForcing::init(PISMVars &vars) {
  PetscErrorCode ierr;
  string filename;
  bool bc_file_set, bc_period_set, bc_ref_year_set;

  enable_time_averaging = false;

  ierr = verbPrintf(2, grid.com,
                    "* Initializing the surface model reading temperature at the top of the ice\n"
                    "  and ice surface mass flux from a file...\n"); CHKERRQ(ierr);

  ierr = PetscOptionsBegin(grid.com, "", "Direct forcing options", ""); CHKERRQ(ierr);
  {
    ierr = PISMOptionsString("-bc_file", "Specifies a file with top-surface boundary conditions",
                             filename, bc_file_set); CHKERRQ(ierr);
    ierr = PISMOptionsReal("-bc_period", "Specifies the length of the climate data period",
                           bc_period, bc_period_set); CHKERRQ(ierr);
    ierr = PISMOptionsReal("-bc_reference_year", "Boundary condition reference year",
                           bc_reference_year, bc_ref_year_set); CHKERRQ(ierr);
    ierr = PISMOptionsReal("-bc_artm_lapse_rate", "Boundary condition temperature lapse rate in Kelvin per kilometer",
                           bc_artm_lapse_rate, bc_artm_lapse_rate_set); CHKERRQ(ierr);
    ierr = PISMOptionsReal("-bc_acab_lapse_rate", "Boundary condition mass balance lapse rate in meters per year per kilometer",
                           bc_acab_lapse_rate, bc_acab_lapse_rate_set); CHKERRQ(ierr);
    ierr = PISMOptionsIsSet("-bc_time_average", "Enable time-averaging of boundary condition data",
                            enable_time_averaging); CHKERRQ(ierr);
  }
  ierr = PetscOptionsEnd(); CHKERRQ(ierr);

  if (bc_file_set == false) {
    PetscPrintf(grid.com, "PISM ERROR: option -bc_file is required.\n");
    PISMEnd();
  }

  if (bc_period_set == false) {
    bc_period = 0;
  }

  if (bc_ref_year_set == false) {
    bc_reference_year = 0;
  }

  ierr = verbPrintf(2,grid.com,
                    "    reading boundary conditions from %s ...\n",
                    filename.c_str()); CHKERRQ(ierr);

  temperature.set_n_records((unsigned int) config.get("climate_forcing_buffer_size"));
  ierr = temperature.create(grid, "artm", false); CHKERRQ(ierr);
  ierr = temperature.set_attrs("climate_forcing",
                               "temperature of the ice at the ice surface but below firn processes",
                               "Kelvin", ""); CHKERRQ(ierr);
  ierr = temperature.init(filename); CHKERRQ(ierr);

  temperature.strict_timestep_limit = ! enable_time_averaging;

  mass_flux.set_n_records((unsigned int) config.get("climate_forcing_buffer_size"));
  ierr = mass_flux.create(grid, "acab", false); CHKERRQ(ierr);
  ierr = mass_flux.set_attrs("climate_forcing",
                             "ice-equivalent surface mass balance (accumulation/ablation) rate",
                             "m s-1", "land_ice_surface_specific_mass_balance"); CHKERRQ(ierr);
  ierr = mass_flux.init(filename); CHKERRQ(ierr);

  mass_flux.strict_timestep_limit = ! enable_time_averaging;

  if (bc_artm_lapse_rate_set) {
    ierr = verbPrintf(2,grid.com,
                      "    reading reference surface from %s ...\n",
                      filename.c_str()); CHKERRQ(ierr);

    bc_surface.set_n_records((unsigned int) config.get("climate_forcing_buffer_size"));
    ierr = bc_surface.create(grid, "usurf", false); CHKERRQ(ierr);
    ierr = bc_surface.set_attrs("climate_forcing",
                               "reference surface for lapse rate",
                               "m", "surface_altitude"); CHKERRQ(ierr);
    ierr = bc_surface.init(filename); CHKERRQ(ierr);

    bc_surface.strict_timestep_limit = ! enable_time_averaging;

    surface = dynamic_cast<IceModelVec2S*>(vars.get("surface_altitude"));    
    if (!surface) SETERRQ(1, "ERROR: 'usurf' is not available or is wrong type in dictionary");

    vH = dynamic_cast<IceModelVec2S*>(vars.get("land_ice_thickness"));    
    if (!vH) SETERRQ(1, "ERROR: 'ice thickness' is not available or is wrong type in dictionary");

  }

  return 0;
}


PetscErrorCode PSDirectForcing::update(PetscReal t_years, PetscReal dt_years) {
  PetscErrorCode ierr;

  // "Periodize" the climate:
  t_years = my_mod(t_years);

  if ((fabs(t_years - t) < 1e-12) &&
      (fabs(dt_years - dt) < 1e-12))
    return 0;

  t  = t_years;
  dt = dt_years;

  ierr = temperature.update(t, dt); CHKERRQ(ierr);
  ierr = bc_surface.update(t, dt); CHKERRQ(ierr);
  ierr = mass_flux.update(t, dt); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PSDirectForcing::max_timestep(PetscReal t_years, PetscReal &dt_years) {
  PetscReal max_dt = -1;

  // "Periodize" the climate:
  t_years = my_mod(t_years);

  dt_years = temperature.max_timestep(t_years);

  max_dt = mass_flux.max_timestep(t_years);

  if (dt_years > 0) {
    if (max_dt > 0)
      dt_years = PetscMin(max_dt, dt_years);
  }
  else dt_years = max_dt;

  // If the user asked for periodized climate, limit time-steps so that PISM
  // never tries to average data over an interval that begins in one period and
  // ends in the next one.
  if (bc_period > 1e-6)
    dt_years = PetscMin(dt_years, bc_period - t_years);

  return 0;
}


void PSDirectForcing::add_vars_to_output(string keyword, set<string> &result) {
  if (keyword != "small") {
    result.insert("artm");
    result.insert("acab");
  }
}


PetscErrorCode PSDirectForcing::define_variables(set<string> vars, const NCTool &nc, nc_type nctype) {
  PetscErrorCode ierr;

  if (set_contains(vars, "artm")) {
    ierr = temperature.define(nc, nctype); CHKERRQ(ierr); 
  }

  if (set_contains(vars, "acab")) {
    ierr = mass_flux.define(nc, nctype); CHKERRQ(ierr);
  }

  return 0;
}



PetscErrorCode PSDirectForcing::write_variables(set<string> vars, string filename) {
  PetscErrorCode ierr;

  if (set_contains(vars, "artm")) {
    if (enable_time_averaging) {
      ierr = temperature.average(t, dt); CHKERRQ(ierr); 
    } else {
      ierr = temperature.get_record_years(t); CHKERRQ(ierr);
    }
    ierr = temperature.write(filename.c_str()); CHKERRQ(ierr);
  }

  if (set_contains(vars, "acab")) {
    if (enable_time_averaging) {
      ierr = mass_flux.average(t, dt); CHKERRQ(ierr); 
    } else {
      ierr = mass_flux.get_record_years(t); CHKERRQ(ierr);
    }
    ierr = mass_flux.write(filename.c_str()); CHKERRQ(ierr);
  }

  return 0;
}


PetscErrorCode PSDirectForcing::ice_surface_mass_flux(PetscReal t_years, PetscReal dt_years,
                                                      IceModelVec2S &result) {
  PetscErrorCode ierr;

  ierr = update(t_years, dt_years); CHKERRQ(ierr); 

  // "Periodize" the climate:
  t_years = my_mod(t_years);

  if (enable_time_averaging) {
    ierr = mass_flux.average(t_years, dt_years); CHKERRQ(ierr);
  } else {
    ierr = mass_flux.get_record_years(t_years); CHKERRQ(ierr);
  }

  ierr = mass_flux.copy_to(result); CHKERRQ(ierr); 

  if (bc_acab_lapse_rate_set) {

    PetscReal delta_acab;
    ierr = vH->begin_access(); CHKERRQ(ierr);
    ierr = result.begin_access();   CHKERRQ(ierr);
    ierr = surface->begin_access();   CHKERRQ(ierr);
    ierr = bc_surface.begin_access();   CHKERRQ(ierr);
    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
	if ((*vH)(i,j) > 0) { // only correct artm if there is ice
	  delta_acab = bc_acab_lapse_rate/1000/secpera * ((*surface)(i,j)-bc_surface(i,j))  ; 
	  result(i,j) += delta_acab;
	  ierr = verbPrintf(4, grid.com,"delta_artm=%f, bc_surf=%f, h=%f, artm=%f\n",delta_acab,bc_surface(i,j),(*surface)(i,j),result(i,j)); CHKERRQ(ierr);
	}
      }
    }
    ierr = bc_surface.end_access();   CHKERRQ(ierr);  
    ierr = surface->end_access();   CHKERRQ(ierr);
    ierr = result.end_access();   CHKERRQ(ierr);    
    ierr = vH->begin_access(); CHKERRQ(ierr);
  }


  return 0;
}


PetscErrorCode PSDirectForcing::ice_surface_temperature(PetscReal t_years, PetscReal dt_years,
                                                        IceModelVec2S &result) {
  PetscErrorCode ierr;

  ierr = update(t_years, dt_years); CHKERRQ(ierr); 

  // "Periodize" the climate:
  t_years = my_mod(t_years);

  if (enable_time_averaging) {
    ierr = temperature.average(t_years, dt_years); CHKERRQ(ierr);
    ierr = bc_surface.average(t_years, dt_years); CHKERRQ(ierr);
  } else {
    ierr = temperature.get_record_years(t_years); CHKERRQ(ierr);
    ierr = bc_surface.get_record_years(t_years); CHKERRQ(ierr);
  }

  ierr = temperature.copy_to(result); CHKERRQ(ierr); 

  if (bc_artm_lapse_rate_set) {

    PetscReal delta_artm;
    ierr = vH->begin_access(); CHKERRQ(ierr);
    ierr = result.begin_access();   CHKERRQ(ierr);
    ierr = surface->begin_access();   CHKERRQ(ierr);
    ierr = bc_surface.begin_access();   CHKERRQ(ierr);
    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
	if ((*vH)(i,j) > 0) { // only correct artm if there is ice
	  delta_artm = bc_artm_lapse_rate/1000 * ((*surface)(i,j)-bc_surface(i,j))  ; 
	  result(i,j) += delta_artm;
	  ierr = verbPrintf(4, grid.com,"delta_artm=%f, bc_surf=%f, h=%f, artm=%f\n",delta_artm,bc_surface(i,j),(*surface)(i,j),result(i,j)); CHKERRQ(ierr);
	}
      }
    }
    ierr = bc_surface.end_access();   CHKERRQ(ierr);  
    ierr = surface->end_access();   CHKERRQ(ierr);
    ierr = result.end_access();   CHKERRQ(ierr);    
    ierr = vH->begin_access(); CHKERRQ(ierr);
  }


  return 0;
}


