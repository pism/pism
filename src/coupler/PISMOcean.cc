// Copyright (C) 2008-2011 Ed Bueler, Constantine Khroulev, Ricarda Winkelmann,
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

#include "PISMOcean.hh"
POConstant::POConstant(IceGrid &g, const NCConfigVariable &conf)
  : PISMOceanModel(g, conf) {

  shelfbmassflux.init_2d("shelfbmassflux", g);
  shelfbmassflux.set_string("pism_intent", "climate_state");
  shelfbmassflux.set_string("long_name",
                            "ice mass flux from ice shelf base (positive flux is loss from ice shelf)");
  shelfbmassflux.set_units("m s-1");
  shelfbmassflux.set_glaciological_units("m year-1");

  shelfbtemp.init_2d("shelfbtemp", g);
  shelfbtemp.set_string("pism_intent", "climate_state");
  shelfbtemp.set_string("long_name",
                        "absolute temperature at ice shelf base");
  shelfbtemp.set_units("Kelvin");
}

PetscErrorCode POConstant::init(PISMVars &vars) {
  PetscErrorCode ierr;

  if (!config.get_flag("is_dry_simulation")) {
    ierr = verbPrintf(2, grid.com, "* Initializing the constant ocean model...\n"); CHKERRQ(ierr);
  }

  ice_thickness = dynamic_cast<IceModelVec2S*>(vars.get("land_ice_thickness"));
  if (!ice_thickness) { SETERRQ(1, "ERROR: ice thickness is not available"); }

  return 0;
}

PetscErrorCode POConstant::sea_level_elevation(PetscReal /*t_years*/, PetscReal /*dt_years*/,
					       PetscReal &result) {
  result = sea_level;
  return 0;
}

PetscErrorCode POConstant::shelf_base_temperature(PetscReal /*t_years*/, PetscReal /*dt_years*/,
						  IceModelVec2S &result) {
  PetscErrorCode ierr;

  const PetscScalar T0 = config.get("water_triple_point_temperature"), // K
    beta_CC_grad = config.get("beta_CC") * config.get("ice_density") * config.get("standard_gravity"), // K m-1
    ice_rho = config.get("ice_density"),
    sea_water_rho = config.get("sea_water_density");

  PetscScalar **H;
  ierr = ice_thickness->get_array(H);   CHKERRQ(ierr);
  ierr = result.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      const PetscScalar shelfbaseelev = - ( ice_rho / sea_water_rho ) * H[i][j]; // FIXME task #7297
      // temp is set to melting point at depth
      result(i,j) = T0 + beta_CC_grad * shelfbaseelev;  // base elev negative here so is below T0
    }
  }
  ierr = ice_thickness->end_access(); CHKERRQ(ierr);
  ierr = result.end_access(); CHKERRQ(ierr);
  
  return 0;                                 
}

//! Computes mass flux in ice-equivalent m s-1, from assumption that basal heat flux rate converts to mass flux.
PetscErrorCode POConstant::shelf_base_mass_flux(PetscReal /*t_years*/, PetscReal /*dt_years*/,
						IceModelVec2S &result) {
  PetscErrorCode ierr;
  PetscReal L = config.get("water_latent_heat_fusion"),
    rho = config.get("ice_density"),
    // following has units:   J m-2 s-1 / (J kg-1 * kg m-3) = m s-1
    meltrate = config.get("ocean_sub_shelf_heat_flux_into_ice") / (L * rho); // m s-1

  ierr = result.set(meltrate); CHKERRQ(ierr);

  return 0;
}

void POConstant::add_vars_to_output(string keyword, set<string> &result) {
  if (keyword != "small") {
    result.insert("shelfbtemp");
    result.insert("shelfbmassflux");
  }
}

PetscErrorCode POConstant::define_variables(set<string> vars, const NCTool &nc,
                                            nc_type nctype) {
  PetscErrorCode ierr;
  int varid;

  if (set_contains(vars, "shelfbtemp")) {
    ierr = shelfbtemp.define(nc, varid, nctype, true); CHKERRQ(ierr);
  }

  if (set_contains(vars, "shelfbmassflux")) {
    ierr = shelfbmassflux.define(nc, varid, nctype, true); CHKERRQ(ierr);
  }

  return 0;
}

PetscErrorCode POConstant::write_variables(set<string> vars, string filename) {
  PetscErrorCode ierr;
  IceModelVec2S tmp;

  if (set_contains(vars, "shelfbtemp")) {
    if (!tmp.was_created()) {
      ierr = tmp.create(grid, "tmp", false); CHKERRQ(ierr);
    }

    ierr = tmp.set_metadata(shelfbtemp, 0); CHKERRQ(ierr);
    ierr = shelf_base_temperature(0, 0, tmp); CHKERRQ(ierr);
    ierr = tmp.write(filename.c_str()); CHKERRQ(ierr);
  }

  if (set_contains(vars, "shelfbmassflux")) {
    if (!tmp.was_created()) {
      ierr = tmp.create(grid, "tmp", false); CHKERRQ(ierr);
    }

    ierr = tmp.set_metadata(shelfbmassflux, 0); CHKERRQ(ierr);
    ierr = shelf_base_mass_flux(0, 0, tmp); CHKERRQ(ierr);
    ierr = tmp.write(filename.c_str()); CHKERRQ(ierr);
  }

  return 0;
}


///// POModifier

void POModifier::attach_input(PISMOceanModel *input) {
  if (input_model != NULL)
    delete input_model;

  input_model = input;
}

///// POForcing

PetscErrorCode POForcing::init(PISMVars &vars) {
  PetscErrorCode ierr;
  bool optSet;
  string dSLfile;

  ierr = input_model->init(vars); CHKERRQ(ierr);

  ierr = verbPrintf(2, grid.com, "* Initializing sea level forcing...\n"); CHKERRQ(ierr);

  ierr = PetscOptionsBegin(grid.com, "", "Sea level forcing options", ""); CHKERRQ(ierr);
  {
    ierr = PISMOptionsString("-dSLforcing", "Sea level forcing file",
			     dSLfile, optSet); CHKERRQ(ierr);
  }
  ierr = PetscOptionsEnd(); CHKERRQ(ierr);

  if (optSet == PETSC_TRUE) {

    dSLforcing = new Timeseries(grid.com, grid.rank, "delta_sea_level", "t");
    ierr = dSLforcing->set_units("m", ""); CHKERRQ(ierr);
    ierr = dSLforcing->set_dimension_units("years", ""); CHKERRQ(ierr);
    ierr = dSLforcing->set_attr("long_name", "sea level elevation offsets"); CHKERRQ(ierr);

    ierr = verbPrintf(2, grid.com, 
		      "  reading delta sea level data from forcing file %s ...\n", 
		      dSLfile.c_str()); CHKERRQ(ierr);
    ierr = dSLforcing->read(dSLfile.c_str()); CHKERRQ(ierr);

    delta_sea_level = new DiagnosticTimeseries(grid.com, grid.rank, "delta_sea_level", "t");
    ierr = delta_sea_level->set_units("m", ""); CHKERRQ(ierr);
    ierr = delta_sea_level->set_dimension_units("years", ""); CHKERRQ(ierr);
    ierr = delta_sea_level->set_attr("long_name", "sea level elevation offsets"); CHKERRQ(ierr);
  } else {
    ierr = verbPrintf(2, grid.com, "NOTE: -dSLforcing option is not set. Forcing is inactive...\n"); CHKERRQ(ierr);
    PISMEnd();
  }

  return 0;
}

PetscErrorCode POForcing::sea_level_elevation(PetscReal t_years, PetscReal dt_years,
					      PetscReal &result) {
  PetscErrorCode ierr;
  double T = t_years + 0.5 * dt_years;

  ierr = input_model->sea_level_elevation(t_years, dt_years, result); CHKERRQ(ierr);

  if (dSLforcing != NULL)
    result += (*dSLforcing)(T);

  return 0;
}


PetscErrorCode POForcing::shelf_base_temperature(PetscReal t_years, PetscReal dt_years,
						 IceModelVec2S &result) {
  PetscErrorCode ierr;

  ierr = input_model->shelf_base_temperature(t_years, dt_years, result); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode POForcing::shelf_base_mass_flux(PetscReal t_years, PetscReal dt_years,
					       IceModelVec2S &result) {
  PetscErrorCode ierr;

  ierr = input_model->shelf_base_mass_flux(t_years, dt_years, result); CHKERRQ(ierr);

  return 0;
}

// PetscErrorCode POForcing::write_diagnostic_fields(PetscReal t_years, PetscReal dt_years,
// 						  string filename) {
//   PetscErrorCode ierr;
//   double T = t_years + 0.5 * dt_years;

//   if (dSLforcing != NULL) {
//     delta_sea_level->output_filename = filename;
//     ierr = delta_sea_level->append(T, (*dSLforcing)(T)); CHKERRQ(ierr);
//     ierr = delta_sea_level->interp(T); CHKERRQ(ierr);
//     ierr = delta_sea_level->flush(); CHKERRQ(ierr);
//   }

//   return 0;
// }
