# include "pgrn_atmosphere.hh"

PA_EISMINT_Greenland::PA_EISMINT_Greenland(IceGrid &g, const NCConfigVariable &conf, PISMVars &vars)
  : PAFausto(g, conf, vars) {
  do_greenhouse_warming = false;
  greenhouse_warming_start_year = 0.0;
}

PetscErrorCode PA_EISMINT_Greenland::mean_annual_temp(PetscReal t_years, PetscReal dt_years,
						      IceModelVec2 &result) {
  PetscErrorCode ierr;
  ierr = update(t_years, dt_years); CHKERRQ(ierr);

  ierr = temp_ma.copy_to(result); CHKERRQ(ierr);
  ierr = result.set_attr("history",
			 "computed using a mean annual air temperature parameterization in " +
			 reference + "\n"); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode PA_EISMINT_Greenland::greenhouse_warming(PetscReal start_year) {
  do_greenhouse_warming = true;
  greenhouse_warming_start_year = start_year;
  return 0;
}

PetscErrorCode PA_EISMINT_Greenland::update(PetscReal t_years, PetscReal dt_years) {
  PetscErrorCode ierr;

  if ((gsl_fcmp(t_years,  t,  1e-4) == 0) &&
      (gsl_fcmp(dt_years, dt, 1e-4) == 0))
    return 0;

  t  = t_years;
  dt = dt_years;

  ierr = surfelev->begin_access();   CHKERRQ(ierr);
  ierr = lat->begin_access(); CHKERRQ(ierr);
  ierr = temp_ma.begin_access();  CHKERRQ(ierr);
  ierr = temp_mj.begin_access();  CHKERRQ(ierr);
  for (PetscInt i = grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j<grid.ys+grid.ym; ++j) {
      PetscReal Z = PetscMax((*surfelev)(i,j), 20 * ((*lat)(i,j) - 65));

      temp_ma(i,j) = 49.13 - 0.007992 * Z - 0.7576 * (*lat)(i,j) + 273.15;
      temp_mj(i,j) = 273.15 + 30.38 - 0.006277 * (*surfelev)(i,j) - 0.3262 * (*lat)(i,j);
    }
  }  
  ierr = surfelev->end_access();   CHKERRQ(ierr);
  ierr = lat->end_access(); CHKERRQ(ierr);
  ierr = temp_ma.end_access();  CHKERRQ(ierr);
  ierr = temp_mj.end_access();  CHKERRQ(ierr);

  if (do_greenhouse_warming == PETSC_TRUE) {
    const PetscScalar shift = greenhouse_shift(t_years, dt_years);
    ierr = temp_ma.shift(shift); CHKERRQ(ierr);
    ierr = temp_mj.shift(shift); CHKERRQ(ierr);
  }

  return 0;
}

PetscErrorCode PA_EISMINT_Greenland::init() {
  PetscErrorCode ierr;
  LocalInterpCtx *lic = NULL;
  bool regrid = false;
  int start = -1;

  ierr = verbPrintf(2, grid.com,
		    "* Initializing Greenland atmosphere model based on the EISMINT Greenland (C. Ritz, 1997)\n"
		    "  air temperature parameterization and using stored time-independent precipitation...\n"); CHKERRQ(ierr);

  reference = "Ritz, C. (1997). EISMINT Intercomparison Experiment: Comparison of existing Greenland models."
    " URL: http://homepages.vub.ac.be/~phuybrec/eismint/greenland.html";

  // Allocate internal IceModelVecs:
  ierr = temp_ma.create(grid, "eismint_temp_ma", false); CHKERRQ(ierr);
  ierr = temp_ma.set_attrs("climate_state",
			   "mean annual near-surface air temperature",
			   "K", 
			   ""); CHKERRQ(ierr);  // no CF standard_name ??
  ierr = temp_ma.set_attr("source", reference);

  ierr = temp_mj.create(grid, "eismint_temp_mj", false); CHKERRQ(ierr);
  ierr = temp_mj.set_attrs("climate_state",
			   "mean July near-surface air temperature",
			   "Kelvin",
			   ""); CHKERRQ(ierr);  // no CF standard_name ??
  ierr = temp_mj.set_attr("source", reference);

  ierr = snowprecip.create(grid, "snowprecip", false); CHKERRQ(ierr);
  ierr = snowprecip.set_attrs("climate_state", 
			      "mean annual ice-equivalent snow precipitation rate",
			      "m s-1", 
			      ""); CHKERRQ(ierr); // no CF standard_name ??
  ierr = snowprecip.set_glaciological_units("m year-1");
  snowprecip.write_in_glaciological_units = true;
  snowprecip.time_independent = true;

  // initialize pointers to fields the parameterization depends on:
  surfelev = dynamic_cast<IceModelVec2*>(variables.get("surface_altitude"));
  if (!surfelev) SETERRQ(1, "ERROR: surface_altitude is not available");

  lat = dynamic_cast<IceModelVec2*>(variables.get("latitude"));
  if (!lat) SETERRQ(1, "ERROR: latitude is not available");

  ierr = find_pism_input(snowprecip_filename, lic, regrid, start); CHKERRQ(ierr);

  // read snow precipitation rate from file
  ierr = verbPrintf(2, grid.com, 
		    "    reading mean annual ice-equivalent snow precipitation rate 'snowprecip'\n"
		    "      from %s ... \n",
		    snowprecip_filename.c_str()); CHKERRQ(ierr); 
  if (regrid) {
    ierr = snowprecip.regrid(snowprecip_filename.c_str(), *lic, true); CHKERRQ(ierr); // fails if not found!
  } else {
    ierr = snowprecip.read(snowprecip_filename.c_str(), start); CHKERRQ(ierr); // fails if not found!
  }
  string snowprecip_history = "read mean annual ice-equivalent snow precipitation rate from " +
    snowprecip_filename + "\n";

  ierr = snowprecip.set_attr("history", snowprecip_history); CHKERRQ(ierr);

  delete lic;

  t = grid.year;
  dt = 0;

  return 0;
}

PetscReal PA_EISMINT_Greenland::greenhouse_shift(PetscReal t_years, PetscReal dt_years) {
  // compute age back to start of GWL3; use midpoint of interval as the time
  PetscScalar age = (t_years + 0.5 * dt_years) - greenhouse_warming_start_year;
  if (age <= 0.0) {
    return 0.0; // before time 0.0, no warming
  } else {
    if (age <= 80.0) {
      return age * 0.035;
    } else if (age <= 500.0) {
      return 2.8 + (age - 80.0) * 0.0017;
    } else { // after 500 years, constant amount
      return 3.514;
    }
  }
}

PetscErrorCode PA_EISMINT_Greenland::temp_snapshot(PetscReal t_years, PetscReal dt_years,
				       IceModelVec2 &result) {
  PetscErrorCode ierr;

  ierr = PAFausto::temp_snapshot(t_years, dt_years, result); CHKERRQ(ierr);

  string history = "computed using the standard (cosine) yearly cycle\n";
  ierr = result.set_attr("history", history); CHKERRQ(ierr);

  return 0;
}
