
#include "PISMAtmosphere.hh"

PetscErrorCode PALapseRates::init(PISMVars &vars) {
  PetscErrorCode ierr;
  LocalInterpCtx *lic = NULL;
  bool regrid = false;
  int start = -1;

  ierr = PAConstant::init(vars); CHKERRQ(ierr);
	    
  ierr = verbPrintf(2, grid.com, "  NOTE: Using a lapse-rate correction for air temperature...\n"); CHKERRQ(ierr);

  usurf = dynamic_cast<IceModelVec2S*>(vars.get("surface_altitude"));
  if (!usurf) { SETERRQ(1, "ERROR: Surface elevation is not available"); }

  ierr = PetscOptionsBegin(grid.com, "", "Air temp. lapse rate model options", ""); CHKERRQ(ierr);
  {
    bool flag;
    ierr = PISMOptionsReal("-lapse_rate", "Air temperature lapse rate, degrees K per meter",
			   gamma, flag); CHKERRQ(ierr);
  }
  ierr = PetscOptionsEnd(); CHKERRQ(ierr);

  ierr = f.create(grid, "usurf", false); CHKERRQ(ierr);
  ierr = f.set_attrs("internal", "ice upper surface elevation",
		     "m", "surface_altitude"); CHKERRQ(ierr);

  // find (again) the PISM input file to read data from:

  ierr = find_pism_input(input_file, lic, regrid, start); CHKERRQ(ierr);

  // read snow precipitation rate and temperatures from file
  ierr = verbPrintf(2, grid.com, 
		    "    reading surface elevation 'usurf' from %s...\n",
		    input_file.c_str()); CHKERRQ(ierr); 
  if (regrid) {
    ierr = f.regrid(input_file.c_str(), *lic, true); CHKERRQ(ierr); // fails if not found!
  } else {
    ierr = f.read(input_file.c_str(), start); CHKERRQ(ierr); // fails if not found!
  }

  ierr = f.set_name("f"); CHKERRQ(ierr);
  ierr = f.set_attrs("internal", "initial condition for the lapse rate equation",
		     "K", ""); CHKERRQ(ierr);

  delete lic;

  ierr = f.scale(gamma); CHKERRQ(ierr);
  ierr = f.add(1.0, temperature); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PALapseRates::mean_annual_temp(PetscReal /*t_years*/, PetscReal /*dt_years*/,
					      IceModelVec2S &result) {
  PetscErrorCode ierr;

  ierr = result.begin_access(); CHKERRQ(ierr);
  ierr = f.begin_access(); CHKERRQ(ierr);
  ierr = usurf->begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      result(i,j) = -gamma * (*usurf)(i,j) + f(i,j);
    }
  }
  ierr = result.end_access(); CHKERRQ(ierr);
  ierr = f.end_access(); CHKERRQ(ierr);
  ierr = usurf->end_access(); CHKERRQ(ierr);

  return 0;
}
 
PetscErrorCode PALapseRates::begin_pointwise_access() {
  PetscErrorCode ierr;
  ierr = PAConstant::begin_pointwise_access(); CHKERRQ(ierr);
  ierr = usurf->begin_access(); CHKERRQ(ierr);
  ierr = f.begin_access(); CHKERRQ(ierr);
  return 0;
}
 
PetscErrorCode PALapseRates::end_pointwise_access() {
  PetscErrorCode ierr;
  ierr = PAConstant::end_pointwise_access(); CHKERRQ(ierr);
  ierr = usurf->end_access(); CHKERRQ(ierr);
  ierr = f.end_access(); CHKERRQ(ierr);
  return 0;
}
   
PetscErrorCode PALapseRates::temp_time_series(int i, int j, int N,
					      PetscReal */*ts*/, PetscReal *values) {
  for (int k = 0; k < N; ++k)
    values[k] = -gamma * (*usurf)(i,j) + f(i,j);

  return 0;
}

PetscErrorCode PALapseRates::write_model_state(PetscReal t_years, PetscReal dt_years,
						string filename) {
  PetscErrorCode ierr;

  ierr = snowprecip.write(filename.c_str()); CHKERRQ(ierr);

  IceModelVec2S temp_ma;
  ierr = temp_ma.create(grid, "airtemp_ma", false); CHKERRQ(ierr); // FIXME! choose the right name
  ierr = temp_ma.set_attrs(
            "climate_state",
            "mean annual near-surface air temperature",
            "K",
	    ""); CHKERRQ(ierr);

  ierr = temp_snapshot(t_years, dt_years, temp_ma); CHKERRQ(ierr);

  ierr = temp_ma.write(filename.c_str()); CHKERRQ(ierr);

  return 0;
}

