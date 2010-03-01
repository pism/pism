
PISMBedDef::PISMBedDef(IceGrid &g, const NCConfigVariable &conf)
  : PISMComponent(g, conf) {

  thk = NULL;
  topg = NULL;

  PetscErrorCode ierr = pismbeddef_allocate();
  if (ierr != 0) {
    PetscPrintf(grid.com, "PISMBedDef::PISMBedDef(...): pismbeddef_allocate() failed\n");
    PetscEnd();
  }
}

PetscErrorCode PISMBedDef::pismbeddef_allocate() {
  PetscErrorCode ierr;
  // These IceModelVecs are automatically de-allocated by the PISMBedDef
  // destructor.

  // "local" because topg is local (this makes some code simpler)
  ierr = dtopgdt.create(grid, "dbdt", true); CHKERRQ(ierr);
  ierr = dtopgdt.set_attrs("model_state", "bedrock uplift rate",
			   "m s-1", "tendency_of_bedrock_altitude"); CHKERRQ(ierr);
  ierr = dtopgdt.set_glaciological_units("m year-1");
  dtopgdt.write_in_glaciological_units = true;
  
  ierr = thk_last(grid, "beddef_thk_last", true); CHKERRQ(ierr);
  // attributes are not set because this field is never read or written

  ierr = topg_last(grid, "topg", true); CHKERRQ(ierr);
  ierr = topg_last.set_attrs("model_state", "bedrock surface elevation",
			     "m", "bedrock_altitude"); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PISMBedDef::write_fields(set<string> vars, PetscReal t_years,
					PetscReal dt_years, string filename) {
  PetscErrorCode ierr;

  ierr = update(t_years, 0); CHKERRQ(ierr);

  if (vars.find("dbdt") != vars.end()) {
    ierr = dtopgdt.write(filename.c_str()); CHKERRQ(ierr);
  }

  if (vars.find("topg") != vars.end()) {
    ierr = topg_last.write(filename.c_str()); CHKERRQ(ierr);
  }

  return 0;
}
