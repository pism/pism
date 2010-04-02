#include "PISMBedDef.hh"

PISMBedDef::PISMBedDef(IceGrid &g, const NCConfigVariable &conf)
  : PISMComponent(g, conf) {

  thk    = NULL;
  topg   = NULL;
  uplift = NULL;

  t_beddef_last = GSL_NAN;

  PetscErrorCode ierr = pismbeddef_allocate();
  if (ierr != 0) {
    PetscPrintf(grid.com, "PISMBedDef::PISMBedDef(...): pismbeddef_allocate() failed\n");
    PetscEnd();
  }
}

PetscErrorCode PISMBedDef::pismbeddef_allocate() {
  PetscErrorCode ierr;
  PetscInt WIDE_STENCIL = 2;

  ierr = topg_last.create(grid, "topg", true, WIDE_STENCIL); CHKERRQ(ierr);
  ierr = topg_last.set_attrs("model_state", "bedrock surface elevation",
			     "m", "bedrock_altitude"); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PISMBedDef::init(PISMVars &vars) {

  t_beddef_last = grid.year;

  thk = dynamic_cast<IceModelVec2S*>(vars.get("land_ice_thickness"));
  if (!thk) SETERRQ(1, "ERROR: thk is not available");

  topg = dynamic_cast<IceModelVec2S*>(vars.get("bedrock_altitude"));
  if (!topg) SETERRQ(1, "ERROR: topg is not available");

  uplift = dynamic_cast<IceModelVec2S*>(vars.get("tendency_of_bedrock_altitude"));
  if (!uplift) SETERRQ(1, "ERROR: uplift is not available");
  
  return 0;
}

//! Compute bed uplift.
PetscErrorCode PISMBedDef::compute_uplift(PetscScalar dt_beddef) {
  PetscErrorCode ierr;

  ierr = topg->add(-1, topg_last, *uplift); CHKERRQ(ierr);
  //! uplift = (topg - topg_last) / dt
  ierr = uplift->scale(1.0 / (dt_beddef * secpera)); CHKERRQ(ierr); 

  return 0;
}
