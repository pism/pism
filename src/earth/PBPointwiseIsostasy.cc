#include "PISMBedDef.hh"

PBPointwiseIsostasy::PBPointwiseIsostasy(IceGrid &g, const NCConfigVariable &conf) 
  : PISMBedDef(g, conf) {
  PetscErrorCode ierr;
  
}

PBPointwiseIsostasy::~PBPointwiseIsostasy() {
}

PetscErrorCode PBPointwiseIsostasy::init(PISMVars &vars) {
  PetscErrorCode ierr;

  thk = dynamic_cast<IceModelVec2*>(vars.get("land_ice_thickness"));
  if (!thk) SETERRQ(1, "ERROR: thk is not available");

  topg = dynamic_cast<IceModelVec2*>(vars.get("bedrock_altitude"));
  if (!topg) SETERRQ(1, "ERROR: topg is not available");

  dt_beddef = GSL_NAN;

  return 0;
}

PetscErrorCode PBPointwiseIsostasy::update(PetscReal t_years, PetscReal dt_years) {
  PetscErrorCode ierr;

  if ((fabs(t_years - t)   < 1e-12) &&
      (fabs(dt_years - dt) < 1e-12))
    return 0;

  t  = t_years;
  dt = dt_years;

  // FIXME
  PetscErrorCode ierr;
  IceModelVec2 vHdiff = vWork2d[0];

  double lithosphere_density = config.get("lithosphere_density");

  const PetscScalar  f = ice->rho / lithosphere_density;
  ierr = vH.add(-1, vHlast, vHdiff); CHKERRQ(ierr);  // Hdiff = H - Hlast
  ierr = vbedlast.add(-f, vHdiff, vbed); CHKERRQ(ierr);  // bed = bedlast - f (Hdiff)
  return 0;

  return 0;
}

PetscErrorCode PBPointwiseIsostasy::write_fields(set<string> vars, PetscReal t_years,
						 PetscReal dt_years, string filename) {
  PetscErrorCode ierr;

  return 0;
}

PetscErrorCode PBPointwiseIsostasy::write_diagnostic_fields(PetscReal t_years, PetscReal dt_years,
							    string filename) {
  PetscErrorCode ierr;

  return 0;
}

PetscErrorCode PBPointwiseIsostasy::write_model_state(PetscReal t_years, PetscReal dt_years,
						      string filename) {
  PetscErrorCode ierr;

  return 0;
}
