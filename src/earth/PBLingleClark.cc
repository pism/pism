#if (PISM_HAVE_FFTW==1)

#include "PISMBedDef.hh"

PBLingleClark::PBLingleClark(IceGrid &g, const NCConfigVariable &conf)
  : PISMBedDef(g, conf) {
  PetscErrorCode ierr;

  if (grid.rank == 0) {
    ierr = bdLC.settings(config,
			 PETSC_FALSE, // turn off elastic model for now
			 grid.Mx, grid.My, grid.dx, grid.dy,
			 //                       2,                 // use Z = 2 for now
			 4,                 // use Z = 4 for now; to reduce global drift?
			 config.get("ice_density"),
			 &Hstartp0, &bedstartp0, &upliftp0, &Hp0, &bedp0);
    if (ierr != 0) {
      PetscPrintf(grid.com, "PBLingleClark::PBLingleClark(...): bdLC.settings() failed.\n");
      PetscEnd();
    }
  }

  ierr = allocate();
  if (ierr != 0) {
    PetscPrintf(grid.com, "PBLingleClark::PBLingleClark(...): allocate() failed\n");
    PetscEnd();
  }
  
}

PBLingleClark::~PBLingleClark() {
  PetscErrorCode ierr;

  ierr = deallocate();
  if (ierr != 0) {
    PetscPrintf(grid.com, "PBLingleClark::~PBLingleClark(...): deallocate() failed\n");
    PetscEnd();
  }

}

/* the following is from the PETSc FAQ page:

How do I collect all the values from a parallel PETSc vector into a vector 
on the zeroth processor?

    * Create the scatter context that will do the communication
          o VecScatterCreateToZero(v,&ctx,&w);
    * Actually do the communication; this can be done repeatedly as needed
          o VecScatterBegin(ctx,v,w,INSERT_VALUES,SCATTER_FORWARD);
          o VecScatterEnd(ctx,v,w,INSERT_VALUES,SCATTER_FORWARD);
    * Remember to free the scatter context when no longer needed
          o VecScatterDestroy(ctx);

Note that this simply concatenates in the parallel ordering of the vector.
If you are using a vector from DACreateGlobalVector() you likely want to 
first call DAGlobalToNaturalBegin/End() to scatter the original vector into 
the natural ordering in a new global vector before calling 
VecScatterBegin/End() to scatter the natural vector onto process 0.
*/

PetscErrorCode PBLingleClark::transfer_to_proc0(IceModelVec2S &source, Vec result) {
  PetscErrorCode ierr;

  ierr = source.copy_to_global(g2);

  ierr = DAGlobalToNaturalBegin(grid.da2, g2, INSERT_VALUES, g2natural); CHKERRQ(ierr);
  ierr =   DAGlobalToNaturalEnd(grid.da2, g2, INSERT_VALUES, g2natural); CHKERRQ(ierr);

  ierr = VecScatterBegin(scatter, g2natural, result, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
  ierr =   VecScatterEnd(scatter, g2natural, result, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PBLingleClark::transfer_from_proc0(Vec source, IceModelVec2S &result) {
  PetscErrorCode ierr;

  ierr = VecScatterBegin(scatter, source, g2natural, INSERT_VALUES, SCATTER_REVERSE); CHKERRQ(ierr);
  ierr =   VecScatterEnd(scatter, source, g2natural, INSERT_VALUES, SCATTER_REVERSE); CHKERRQ(ierr);

  ierr = DANaturalToGlobalBegin(grid.da2, g2natural, INSERT_VALUES, g2); CHKERRQ(ierr);
  ierr =   DANaturalToGlobalEnd(grid.da2, g2natural, INSERT_VALUES, g2); CHKERRQ(ierr);

  ierr = result.copy_from_global(g2); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PBLingleClark::allocate() {
  PetscErrorCode ierr;

  ierr = DACreateGlobalVector(grid.da2, &g2); CHKERRQ(ierr);

  // note we want a global Vec but reordered in the natural ordering so when it is
  // scattered to proc zero it is not all messed up; see above
  ierr = DACreateNaturalVector(grid.da2, &g2natural); CHKERRQ(ierr);
  // next get context *and* allocate samplep0 (on proc zero only, naturally)
  ierr = VecScatterCreateToZero(g2natural, &scatter, &Hp0); CHKERRQ(ierr);

  ierr = VecDuplicate(Hp0,&bedp0); CHKERRQ(ierr);
  ierr = VecDuplicate(Hp0,&Hstartp0); CHKERRQ(ierr);
  ierr = VecDuplicate(Hp0,&bedstartp0); CHKERRQ(ierr);
  ierr = VecDuplicate(Hp0,&upliftp0); CHKERRQ(ierr);

  if (grid.rank == 0) {
    ierr = bdLC.alloc(); CHKERRQ(ierr);
  }

  return 0;
}

PetscErrorCode PBLingleClark::deallocate() {
  PetscErrorCode ierr;

  ierr = VecDestroy(g2); CHKERRQ(ierr);
  ierr = VecDestroy(g2natural); CHKERRQ(ierr);
  ierr = VecScatterDestroy(scatter); CHKERRQ(ierr);

  ierr = VecDestroy(Hp0); CHKERRQ(ierr);
  ierr = VecDestroy(bedp0); CHKERRQ(ierr);
  ierr = VecDestroy(Hstartp0); CHKERRQ(ierr);
  ierr = VecDestroy(bedstartp0); CHKERRQ(ierr);
  ierr = VecDestroy(upliftp0); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PBLingleClark::init(PISMVars &vars) {
  PetscErrorCode ierr;

  thk = dynamic_cast<IceModelVec2S*>(vars.get("land_ice_thickness"));
  if (!thk) SETERRQ(1, "ERROR: thk is not available");

  topg = dynamic_cast<IceModelVec2S*>(vars.get("bedrock_altitude"));
  if (!topg) SETERRQ(1, "ERROR: topg is not available");

  dt_beddef = GSL_NAN;

  ierr = transfer_to_proc0(*thk,  Hstartp0); CHKERRQ(ierr);
  ierr = transfer_to_proc0(*topg, bedstartp0); CHKERRQ(ierr);
  ierr = transfer_to_proc0(dbdt,  upliftp0); CHKERRQ(ierr);

  ierr = bdLC.uplift_init(); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PBLingleClark::update(PetscReal t_years, PetscReal /*dt_years*/) {
  PetscErrorCode ierr;

  double bed_def_interval_years = config.get("bed_def_interval_years");

  // FIXME!
  if ((fabs(t_years - t)   < 1e-12) &&
      (fabs(dt_years - dt) < 1e-12))
    return 0;

  

  t  = t_years;
  dt = dt_years;

  return 0;
}

#else // PISM_HAVE_FFTW
#error "PISM build system error: Lingle and Clark bed deformation model requires FFTW3."
#endif
