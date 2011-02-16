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

#if (PISM_HAVE_FFTW==1)

#include "PISMBedDef.hh"
#include "PISMIO.hh"

PBLingleClark::PBLingleClark(IceGrid &g, const NCConfigVariable &conf)
  : PISMBedDef(g, conf) {
  PetscErrorCode ierr;

  ierr = allocate();
  if (ierr != 0) {
    PetscPrintf(grid.com, "PBLingleClark::PBLingleClark(...): allocate() failed\n");
    PISMEnd();
  }
  
}

PBLingleClark::~PBLingleClark() {
  PetscErrorCode ierr;

  ierr = deallocate();
  if (ierr != 0) {
    PetscPrintf(grid.com, "PBLingleClark::~PBLingleClark(...): deallocate() failed\n");
    PISMEnd();
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

PetscErrorCode PBLingleClark::transfer_to_proc0(IceModelVec2S *source, Vec result) {
  PetscErrorCode ierr;

  ierr = source->copy_to(g2);

  ierr = DAGlobalToNaturalBegin(grid.da2, g2, INSERT_VALUES, g2natural); CHKERRQ(ierr);
  ierr =   DAGlobalToNaturalEnd(grid.da2, g2, INSERT_VALUES, g2natural); CHKERRQ(ierr);

  ierr = VecScatterBegin(scatter, g2natural, result, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
  ierr =   VecScatterEnd(scatter, g2natural, result, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PBLingleClark::transfer_from_proc0(Vec source, IceModelVec2S *result) {
  PetscErrorCode ierr;

  ierr = VecScatterBegin(scatter, source, g2natural, INSERT_VALUES, SCATTER_REVERSE); CHKERRQ(ierr);
  ierr =   VecScatterEnd(scatter, source, g2natural, INSERT_VALUES, SCATTER_REVERSE); CHKERRQ(ierr);

  ierr = DANaturalToGlobalBegin(grid.da2, g2natural, INSERT_VALUES, g2); CHKERRQ(ierr);
  ierr =   DANaturalToGlobalEnd(grid.da2, g2natural, INSERT_VALUES, g2); CHKERRQ(ierr);

  ierr = result->copy_from(g2); CHKERRQ(ierr);

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
    ierr = bdLC.settings(config,
			 PETSC_FALSE, // turn off elastic model for now
			 grid.Mx, grid.My, grid.dx, grid.dy,
			 //                       2,                 // use Z = 2 for now
			 4,                 // use Z = 4 for now; to reduce global drift?
			 config.get("ice_density"),
			 &Hstartp0, &bedstartp0, &upliftp0, &Hp0, &bedp0);
    CHKERRQ(ierr);

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

//! Initialize the Lingle-Clark bed deformation model using uplift.
PetscErrorCode PBLingleClark::init(PISMVars &vars) {
  PetscErrorCode ierr;

  ierr = PISMBedDef::init(vars); CHKERRQ(ierr);

  ierr = verbPrintf(2, grid.com,
		    "* Initializing the Lingle-Clark bed deformation model...\n"); CHKERRQ(ierr);

  ierr = correct_topg(); CHKERRQ(ierr);

  ierr = transfer_to_proc0(thk,    Hstartp0); CHKERRQ(ierr);
  ierr = transfer_to_proc0(topg,   bedstartp0); CHKERRQ(ierr);
  ierr = transfer_to_proc0(uplift, upliftp0); CHKERRQ(ierr);

  if (grid.rank == 0) {
    ierr = bdLC.uplift_init(); CHKERRQ(ierr);
  }

  return 0;
}

PetscErrorCode PBLingleClark::correct_topg() {
  PetscErrorCode ierr;
  bool use_special_regrid_semantics, regrid_file_set, boot_file_set,
    topg_exists, topg_initial_exists, regrid_vars_set;
  string boot_filename, regrid_filename;
  PISMIO nc(&grid);

  ierr = PISMOptionsIsSet("-regrid_bed_special",
                          "Correct topg when switching to a different grid",
                          use_special_regrid_semantics); CHKERRQ(ierr);

  // Stop if topg correction was not requiested.
  if (!use_special_regrid_semantics) return 0;

  ierr = PISMOptionsString("-regrid_file", "Specifies the name of a file to regrid from",
                           regrid_filename, regrid_file_set); CHKERRQ(ierr);

  ierr = PISMOptionsString("-boot_file", "Specifies the name of the file to bootstrap from",
                           boot_filename, boot_file_set); CHKERRQ(ierr);

  // Stop if it was requested, but we're not bootstrapping *and* regridding.
  if (! (regrid_file_set && boot_file_set) ) return 0;

  ierr = nc.open_for_reading(regrid_filename.c_str()); CHKERRQ(ierr);

  ierr = nc.find_variable("topg_initial", NULL, topg_initial_exists); CHKERRQ(ierr);
  ierr = nc.find_variable("topg", NULL, topg_exists); CHKERRQ(ierr);

  // Stop if the regridding file does not have both topg and topg_initial.
  if ( !(topg_initial_exists && topg_exists) ) {
    ierr = nc.close(); CHKERRQ(ierr);
    return 0;
  }

  // Stop if the user asked to regrid topg (in this case no correction is necessary).
  vector<string> regrid_vars;
  ierr = PISMOptionsStringArray("-regrid_vars", "Specifies regridding variables", "",
                                regrid_vars, regrid_vars_set); CHKERRQ(ierr);

  if (regrid_vars_set) {
    for (unsigned int i = 0; i < regrid_vars.size(); ++i) {
      if (regrid_vars[i] == "topg") {
        ierr = verbPrintf(2, grid.com,
                          "  Bed elevation correction requested, but -regrid_vars contains topg...\n"); CHKERRQ(ierr); 
        ierr = nc.close(); CHKERRQ(ierr);
        return 0;
      }
    }
  }

  ierr = verbPrintf(2, grid.com,
                    "  Correcting topg from the bootstrapping file '%s' by adding the effect\n"
                    "  of the bed deformation from '%s'...\n",
                    boot_filename.c_str(), regrid_filename.c_str()); CHKERRQ(ierr); 
  
  IceModelVec2S topg_tmp;       // will be de-allocated at 'return 0' below.
  int WIDE_STENCIL = grid.max_stencil_width;
  ierr = topg_tmp.create(grid, "topg", true, WIDE_STENCIL); CHKERRQ(ierr);
  ierr = topg_tmp.set_attrs("model_state", "bedrock surface elevation (at the end of the previous run)",
                            "m", "bedrock_altitude"); CHKERRQ(ierr);
  grid_info g;
  ierr = nc.get_grid_info_2d(g); CHKERRQ(ierr);
  ierr = nc.close(); CHKERRQ(ierr);

  // Get topg and topg_initial from the regridding file.
  ierr = topg_initial.regrid(regrid_filename.c_str(), true); CHKERRQ(ierr);
  ierr =     topg_tmp.regrid(regrid_filename.c_str(), true); CHKERRQ(ierr);

  // After bootstrapping, topg contains the bed elevation field from
  // -boot_file.

  ierr = topg_tmp.add(-1.0, topg_initial); CHKERRQ(ierr); 
  // Now topg_tmp contains the change in bed elevation computed during the run
  // that produced -regrid_file.

  // Apply this change to topg from -boot_file:
  ierr = topg->add(1.0, topg_tmp); CHKERRQ(ierr);

  // Store the corrected topg as the new "topg_initial".
  ierr = topg->copy_to(topg_initial); CHKERRQ(ierr);

  return 0;
}


//! Update the Lingle-Clark bed deformation model.
PetscErrorCode PBLingleClark::update(PetscReal t_years, PetscReal dt_years) {
  PetscErrorCode ierr;

  if ((fabs(t_years - t)   < 1e-12) &&
      (fabs(dt_years - dt) < 1e-12))
    return 0;

  t  = t_years;
  dt = dt_years;

  // Check if it's time to update:
  PetscScalar dt_beddef = t_years - t_beddef_last;
  if (dt_beddef < config.get("bed_def_interval_years"))
    return 0;

  t_beddef_last = t_years;

  ierr = transfer_to_proc0(thk,  Hp0);   CHKERRQ(ierr);
  ierr = transfer_to_proc0(topg, bedp0); CHKERRQ(ierr);

  if (grid.rank == 0) {  // only processor zero does the step
    ierr = bdLC.step(dt_beddef, t_years - grid.start_year); CHKERRQ(ierr);
  }

  ierr = transfer_from_proc0(bedp0, topg); CHKERRQ(ierr);

  //! Finally, we need to update bed uplift and topg_last.
  ierr = compute_uplift(dt_beddef); CHKERRQ(ierr);
  ierr = topg->copy_to(topg_last); CHKERRQ(ierr);

  return 0;
}

#else // PISM_HAVE_FFTW
#error "PISM build system error: Lingle and Clark bed deformation model requires FFTW3."
#endif
