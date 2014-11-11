// Copyright (C) 2010, 2011, 2012, 2013, 2014 Constantine Khroulev
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
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

#include "PBLingleClark.hh"
#include "PIO.hh"
#include "PISMTime.hh"
#include "IceGrid.hh"
#include "pism_options.hh"
#include "PISMConfig.hh"

#include <stdexcept>
#include "error_handling.hh"

namespace pism {

PBLingleClark::PBLingleClark(IceGrid &g, const Config &conf)
  : BedDef(g, conf) {

  if (allocate() != 0) {
    throw std::runtime_error("PBLingleClark allocation failed");
  }

}

PBLingleClark::~PBLingleClark() {

  if (deallocate() != 0) {
    PetscPrintf(grid.com, "PBLingleClark::~PBLingleClark(...): deallocate() failed\n");
  }

}

PetscErrorCode PBLingleClark::allocate() {
  PetscErrorCode ierr;

  ierr = topg_initial.allocate_proc0_copy(Hp0); CHKERRQ(ierr);
  ierr = topg_initial.allocate_proc0_copy(bedp0); CHKERRQ(ierr);
  ierr = topg_initial.allocate_proc0_copy(Hstartp0); CHKERRQ(ierr);
  ierr = topg_initial.allocate_proc0_copy(bedstartp0); CHKERRQ(ierr);
  ierr = topg_initial.allocate_proc0_copy(upliftp0); CHKERRQ(ierr);

  bool use_elastic_model = config.get_flag("bed_def_lc_elastic_model");

  if (grid.rank == 0) {
    ierr = bdLC.settings(config, use_elastic_model,
                         grid.Mx, grid.My, grid.dx, grid.dy,
                         4,     // use Z = 4 for now; to reduce global drift?
                         &Hstartp0, &bedstartp0, &upliftp0, &Hp0, &bedp0);
    CHKERRQ(ierr);

    ierr = bdLC.alloc(); CHKERRQ(ierr);
  }

  return 0;
}

PetscErrorCode PBLingleClark::deallocate() {
  PetscErrorCode ierr;

  ierr = VecDestroy(&Hp0);
  PISM_PETSC_CHK(ierr, "VecDestroy");
  ierr = VecDestroy(&bedp0);
  PISM_PETSC_CHK(ierr, "VecDestroy");
  ierr = VecDestroy(&Hstartp0);
  PISM_PETSC_CHK(ierr, "VecDestroy");
  ierr = VecDestroy(&bedstartp0);
  PISM_PETSC_CHK(ierr, "VecDestroy");
  ierr = VecDestroy(&upliftp0);
  PISM_PETSC_CHK(ierr, "VecDestroy");

  return 0;
}

//! Initialize the Lingle-Clark bed deformation model using uplift.
PetscErrorCode PBLingleClark::init(Vars &vars) {
  PetscErrorCode ierr;

  ierr = BedDef::init(vars); CHKERRQ(ierr);

  ierr = verbPrintf(2, grid.com,
                    "* Initializing the Lingle-Clark bed deformation model...\n"); CHKERRQ(ierr);

  ierr = correct_topg(); CHKERRQ(ierr);

  ierr = topg->copy_to(topg_last); CHKERRQ(ierr);

  ierr = thk->put_on_proc0(Hstartp0); CHKERRQ(ierr);
  ierr = topg->put_on_proc0(bedstartp0); CHKERRQ(ierr);
  ierr = uplift->put_on_proc0(upliftp0); CHKERRQ(ierr);

  if (grid.rank == 0) {
    ierr = bdLC.init(); CHKERRQ(ierr);
    ierr = bdLC.uplift_init(); CHKERRQ(ierr);
  }

  return 0;
}

PetscErrorCode PBLingleClark::correct_topg() {
  PetscErrorCode ierr;
  bool use_special_regrid_semantics, regrid_file_set, boot_file_set,
    topg_exists, topg_initial_exists, regrid_vars_set;
  std::string boot_filename, regrid_filename;
  PIO nc(grid, "guess_mode");

  ierr = OptionsIsSet("-regrid_bed_special",
                          "Correct topg when switching to a different grid",
                          use_special_regrid_semantics); CHKERRQ(ierr);

  // Stop if topg correction was not requiested.
  if (not use_special_regrid_semantics) {
    return 0;
  }

  ierr = OptionsString("-regrid_file", "Specifies the name of a file to regrid from",
                           regrid_filename, regrid_file_set); CHKERRQ(ierr);

  ierr = OptionsString("-boot_file", "Specifies the name of the file to bootstrap from",
                           boot_filename, boot_file_set); CHKERRQ(ierr);

  // Stop if it was requested, but we're not bootstrapping *and* regridding.
  if (not (regrid_file_set && boot_file_set)) {
    return 0;
  }

  nc.open(regrid_filename, PISM_READONLY);

  topg_initial_exists = nc.inq_var("topg_initial");
  topg_exists = nc.inq_var("topg");
  nc.close();

  // Stop if the regridding file does not have both topg and topg_initial.
  if (!(topg_initial_exists && topg_exists)) {
    return 0;
  }

  // Stop if the user asked to regrid topg (in this case no correction is necessary).
  std::set<std::string> regrid_vars;
  ierr = OptionsStringSet("-regrid_vars", "Specifies regridding variables", "",
                              regrid_vars, regrid_vars_set); CHKERRQ(ierr);

  if (regrid_vars_set) {
    if (set_contains(regrid_vars, "topg")) {
      ierr = verbPrintf(2, grid.com,
                        "  Bed elevation correction requested, but -regrid_vars contains topg...\n"); CHKERRQ(ierr);
      return 0;
    }
  }

  ierr = verbPrintf(2, grid.com,
                    "  Correcting topg from the bootstrapping file '%s' by adding the effect\n"
                    "  of the bed deformation from '%s'...\n",
                    boot_filename.c_str(), regrid_filename.c_str()); CHKERRQ(ierr);

  IceModelVec2S topg_tmp;       // will be de-allocated at 'return 0' below.
  const unsigned int WIDE_STENCIL = config.get("grid_max_stencil_width");
  ierr = topg_tmp.create(grid, "topg", WITH_GHOSTS, WIDE_STENCIL); CHKERRQ(ierr);
  ierr = topg_tmp.set_attrs("model_state", "bedrock surface elevation (at the end of the previous run)",
                            "m", "bedrock_altitude"); CHKERRQ(ierr);

  // Get topg and topg_initial from the regridding file.
  ierr = topg_initial.regrid(regrid_filename, CRITICAL); CHKERRQ(ierr);
  ierr =     topg_tmp.regrid(regrid_filename, CRITICAL); CHKERRQ(ierr);

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
PetscErrorCode PBLingleClark::update(double my_t, double my_dt) {
  PetscErrorCode ierr;

  if ((fabs(my_t - m_t)   < 1e-12) &&
      (fabs(my_dt - m_dt) < 1e-12))
    return 0;

  m_t  = my_t;
  m_dt = my_dt;

  double t_final = m_t + m_dt;

  // Check if it's time to update:
  double dt_beddef = t_final - t_beddef_last; // in seconds
  if ((dt_beddef < config.get("bed_def_interval_years", "years", "seconds") &&
       t_final < grid.time->end()) ||
      dt_beddef < 1e-12)
    return 0;

  t_beddef_last = t_final;

  ierr = thk->put_on_proc0(Hp0); CHKERRQ(ierr);
  ierr = topg->put_on_proc0(bedp0); CHKERRQ(ierr);

  if (grid.rank == 0) {  // only processor zero does the step
    ierr = bdLC.step(dt_beddef, // time step, in seconds
                     t_final - grid.time->start()); // time since the start of the run, in seconds
    CHKERRQ(ierr);
  }

  ierr = topg->get_from_proc0(bedp0); CHKERRQ(ierr);

  //! Finally, we need to update bed uplift and topg_last.
  ierr = compute_uplift(dt_beddef); CHKERRQ(ierr);
  ierr = topg->copy_to(topg_last); CHKERRQ(ierr);

  //! Increment the topg state counter. SIAFD relies on this!
  topg->inc_state_counter();

  return 0;
}

} // end of namespace pism
