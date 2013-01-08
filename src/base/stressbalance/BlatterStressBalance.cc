// Copyright (C) 2010-2013 Ed Bueler and Constantine Khroulev
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

#include "BlatterStressBalance.hh"
#include "PISMOcean.hh"
#include "IceGrid.hh"
#include "PISMVars.hh"

/*!
 * FIXMEs:
 *
 * This setup does not allow depth-dependent anything.
 *
 * Periodic boundary conditions: we need a strength extension or something similar.
 *
 * We need to allow spatially-variable (and depth-dependent) ice hardness.
 *
 * We need to compute the basal frictional heating and the volumetric strain heating.
 *
 * We need to compute drag from tauc instead of just plugging it in (which is
 * wrong). Alternatively, we can adjust THIFriction(...) so that it has tauc as
 * an input (instead of beta0).
 *
 * We need to get rid of thi->alpha.
 *
 * The mesh_to_regular_grid() method needs to be fixed: it needs to interpolate
 * u and v from mesh nodes to regular grid locations.
 *
 * We need to coarsen instead of refining. (The finest multigrid level is the
 * one that needs to match the PISM grid.)
 *
 * We need to interpolate b (bed elevation), s (surface elevation) and beta0
 * (tauc) into meshes used in all the multigrid levels. (So far we can use only
 * one level.)
 *
 * We need to switch from using DMMG to SNESSetDM(...). See
 * src/snes/examples/tutorials/ex{5,19,50}.c.
 *
 * We need to get the scaling straight. (From MKS to scaled variables in
 * BlatterStressBalance::setup() and back in mesh_to_regular_grid().
 *
 * Right now the code in THI.cc has the no-slip boundary condition at the base
 * hard-wired. This needs to change.
 */


BlatterStressBalance::BlatterStressBalance(IceGrid &g, PISMOceanModel *ocean_model, const NCConfigVariable &conf)
  : PISMStressBalance(g,NULL,NULL,ocean_model,conf)
{
  if (allocate_blatter() != 0) {
    PetscPrintf(grid.com, "FATAL ERROR: BlatterStressBalance allocation failed.\n");
    PISMEnd();
  }
}

BlatterStressBalance::~BlatterStressBalance()
{
  if (deallocate_blatter() != 0) {
    PetscPrintf(grid.com, "FATAL ERROR: BlatterStressBalance deallocation failed.\n");
    PISMEnd();
  }
}

PetscErrorCode BlatterStressBalance::allocate_blatter() {
  PetscErrorCode ierr;
  DM da2;
  ierr = grid.get_dm(1, grid.max_stencil_width, da2); CHKERRQ(ierr);

  ierr = THICreate(grid.com,&thi);CHKERRQ(ierr);
  ierr = THISetup(grid.com, da2, 2*grid.Lx, 2*grid.Ly, thi, &dmmg); CHKERRQ(ierr);

  // now allocate u, v, and w:

  ierr =     u.create(grid, "uvel", true); CHKERRQ(ierr);
  ierr =     u.set_attrs("diagnostic", "horizontal velocity of ice in the X direction",
			  "m s-1", "land_ice_x_velocity"); CHKERRQ(ierr);
  ierr =     u.set_glaciological_units("m year-1"); CHKERRQ(ierr);
  u.write_in_glaciological_units = true;

  ierr =     v.create(grid, "vvel", true); CHKERRQ(ierr);
  ierr =     v.set_attrs("diagnostic", "horizontal velocity of ice in the Y direction",
			  "m s-1", "land_ice_y_velocity"); CHKERRQ(ierr);
  ierr =     v.set_glaciological_units("m year-1"); CHKERRQ(ierr);
  v.write_in_glaciological_units = true;

  ierr = vertically_averaged_velocity.create(grid, "bar", true); CHKERRQ(ierr); // components ubar, vbar
  ierr = vertically_averaged_velocity.set_attrs("model_state",
                                                "thickness-advective ice velocity (x-component)",
                                                "m s-1", "", 0); CHKERRQ(ierr);
  ierr = vertically_averaged_velocity.set_attrs("model_state",
                                                "thickness-advective ice velocity (y-component)",
                                                "m s-1", "", 1); CHKERRQ(ierr);
  ierr = vertically_averaged_velocity.set_glaciological_units("m year-1"); CHKERRQ(ierr);
  vertically_averaged_velocity.write_in_glaciological_units = true;

  ierr = basal_frictional_heating.create(grid, "bfrict", false); CHKERRQ(ierr);
  ierr = basal_frictional_heating.set_attrs("diagnostic",
                                            "basal frictional heating",
                                            "W m-2", ""); CHKERRQ(ierr);
  ierr = basal_frictional_heating.set_glaciological_units("mW m-2"); CHKERRQ(ierr);
  basal_frictional_heating.write_in_glaciological_units = true;

  // Sigma
  ierr = Sigma.create(grid, "strainheat", false); CHKERRQ(ierr); // never diff'ed in hor dirs
  ierr = Sigma.set_attrs("internal",
                          "rate of strain heating in ice (dissipation heating)",
	        	  "W m-3", ""); CHKERRQ(ierr);
  ierr = Sigma.set_glaciological_units("mW m-3"); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode BlatterStressBalance::deallocate_blatter() {
  PetscErrorCode ierr;

  ierr = DMMGDestroy(dmmg);CHKERRQ(ierr);
  ierr = THIDestroy(thi);CHKERRQ(ierr);

  return 0;
}


PetscErrorCode BlatterStressBalance::init(PISMVars &vars) {
  //PetscErrorCode ierr;

  variables = &vars;

  topg = dynamic_cast<IceModelVec2S*>(vars.get("bedrock_altitude"));
  if (topg == NULL) SETERRQ(grid.com, 1, "bedrock_altitude is not available");

  usurf = dynamic_cast<IceModelVec2S*>(vars.get("surface_altitude"));
  if (usurf == NULL) SETERRQ(grid.com, 1, "surface_altitude is not available");

  tauc = dynamic_cast<IceModelVec2S*>(vars.get("tauc"));
  if (tauc == NULL) SETERRQ(grid.com, 1, "tauc is not available");

  return 0;
}

PetscErrorCode BlatterStressBalance::update(bool fast) {
  PetscErrorCode ierr;

  if (fast) {
    ierr = verbPrintf(1,grid.com,
       "PISM ERROR:  'fast' mode not meaningful for BlatterStressBalance\n"
       "  ENDING ...\n\n"); CHKERRQ(ierr);
    PISMEnd();
  }

  // setup
  ierr =  setup(); CHKERRQ(ierr);

  // solve
  ierr = DMMGSolve(dmmg);CHKERRQ(ierr);

  // Transfer solution from the FEM mesh to the regular grid used in the rest
  // of PISM and compute the vertically-averaged velocity.
  ierr = mesh_to_regular_grid(); CHKERRQ(ierr);

  ierr = compute_vertical_velocity(&u, &v, basal_melt_rate, w); CHKERRQ(ierr);

  // ierr =  compute_basal_frictional_heating(); CHKERRQ(ierr);

  // ierr =  compute_volumetric_strain_heating(); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode BlatterStressBalance::setup() {
  PetscErrorCode ierr;
  PrmNode **parameters;
  DM da = dmmg[0]->dm;
  ierr =  DMDAGetPrmNodeArray(da, &parameters); CHKERRQ(ierr);

  ierr = topg->begin_access(); CHKERRQ(ierr);
  ierr = usurf->begin_access(); CHKERRQ(ierr);
  ierr = tauc->begin_access(); CHKERRQ(ierr);

  // So far this is done on the coarsest level only.
  for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
      parameters[i][j].b = (*topg)(i,j);
      parameters[i][j].h = (*usurf)(i,j);
      parameters[i][j].beta2 = (*tauc)(i,j);
    }
  }

  ierr = topg->end_access(); CHKERRQ(ierr);
  ierr = usurf->end_access(); CHKERRQ(ierr);
  ierr = tauc->end_access(); CHKERRQ(ierr);

  ierr = DMDARestorePrmNodeArray(da, &parameters); CHKERRQ(ierr);

  ierr = DMDAPrmNodeArrayCommBegin(da); CHKERRQ(ierr);
  ierr = DMDAPrmNodeArrayCommEnd(da); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode BlatterStressBalance::mesh_to_regular_grid() {
  PetscErrorCode ierr;

  PISMVector2 ***U;
  PetscScalar *u_ij, *v_ij;
  ierr = DMDAVecGetArray(dmmg[0]->dm, dmmg[0]->x, &U); CHKERRQ(ierr);

  ierr = u.begin_access(); CHKERRQ(ierr);
  ierr = v.begin_access(); CHKERRQ(ierr);

  for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
      ierr = u.getInternalColumn(i,j,&u_ij); CHKERRQ(ierr);
      ierr = v.getInternalColumn(i,j,&v_ij); CHKERRQ(ierr);

      for (int k = 0; k < grid.Mz; ++k) {
        u_ij[k] = U[i][j][1].u;
        v_ij[k] = U[i][j][1].v;
      }
    }
  }

  ierr = v.end_access(); CHKERRQ(ierr);
  ierr = u.end_access(); CHKERRQ(ierr);

  ierr = DMDAVecRestoreArray(dmmg[0]->dm, dmmg[0]->x, &U); CHKERRQ(ierr);

  ierr = u.update_ghosts(); CHKERRQ(ierr);
  ierr = v.update_ghosts(); CHKERRQ(ierr);

  return 0;
}


PetscErrorCode BlatterStressBalance::get_max_2d_velocity(PetscReal &/*maxu*/, PetscReal &/*maxv*/) {
  // compute in mesh_to_regular_grid()
  return 0;
}

PetscErrorCode BlatterStressBalance::get_3d_velocity(IceModelVec3* &u_out, IceModelVec3* &v_out, IceModelVec3* &w_out) {
  u_out = &u;
  v_out = &v;
  w_out = &w;
  return 0;
}

PetscErrorCode BlatterStressBalance::get_max_3d_velocity(PetscReal &maxu, PetscReal &maxv, PetscReal &maxw) {
  // FIXME: this is probably not right *NOR* efficient.
  // revise for correctness by looking at FEM nodal values of horizontal velocity
  //   grid values on IceModel grid for vertical velocity
  // revise for efficiency if: finding this max could be done earlier or as part of update() or
  //   with less communication or better memory locality

  PetscErrorCode ierr;

  ierr = u.begin_access(); CHKERRQ(ierr);
  ierr = v.begin_access(); CHKERRQ(ierr);
  ierr = w.begin_access(); CHKERRQ(ierr);
  PetscReal my_umax = 0, my_vmax = 0, my_wmax = 0;
  PetscReal *ucol,*vcol,*wcol;
  for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
      ierr = u.getInternalColumn(i, j, &ucol); CHKERRQ(ierr);
      ierr = v.getInternalColumn(i, j, &vcol); CHKERRQ(ierr);
      ierr = w.getInternalColumn(i, j, &wcol); CHKERRQ(ierr);
      for (PetscInt k = 0; k < grid.Mz; ++k) {
        my_umax = PetscMax(my_umax, PetscAbs(ucol[k]));
        my_vmax = PetscMax(my_vmax, PetscAbs(vcol[k]));
        my_wmax = PetscMax(my_wmax, PetscAbs(wcol[k]));
      }
    }
  }
  ierr = w.end_access(); CHKERRQ(ierr);
  ierr = v.end_access(); CHKERRQ(ierr);
  ierr = u.end_access(); CHKERRQ(ierr);

  ierr = PISMGlobalMax(&my_umax, &maxu, grid.com); CHKERRQ(ierr);
  ierr = PISMGlobalMax(&my_vmax, &maxv, grid.com); CHKERRQ(ierr);
  ierr = PISMGlobalMax(&my_wmax, &maxw, grid.com); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode BlatterStressBalance::extend_the_grid(PetscInt old_Mz) {
  PetscErrorCode ierr;
  ierr = u.extend_vertically(old_Mz, 0.0); CHKERRQ(ierr);
  ierr = v.extend_vertically(old_Mz, 0.0); CHKERRQ(ierr);
  ierr = w.extend_vertically(old_Mz, 0.0); CHKERRQ(ierr);
  ierr = Sigma.extend_vertically(old_Mz, 0.0); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode BlatterStressBalance::compute_basal_frictional_heating() {
  SETERRQ(grid.com, 1, "not implemented");
  return 0;
}

PetscErrorCode BlatterStressBalance::compute_volumetric_strain_heating() {
  SETERRQ(grid.com, 1, "not implemented");
  return 0;
}
