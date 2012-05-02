// Copyright (C) 2010, 2011, 2012 Constantine Khroulev
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

#include "PISMStressBalance_diagnostics.hh"
#include "Mask.hh"
#include "ShallowStressBalance.hh"
#include "SSB_Modifier.hh"
#include "PISMVars.hh"

void PISMStressBalance::get_diagnostics(map<string, PISMDiagnostic*> &dict) {

  dict["bfrict"] = new PSB_bfrict(this, grid, *variables);

  dict["cbar"]     = new PSB_cbar(this,     grid, *variables);
  dict["cflx"]     = new PSB_cflx(this,     grid, *variables);
  dict["cbase"]    = new PSB_cbase(this,    grid, *variables);
  dict["csurf"]    = new PSB_csurf(this,    grid, *variables);

  dict["uvel"]     = new PSB_uvel(this, grid, *variables);
  dict["vvel"]     = new PSB_vvel(this, grid, *variables);

  dict["strainheat"] = new PSB_strainheat(this, grid, *variables);

  dict["velbar"]   = new PSB_velbar(this,   grid, *variables);
  dict["velbase"]  = new PSB_velbase(this,  grid, *variables);
  dict["velsurf"]  = new PSB_velsurf(this,  grid, *variables);

  dict["wvel"]     = new PSB_wvel(this,     grid, *variables);
  dict["wvelbase"] = new PSB_wvelbase(this, grid, *variables);
  dict["wvelsurf"] = new PSB_wvelsurf(this, grid, *variables);
  dict["wvel_rel"] = new PSB_wvel_rel(this, grid, *variables);
  dict["taud_mag"] = new PSB_taud_mag(this, grid, *variables);


  stress_balance->get_diagnostics(dict);
  modifier->get_diagnostics(dict);
}

PSB_velbar::PSB_velbar(PISMStressBalance *m, IceGrid &g, PISMVars &my_vars)
  : PISMDiag<PISMStressBalance>(m, g, my_vars) {

  dof = 2;
  vars.resize(dof);

  // set metadata:
  vars[0].init_2d("ubar", grid);
  vars[1].init_2d("vbar", grid);

  set_attrs("vertical mean of horizontal ice velocity in the X direction",
            "land_ice_vertical_mean_x_velocity",
            "m s-1", "m year-1", 0);
  set_attrs("vertical mean of horizontal ice velocity in the Y direction",
            "land_ice_vertical_mean_y_velocity",
            "m s-1", "m year-1", 1);
}

PetscErrorCode PSB_velbar::compute(IceModelVec* &output) {
  PetscErrorCode ierr;

  IceModelVec3 *u3, *v3, *w3;
  IceModelVec2S *thickness;
  IceModelVec2V *result;
  PetscScalar *u_ij, *v_ij;

  result = new IceModelVec2V;
  ierr = result->create(grid, "bar", false); CHKERRQ(ierr);
  ierr = result->set_metadata(vars[0], 0); CHKERRQ(ierr);
  ierr = result->set_metadata(vars[1], 1); CHKERRQ(ierr);

  thickness = dynamic_cast<IceModelVec2S*>(variables.get("land_ice_thickness"));
  if (thickness == NULL) SETERRQ(grid.com, 1, "land_ice_thickness is not available");

  ierr = model->get_3d_velocity(u3, v3, w3); CHKERRQ(ierr);

  ierr = u3->begin_access(); CHKERRQ(ierr);
  ierr = v3->begin_access(); CHKERRQ(ierr);
  ierr = thickness->begin_access(); CHKERRQ(ierr);
  ierr = result->begin_access(); CHKERRQ(ierr);

  for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
      PetscReal u_sum = 0, v_sum = 0,
        thk = (*thickness)(i,j);
      PetscInt ks = grid.kBelowHeight(thk);

      // an "ice-free" cell:
      if (thk < 0.1) {
        (*result)(i,j).u = 0;
        (*result)(i,j).v = 0;
        continue;
      }

      // an ice-filled cell:
      ierr = u3->getInternalColumn(i, j, &u_ij); CHKERRQ(ierr);
      ierr = v3->getInternalColumn(i, j, &v_ij); CHKERRQ(ierr);

      if (thk <= grid.zlevels[1]) {
        (*result)(i,j).u = u_ij[0];
        (*result)(i,j).v = v_ij[0];
        continue;
      }

      for (int k = 1; k <= ks; ++k) {
        u_sum += (grid.zlevels[k] - grid.zlevels[k-1]) * (u_ij[k] + u_ij[k-1]);
        v_sum += (grid.zlevels[k] - grid.zlevels[k-1]) * (v_ij[k] + v_ij[k-1]);
      }

      // Finish the trapezoidal rule integration (times 1/2) and turn this
      // integral into a vertical average. Note that we ignore the ice between
      // zlevels[ks] and the surface, so in order to have a true average we
      // divide by zlevels[ks] and not thk.
      (*result)(i,j).u = 0.5 * u_sum / grid.zlevels[ks];
      (*result)(i,j).v = 0.5 * v_sum / grid.zlevels[ks];
    }
  }

  ierr = result->end_access(); CHKERRQ(ierr);
  ierr = thickness->end_access(); CHKERRQ(ierr);
  ierr = v3->end_access(); CHKERRQ(ierr);
  ierr = u3->end_access(); CHKERRQ(ierr);

  output = result;
  return 0;
}

PSB_cbar::PSB_cbar(PISMStressBalance *m, IceGrid &g, PISMVars &my_vars)
  : PISMDiag<PISMStressBalance>(m, g, my_vars) {

  // set metadata:
  vars[0].init_2d("cbar", grid);

  set_attrs("magnitude of vertically-integrated horizontal velocity of ice", "",
            "m s-1", "m year-1", 0);
  vars[0].set("_FillValue", convert(-0.01, "m/year", "m/second"));
  vars[0].set("valid_min", 0.0);
}

PetscErrorCode PSB_cbar::compute(IceModelVec* &output) {
  PetscErrorCode ierr;
  IceModelVec *tmp;
  IceModelVec2V *velbar_vec;
  IceModelVec2S *thickness, *result;

  // get the thickness
  thickness = dynamic_cast<IceModelVec2S*>(variables.get("land_ice_thickness"));
  if (thickness == NULL) SETERRQ(grid.com, 1, "land_ice_thickness is not available");

  result = new IceModelVec2S;
  ierr = result->create(grid, "cbar", false);
  ierr = result->set_metadata(vars[0], 0); CHKERRQ(ierr);

  // compute vertically-averaged horizontal velocity:
  PSB_velbar velbar(model, grid, variables);
  ierr = velbar.compute(tmp); CHKERRQ(ierr);

  velbar_vec = dynamic_cast<IceModelVec2V*>(tmp);
  if (velbar_vec == NULL) SETERRQ(grid.com, 1, "dynamic cast failure");

  // compute its magnitude:
  ierr = velbar_vec->magnitude(*result); CHKERRQ(ierr);

  // mask out ice-free areas:
  ierr = result->mask_by(*thickness, convert(-0.01, "m/year", "m/second")); CHKERRQ(ierr);

  delete tmp;
  output = result;
  return 0;
}

PSB_cflx::PSB_cflx(PISMStressBalance *m, IceGrid &g, PISMVars &my_vars)
  : PISMDiag<PISMStressBalance>(m, g, my_vars) {

  // set metadata:
  vars[0].init_2d("cflx", grid);

  set_attrs("magnitude of vertically-integrated horizontal flux of ice", "",
            "m2 s-1", "m2 year-1", 0);
  vars[0].set("_FillValue", convert(-0.01, "m2/year", "m2/second"));
  vars[0].set("valid_min", 0.0);
}

PetscErrorCode PSB_cflx::compute(IceModelVec* &output) {
  PetscErrorCode ierr;
  IceModelVec2S *thickness, *result;
  IceModelVec *tmp;

  // get the thickness
  thickness = dynamic_cast<IceModelVec2S*>(variables.get("land_ice_thickness"));
  if (thickness == NULL) SETERRQ(grid.com, 1, "land_ice_thickness is not available");

  // Compute the vertically-average horizontal ice velocity:
  PSB_cbar cbar(model, grid, variables);
  ierr = cbar.compute(tmp); CHKERRQ(ierr);
  // NB: the call above allocates memory

  result = dynamic_cast<IceModelVec2S*>(tmp);
  if (result == NULL) SETERRQ(grid.com, 1, "dynamic_cast failure");

  ierr = thickness->begin_access(); CHKERRQ(ierr);
  ierr = result->begin_access(); CHKERRQ(ierr);

  for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i)
    for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j)
      (*result)(i,j) *= (*thickness)(i,j);

  ierr = result->end_access(); CHKERRQ(ierr);
  ierr = thickness->end_access(); CHKERRQ(ierr);

  ierr = result->mask_by(*thickness, convert(-0.01, "m2/year", "m2/second")); CHKERRQ(ierr);

  ierr = result->set_metadata(vars[0], 0); CHKERRQ(ierr);

  output = result;
  return 0;
}

PSB_cbase::PSB_cbase(PISMStressBalance *m, IceGrid &g, PISMVars &my_vars)
 : PISMDiag<PISMStressBalance>(m, g, my_vars) {

  // set metadata:
  vars[0].init_2d("cbase", grid);

  set_attrs("magnitude of horizontal velocity of ice at base of ice", "",
            "m s-1", "m year-1", 0);
  vars[0].set("_FillValue", convert(-0.01, "m/year", "m/second"));
  vars[0].set("valid_min", 0.0);
}

PetscErrorCode PSB_cbase::compute(IceModelVec* &output) {
  PetscErrorCode ierr;
  PetscScalar fill_value = convert(-0.01, "m/year", "m/second");
  IceModelVec3 *u3, *v3, *w3;
  IceModelVec2S tmp, *result, *thickness;

  ierr = tmp.create(grid, "tmp", false); CHKERRQ(ierr);

  result = new IceModelVec2S;
  ierr = result->create(grid, "cbase", false); CHKERRQ(ierr);
  ierr = result->set_metadata(vars[0], 0); CHKERRQ(ierr);

  ierr = model->get_3d_velocity(u3, v3, w3); CHKERRQ(ierr);

  thickness = dynamic_cast<IceModelVec2S*>(variables.get("land_ice_thickness"));
  if (thickness == NULL) SETERRQ(grid.com, 1, "land_ice_thickness is not available");

  ierr = u3->getHorSlice(*result, 0.0); CHKERRQ(ierr); // result = u_{z=0}
  ierr = v3->getHorSlice(tmp, 0.0); CHKERRQ(ierr);    // tmp = v_{z=0}

  ierr = result->set_to_magnitude(*result, tmp); CHKERRQ(ierr);

  ierr = result->mask_by(*thickness, fill_value); CHKERRQ(ierr); // mask out ice-free areas

  output = result;
  return 0;
}

PSB_csurf::PSB_csurf(PISMStressBalance *m, IceGrid &g, PISMVars &my_vars)
  : PISMDiag<PISMStressBalance>(m, g, my_vars) {
  PetscReal fill_value = convert(-0.01, "m/year", "m/second");
  // set metadata:
  vars[0].init_2d("csurf", grid);

  set_attrs("magnitude of horizontal velocity of ice at ice surface", "",
            "m s-1", "m year-1", 0);
  vars[0].set("_FillValue", fill_value);
  vars[0].set("valid_min",  0.0);
}

PetscErrorCode PSB_csurf::compute(IceModelVec* &output) {
  PetscErrorCode ierr;
  PetscReal fill_value = convert(-0.01, "m/year", "m/second");

  IceModelVec3 *u3, *v3, *w3;
  IceModelVec2S tmp, *result, *thickness;

  ierr = tmp.create(grid, "tmp", false); CHKERRQ(ierr);

  result = new IceModelVec2S;
  ierr = result->create(grid, "csurf", false); CHKERRQ(ierr);
  ierr = result->set_metadata(vars[0], 0); CHKERRQ(ierr);

  ierr = model->get_3d_velocity(u3, v3, w3); CHKERRQ(ierr);

  thickness = dynamic_cast<IceModelVec2S*>(variables.get("land_ice_thickness"));
  if (thickness == NULL) SETERRQ(grid.com, 1, "land_ice_thickness is not available");

  ierr = u3->getSurfaceValues(*result, *thickness); CHKERRQ(ierr);
  ierr = v3->getSurfaceValues(tmp, *thickness); CHKERRQ(ierr);

  ierr = result->set_to_magnitude(*result, tmp); CHKERRQ(ierr);

  ierr = result->mask_by(*thickness, fill_value); CHKERRQ(ierr); // mask out ice-free areas

  output = result;
  return 0;
}


PSB_velsurf::PSB_velsurf(PISMStressBalance *m, IceGrid &g, PISMVars &my_vars)
  : PISMDiag<PISMStressBalance>(m, g, my_vars) {

  // set metadata:
  dof = 2;
  vars.resize(dof);

  vars[0].init_2d("uvelsurf", grid);
  vars[1].init_2d("vvelsurf", grid);

  set_attrs("x-component of the horizontal velocity of ice at ice surface", "",
            "m s-1", "m year-1", 0);
  set_attrs("y-component of the horizontal velocity of ice at ice surface", "",
            "m s-1", "m year-1", 1);

  vars[0].set("valid_min", convert(-1e6, "m/year", "m/second"));
  vars[0].set("valid_max", convert(1e6, "m/year", "m/second"));
  vars[0].set("_FillValue", convert(2e6, "m/year", "m/second"));

  vars[1].set("valid_min", convert(-1e6, "m/year", "m/second"));
  vars[1].set("valid_max", convert(1e6, "m/year", "m/second"));
  vars[1].set("_FillValue", convert(2e6, "m/year", "m/second"));
}

PetscErrorCode PSB_velsurf::compute(IceModelVec* &output) {
  PetscErrorCode ierr;
  IceModelVec2V *result;
  IceModelVec3 *u3, *v3, *w3;
  IceModelVec2S *thickness, tmp;
  PetscScalar fill_value = convert(2e6, "m/year", "m/second");

  result = new IceModelVec2V;
  ierr = result->create(grid, "surf", false); CHKERRQ(ierr);
  ierr = result->set_metadata(vars[0], 0); CHKERRQ(ierr);
  ierr = result->set_metadata(vars[1], 1); CHKERRQ(ierr);

  ierr = tmp.create(grid, "tmp", false); CHKERRQ(ierr);

  ierr = model->get_3d_velocity(u3, v3, w3); CHKERRQ(ierr);

  thickness = dynamic_cast<IceModelVec2S*>(variables.get("land_ice_thickness"));
  if (thickness == NULL) SETERRQ(grid.com, 1, "land_ice_thickness is not available");

  ierr = u3->getSurfaceValues(tmp, *thickness); CHKERRQ(ierr);
  ierr = result->set_component(0, tmp); CHKERRQ(ierr);

  ierr = v3->getSurfaceValues(tmp, *thickness); CHKERRQ(ierr);
  ierr = result->set_component(1, tmp); CHKERRQ(ierr);

  IceModelVec2Int *mask = dynamic_cast<IceModelVec2Int*>(variables.get("mask"));
  if (mask == NULL) SETERRQ(grid.com, 1, "mask is not available");

  MaskQuery M(*mask);
  ierr = mask->begin_access(); CHKERRQ(ierr);
  ierr = result->begin_access(); CHKERRQ(ierr);

  for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
      if (M.ice_free(i, j)) {
        (*result)(i, j).u = fill_value;
        (*result)(i, j).v = fill_value;
      }
    }
  }

  ierr = result->end_access(); CHKERRQ(ierr);
  ierr = mask->end_access(); CHKERRQ(ierr);

  output = result;
  return 0;
}

PSB_wvel::PSB_wvel(PISMStressBalance *m, IceGrid &g, PISMVars &my_vars)
  : PISMDiag<PISMStressBalance>(m, g, my_vars) {

  // set metadata:
  vars[0].init_3d("wvel", grid, g.zlevels);

  set_attrs("vertical velocity of ice, relative to geoid", "",
            "m s-1", "m year-1", 0);
  vars[0].set("valid_min", convert(-1e6, "m/year", "m/second"));
  vars[0].set("valid_max", convert(1e6, "m/year", "m/second"));
}

PetscErrorCode PSB_wvel::compute(IceModelVec* &output) {
  PetscErrorCode ierr;
  IceModelVec3 *result, *u3, *v3, *w3;
  IceModelVec2S *bed, *uplift;
  PetscScalar *u, *v, *w, *res;

  result = new IceModelVec3;
  ierr = result->create(grid, "wvel", false); CHKERRQ(ierr);
  ierr = result->set_metadata(vars[0], 0); CHKERRQ(ierr);

  bed = dynamic_cast<IceModelVec2S*>(variables.get("bedrock_altitude"));
  if (bed == NULL) SETERRQ(grid.com, 1, "bedrock_altitude is not available");

  uplift = dynamic_cast<IceModelVec2S*>(variables.get("tendency_of_bedrock_altitude"));
  if (uplift == NULL) SETERRQ(grid.com, 1, "tendency_of_bedrock_altitude is not available");

  ierr = model->get_3d_velocity(u3, v3, w3); CHKERRQ(ierr);

  ierr = bed->begin_access(); CHKERRQ(ierr);
  ierr = u3->begin_access(); CHKERRQ(ierr);
  ierr = v3->begin_access(); CHKERRQ(ierr);
  ierr = w3->begin_access(); CHKERRQ(ierr);
  ierr = uplift->begin_access(); CHKERRQ(ierr);
  ierr = result->begin_access(); CHKERRQ(ierr);

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      ierr = u3->getInternalColumn(i, j, &u); CHKERRQ(ierr);
      ierr = v3->getInternalColumn(i, j, &v); CHKERRQ(ierr);
      ierr = w3->getInternalColumn(i, j, &w); CHKERRQ(ierr);
      ierr = result->getInternalColumn(i, j, &res); CHKERRQ(ierr);

      for (PetscInt k = 0; k < grid.Mz; ++k)
	res[k] = w[k] + (*uplift)(i,j) + u[k] * bed->diff_x_p(i,j) + v[k] * bed->diff_y_p(i,j);
    }
  }

  ierr = result->end_access(); CHKERRQ(ierr);
  ierr = uplift->end_access(); CHKERRQ(ierr);
  ierr = w3->end_access(); CHKERRQ(ierr);
  ierr = v3->end_access(); CHKERRQ(ierr);
  ierr = u3->end_access(); CHKERRQ(ierr);
  ierr = bed->end_access(); CHKERRQ(ierr);

  output = result;
  return 0;
}

PSB_wvelsurf::PSB_wvelsurf(PISMStressBalance *m, IceGrid &g, PISMVars &my_vars)
  : PISMDiag<PISMStressBalance>(m, g, my_vars) {

  // set metadata:
  vars[0].init_2d("wvelsurf", grid);

  set_attrs("vertical velocity of ice at ice surface, relative to the geoid", "",
            "m s-1", "m year-1", 0);
  vars[0].set("valid_min", convert(-1e6, "m/year", "m/second"));
  vars[0].set("valid_max", convert(1e6, "m/year", "m/second"));
  vars[0].set("_FillValue", convert(2e6, "m/year", "m/second"));
}

PetscErrorCode PSB_wvelsurf::compute(IceModelVec* &output) {
  PetscErrorCode ierr;
  IceModelVec *tmp;
  IceModelVec3 *w3;
  IceModelVec2S *result, *thickness;
  PetscScalar fill_value = convert(2e6, "m/year", "m/second");

  result = new IceModelVec2S;
  ierr = result->create(grid, "wvelsurf", false); CHKERRQ(ierr);
  ierr = result->set_metadata(vars[0], 0); CHKERRQ(ierr);

  PSB_wvel wvel(model, grid, variables);

  ierr = wvel.compute(tmp); CHKERRQ(ierr);

  w3 = dynamic_cast<IceModelVec3*>(tmp);
  if (tmp == NULL) SETERRQ(grid.com, 1, "dynamic_cast failure");

  thickness = dynamic_cast<IceModelVec2S*>(variables.get("land_ice_thickness"));
  if (thickness == NULL) SETERRQ(grid.com, 1, "land_ice_thickness is not available");

  ierr = w3->getSurfaceValues(*result, *thickness); CHKERRQ(ierr);

  IceModelVec2Int *mask = dynamic_cast<IceModelVec2Int*>(variables.get("mask"));
  if (mask == NULL) SETERRQ(grid.com, 1, "mask is not available");

  MaskQuery M(*mask);
  ierr = mask->begin_access(); CHKERRQ(ierr);
  ierr = result->begin_access(); CHKERRQ(ierr);

  for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
      if (M.ice_free(i, j))
        (*result)(i, j) = fill_value;
    }
  }

  ierr = result->end_access(); CHKERRQ(ierr);
  ierr = mask->end_access(); CHKERRQ(ierr);

  delete tmp;
  output = result;
  return 0;
}

PSB_wvelbase::PSB_wvelbase(PISMStressBalance *m, IceGrid &g, PISMVars &my_vars)
  : PISMDiag<PISMStressBalance>(m, g, my_vars) {

  // set metadata:
  vars[0].init_2d("wvelbase", grid);

  set_attrs("vertical velocity of ice at the base of ice, relative to the geoid", "",
            "m s-1", "m year-1", 0);
  vars[0].set("valid_min", convert(-1e6, "m/year", "m/second"));
  vars[0].set("valid_max", convert(1e6, "m/year", "m/second"));
  vars[0].set("_FillValue", convert(2e6, "m/year", "m/second"));
}

PetscErrorCode PSB_wvelbase::compute(IceModelVec* &output) {
  PetscErrorCode ierr;
  IceModelVec *tmp;
  IceModelVec3 *w3;
  IceModelVec2S *result;
  PetscScalar fill_value = convert(2e6, "m/year", "m/second");

  result = new IceModelVec2S;
  ierr = result->create(grid, "wvelbase", false); CHKERRQ(ierr);
  ierr = result->set_metadata(vars[0], 0); CHKERRQ(ierr);

  PSB_wvel wvel(model, grid, variables);

  ierr = wvel.compute(tmp); CHKERRQ(ierr);

  w3 = dynamic_cast<IceModelVec3*>(tmp);
  if (tmp == NULL) SETERRQ(grid.com, 1, "dynamic_cast failure");

  ierr = w3->getHorSlice(*result, 0.0); CHKERRQ(ierr);

  IceModelVec2Int *mask = dynamic_cast<IceModelVec2Int*>(variables.get("mask"));
  if (mask == NULL) SETERRQ(grid.com, 1, "mask is not available");

  MaskQuery M(*mask);
  ierr = mask->begin_access(); CHKERRQ(ierr);
  ierr = result->begin_access(); CHKERRQ(ierr);

  for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
      if (M.ice_free(i, j))
        (*result)(i, j) = fill_value;
    }
  }

  ierr = result->end_access(); CHKERRQ(ierr);
  ierr = mask->end_access(); CHKERRQ(ierr);

  delete tmp;
  output = result;
  return 0;
}

PSB_velbase::PSB_velbase(PISMStressBalance *m, IceGrid &g, PISMVars &my_vars)
  : PISMDiag<PISMStressBalance>(m, g, my_vars) {

  // set metadata:
  dof = 2;
  vars.resize(dof);

  vars[0].init_2d("uvelbase", grid);
  vars[1].init_2d("vvelbase", grid);

  set_attrs("x-component of the horizontal velocity of ice at the base of ice", "",
            "m s-1", "m year-1", 0);
  set_attrs("y-component of the horizontal velocity of ice at the base of ice", "",
            "m s-1", "m year-1", 1);

  vars[0].set("valid_min", convert(-1e6, "m/year", "m/second"));
  vars[0].set("valid_max", convert(1e6, "m/year", "m/second"));
  vars[0].set("_FillValue", convert(2e6, "m/year", "m/second"));

  vars[1].set("valid_min", convert(-1e6, "m/year", "m/second"));
  vars[1].set("valid_max", convert(1e6, "m/year", "m/second"));
  vars[1].set("_FillValue", convert(2e6, "m/year", "m/second"));
}

PetscErrorCode PSB_velbase::compute(IceModelVec* &output) {
  PetscErrorCode ierr;
  IceModelVec2V *result;
  IceModelVec3 *u3, *v3, *w3;
  IceModelVec2S tmp;            // will be de-allocated automatically
  PetscScalar fill_value = convert(2e6, "m/year", "m/second");

  result = new IceModelVec2V;
  ierr = result->create(grid, "base", false); CHKERRQ(ierr);
  ierr = result->set_metadata(vars[0], 0); CHKERRQ(ierr);
  ierr = result->set_metadata(vars[1], 1); CHKERRQ(ierr);

  ierr = tmp.create(grid, "tmp", false); CHKERRQ(ierr);

  ierr = model->get_3d_velocity(u3, v3, w3); CHKERRQ(ierr);

  ierr = u3->getHorSlice(tmp, 0.0); CHKERRQ(ierr);
  ierr = result->set_component(0, tmp); CHKERRQ(ierr);

  ierr = v3->getHorSlice(tmp, 0.0); CHKERRQ(ierr);
  ierr = result->set_component(1, tmp); CHKERRQ(ierr);

  IceModelVec2Int *mask = dynamic_cast<IceModelVec2Int*>(variables.get("mask"));
  if (mask == NULL) SETERRQ(grid.com, 1, "mask is not available");

  MaskQuery M(*mask);
  ierr = mask->begin_access(); CHKERRQ(ierr);
  ierr = result->begin_access(); CHKERRQ(ierr);

  for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
      if (M.ice_free(i, j)) {
        (*result)(i, j).u = fill_value;
        (*result)(i, j).v = fill_value;
      }
    }
  }

  ierr = result->end_access(); CHKERRQ(ierr);
  ierr = mask->end_access(); CHKERRQ(ierr);

  output = result;
  return 0;
}


PSB_bfrict::PSB_bfrict(PISMStressBalance *m, IceGrid &g, PISMVars &my_vars)
  : PISMDiag<PISMStressBalance>(m, g, my_vars) {

  // set metadata:
  vars[0].init_2d("bfrict", grid);

  set_attrs("basal frictional heating", "",
            "W m-2", "W m-2", 0);
}

PetscErrorCode PSB_bfrict::compute(IceModelVec* &output) {
  PetscErrorCode ierr;

  IceModelVec2S *result = new IceModelVec2S;
  ierr = result->create(grid, "bfrict", false); CHKERRQ(ierr);
  ierr = result->set_metadata(vars[0], 0); CHKERRQ(ierr);

  IceModelVec2S *bfrict;
  ierr = model->get_basal_frictional_heating(bfrict); CHKERRQ(ierr);

  ierr = bfrict->copy_to(*result); CHKERRQ(ierr);

  output = result;
  return 0;
}


PSB_uvel::PSB_uvel(PISMStressBalance *m, IceGrid &g, PISMVars &my_vars)
  : PISMDiag<PISMStressBalance>(m, g, my_vars) {

  // set metadata:
  vars[0].init_3d("uvel", grid, g.zlevels);

  set_attrs("horizontal velocity of ice in the X direction", "land_ice_x_velocity",
            "m s-1", "m year-1", 0);
}

PetscErrorCode PSB_uvel::compute(IceModelVec* &output) {
  PetscErrorCode ierr;

  IceModelVec3 *result = new IceModelVec3;
  ierr = result->create(grid, "uvel", false); CHKERRQ(ierr);
  ierr = result->set_metadata(vars[0], 0); CHKERRQ(ierr);

  IceModelVec2S *thickness;
  thickness = dynamic_cast<IceModelVec2S*>(variables.get("land_ice_thickness"));
  if (thickness == NULL) SETERRQ(grid.com, 1, "land_ice_thickness is not available");

  IceModelVec3 *u3, *v3, *w3;
  ierr = model->get_3d_velocity(u3, v3, w3); CHKERRQ(ierr);

  ierr = u3->begin_access(); CHKERRQ(ierr);
  ierr = result->begin_access(); CHKERRQ(ierr);
  ierr = thickness->begin_access(); CHKERRQ(ierr);

  PetscScalar *u_ij, *u_out_ij;
  for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
      int ks = grid.kBelowHeight((*thickness)(i,j));

      ierr = u3->getInternalColumn(i,j,&u_ij); CHKERRQ(ierr);
      ierr = result->getInternalColumn(i,j,&u_out_ij); CHKERRQ(ierr);

      // in the ice:
      for (int k = 0; k <= ks ; k++) {
        u_out_ij[k] = u_ij[k];
      }
      // above the ice:
      for (int k = ks+1; k < grid.Mz ; k++) {
        u_out_ij[k] = 0.0;
      }
    }
  }

  ierr = thickness->end_access(); CHKERRQ(ierr);
  ierr = result->end_access(); CHKERRQ(ierr);
  ierr = u3->end_access(); CHKERRQ(ierr);

  output = result;
  return 0;
}

PSB_vvel::PSB_vvel(PISMStressBalance *m, IceGrid &g, PISMVars &my_vars)
  : PISMDiag<PISMStressBalance>(m, g, my_vars) {

  // set metadata:
  vars[0].init_3d("vvel", grid, g.zlevels);

  set_attrs("horizontal velocity of ice in the Y direction", "land_ice_y_velocity",
            "m s-1", "m year-1", 0);
}

PetscErrorCode PSB_vvel::compute(IceModelVec* &output) {
  PetscErrorCode ierr;

  IceModelVec3 *result = new IceModelVec3;
  ierr = result->create(grid, "vvel", false); CHKERRQ(ierr);
  ierr = result->set_metadata(vars[0], 0); CHKERRQ(ierr);

  IceModelVec2S *thickness;
  thickness = dynamic_cast<IceModelVec2S*>(variables.get("land_ice_thickness"));
  if (thickness == NULL) SETERRQ(grid.com, 1, "land_ice_thickness is not available");

  IceModelVec3 *u3, *v3, *w3;
  ierr = model->get_3d_velocity(u3, v3, w3); CHKERRQ(ierr);

  ierr = v3->begin_access(); CHKERRQ(ierr);
  ierr = result->begin_access(); CHKERRQ(ierr);
  ierr = thickness->begin_access(); CHKERRQ(ierr);

  PetscScalar *v_ij, *v_out_ij;
  for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
      int ks = grid.kBelowHeight((*thickness)(i,j));

      ierr = v3->getInternalColumn(i,j,&v_ij); CHKERRQ(ierr);
      ierr = result->getInternalColumn(i,j,&v_out_ij); CHKERRQ(ierr);

      // in the ice:
      for (int k = 0; k <= ks ; k++) {
        v_out_ij[k] = v_ij[k];
      }
      // above the ice:
      for (int k = ks+1; k < grid.Mz ; k++) {
        v_out_ij[k] = 0.0;
      }
    }
  }

  ierr = thickness->end_access(); CHKERRQ(ierr);
  ierr = result->end_access(); CHKERRQ(ierr);
  ierr = v3->end_access(); CHKERRQ(ierr);

  output = result;
  return 0;
}

PSB_wvel_rel::PSB_wvel_rel(PISMStressBalance *m, IceGrid &g, PISMVars &my_vars)
  : PISMDiag<PISMStressBalance>(m, g, my_vars) {

  // set metadata:
  vars[0].init_3d("wvel_rel", grid, g.zlevels);

  set_attrs("vertical velocity of ice, relative to base of ice directly below", "",
            "m s-1", "m year-1", 0);
}

PetscErrorCode PSB_wvel_rel::compute(IceModelVec* &output) {
  PetscErrorCode ierr;

  IceModelVec3 *result = new IceModelVec3;
  ierr = result->create(grid, "wvel_rel", false); CHKERRQ(ierr);
  ierr = result->set_metadata(vars[0], 0); CHKERRQ(ierr);

  IceModelVec2S *thickness;
  thickness = dynamic_cast<IceModelVec2S*>(variables.get("land_ice_thickness"));
  if (thickness == NULL) SETERRQ(grid.com, 1, "land_ice_thickness is not available");

  IceModelVec3 *u3, *v3, *w3;
  ierr = model->get_3d_velocity(u3, v3, w3); CHKERRQ(ierr);

  ierr = w3->begin_access(); CHKERRQ(ierr);
  ierr = result->begin_access(); CHKERRQ(ierr);
  ierr = thickness->begin_access(); CHKERRQ(ierr);

  PetscScalar *w_ij, *w_out_ij;
  for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
      int ks = grid.kBelowHeight((*thickness)(i,j));

      ierr = w3->getInternalColumn(i,j,&w_ij); CHKERRQ(ierr);
      ierr = result->getInternalColumn(i,j,&w_out_ij); CHKERRQ(ierr);

      // in the ice:
      for (int k = 0; k <= ks ; k++) {
        w_out_ij[k] = w_ij[k];
      }
      // above the ice:
      for (int k = ks+1; k < grid.Mz ; k++) {
        w_out_ij[k] = 0.0;
      }
    }
  }

  ierr = thickness->end_access(); CHKERRQ(ierr);
  ierr = result->end_access(); CHKERRQ(ierr);
  ierr = w3->end_access(); CHKERRQ(ierr);

  output = result;
  return 0;
}

PSB_taud_mag::PSB_taud_mag(PISMStressBalance *m, IceGrid &g, PISMVars &my_vars)
  : PISMDiag<PISMStressBalance>(m, g, my_vars) {

  // set metadata:
  vars[0].init_2d("taud_mag", grid);

  set_attrs("magnitude of the driving shear stress at the base of ice", "",
            "Pa", "Pa", 0);
}

PetscErrorCode PSB_taud_mag::compute(IceModelVec* &output) {
  PetscErrorCode ierr;

  // Allocate memory:
  IceModelVec2S *result = new IceModelVec2S;
  ierr = result->create(grid, "taud_mag", false); CHKERRQ(ierr);
  ierr = result->set_metadata(vars[0], 0); CHKERRQ(ierr);
  result->write_in_glaciological_units = true;

  IceModelVec2S *thickness, *surface, *bed;
  IceModelVec2Int *mask;

  thickness = dynamic_cast<IceModelVec2S*>(variables.get("land_ice_thickness"));
  if (thickness == NULL) SETERRQ(grid.com, 1, "land_ice_thickness is not available");

  surface = dynamic_cast<IceModelVec2S*>(variables.get("surface_altitude"));
  if (surface == NULL) SETERRQ(grid.com, 1, "surface_altitude is not available");

  bed = dynamic_cast<IceModelVec2S*>(variables.get("bedrock_altitude"));
  if (bed == NULL) SETERRQ(grid.com, 1, "bedrock_altitude is not available");

  mask = dynamic_cast<IceModelVec2Int*>(variables.get("mask"));
  if (mask == NULL) SETERRQ(grid.com, 1, "mask is not available");

  IceModelVec2S &thk = *thickness; // to improve readability (below)

  const PetscScalar n       = model->config.get("Glen_exponent"), // frequently n = 3
                    etapow  = (2.0 * n + 2.0)/n,  // = 8/3 if n = 3
                    invpow  = 1.0 / etapow,  // = 3/8
                    dinvpow = (- n - 2.0) / (2.0 * n + 2.0); // = -5/8
  const PetscScalar minThickEtaTransform = 5.0; // m
  const PetscScalar dx=grid.dx, dy=grid.dy;

  PetscReal standard_gravity = model->config.get("standard_gravity"),
    ice_density = model->config.get("ice_density");
  bool use_eta = (model->config.get_string("surface_gradient_method") == "eta");

  MaskQuery M(*mask);

  ierr =   surface->begin_access();    CHKERRQ(ierr);
  ierr =       bed->begin_access();  CHKERRQ(ierr);
  ierr =      mask->begin_access();  CHKERRQ(ierr);
  ierr =        thk.begin_access();  CHKERRQ(ierr);

  ierr = result->begin_access(); CHKERRQ(ierr);

  PetscReal result_ij_u, result_ij_v;
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      const PetscScalar pressure = ice_density * standard_gravity * thk(i,j);
      if (pressure <= 0.0) {
        result_ij_u = 0.0;
        result_ij_v = 0.0;
      } else {
        PetscScalar h_x = 0.0, h_y = 0.0;
        // FIXME: we need to handle grid periodicity correctly.
        if (M.grounded(i,j) && (use_eta == true)) {
	        // in grounded case, differentiate eta = H^{8/3} by chain rule
          if (thk(i,j) > 0.0) {
            const PetscScalar myH = (thk(i,j) < minThickEtaTransform ?
                                     minThickEtaTransform : thk(i,j));
            const PetscScalar eta = pow(myH, etapow), factor = invpow * pow(eta, dinvpow);
            h_x = factor * (pow(thk(i+1,j),etapow) - pow(thk(i-1,j),etapow)) / (2*dx);
            h_y = factor * (pow(thk(i,j+1),etapow) - pow(thk(i,j-1),etapow)) / (2*dy);
          }
          // now add bed slope to get actual h_x,h_y
          // FIXME: there is no reason to assume user's bed is periodized
          h_x += bed->diff_x(i,j);
          h_y += bed->diff_y(i,j);
        } else {  // floating or eta transformation is not used
          h_x = surface->diff_x_p(i,j);
          h_y = surface->diff_y_p(i,j);
        }

        result_ij_u = - pressure * h_x;
        result_ij_v = - pressure * h_y;
      }

      (*result)(i,j) = sqrt(PetscSqr(result_ij_u) + PetscSqr(result_ij_v));
    }
  }

  ierr =        thk.end_access(); CHKERRQ(ierr);
  ierr =       bed->end_access(); CHKERRQ(ierr);
  ierr =   surface->end_access(); CHKERRQ(ierr);
  ierr =      mask->end_access(); CHKERRQ(ierr);
  ierr =     result->end_access(); CHKERRQ(ierr);

  output = result;
  return 0;
}

PSB_strainheat::PSB_strainheat(PISMStressBalance *m, IceGrid &g, PISMVars &my_vars)
  : PISMDiag<PISMStressBalance>(m, g, my_vars) {

  // set metadata:
  vars[0].init_3d("strainheat", grid, grid.zlevels);

  set_attrs("rate of strain heating in ice (dissipation heating)", "",
            "W m-3", "mW m-3", 0);
}

PetscErrorCode PSB_strainheat::compute(IceModelVec* &output) {
  PetscErrorCode ierr;

  IceModelVec3 *result = new IceModelVec3;
  ierr = result->create(grid, "strainheat", false); CHKERRQ(ierr);
  ierr = result->set_metadata(vars[0], 0); CHKERRQ(ierr);
  result->write_in_glaciological_units = true;

  IceModelVec3 *tmp;
  ierr = model->get_volumetric_strain_heating(tmp); CHKERRQ(ierr);

  ierr = tmp->copy_to(*result); CHKERRQ(ierr);

  output = result;
  return 0;
}

