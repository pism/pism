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

  dict["bfrict"]   = new PSB_bfrict(this, grid, *variables);

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
  dict["strain_rates"] = new PSB_strain_rates(this, grid, *variables);
  dict["deviatoric_stresses"] = new PSB_deviatoric_stresses(this, grid, *variables);

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
  vars[0].set("_FillValue", convert(grid.config.get("fill_value"), "m/year", "m/s"));
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
  ierr = result->mask_by(*thickness, convert(grid.config.get("fill_value"), "m/year", "m/s")); CHKERRQ(ierr);

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
  vars[0].set("_FillValue", convert(grid.config.get("fill_value"), "m2/year", "m2/s"));
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

  ierr = result->mask_by(*thickness, convert(grid.config.get("fill_value"), "m/year", "m/s")); CHKERRQ(ierr);

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
  vars[0].set("_FillValue", convert(grid.config.get("fill_value"), "m/year", "m/s"));
  vars[0].set("valid_min", 0.0);
}

PetscErrorCode PSB_cbase::compute(IceModelVec* &output) {
  PetscErrorCode ierr;
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

  ierr = result->mask_by(*thickness, convert(grid.config.get("fill_value"), "m/year", "m/s")); CHKERRQ(ierr); // mask out ice-free areas

  output = result;
  return 0;
}

PSB_csurf::PSB_csurf(PISMStressBalance *m, IceGrid &g, PISMVars &my_vars)
  : PISMDiag<PISMStressBalance>(m, g, my_vars) {
  // set metadata:
  vars[0].init_2d("csurf", grid);

  set_attrs("magnitude of horizontal velocity of ice at ice surface", "",
            "m s-1", "m year-1", 0);
  vars[0].set("_FillValue", convert(grid.config.get("fill_value"), "m/year", "m/s"));
  vars[0].set("valid_min",  0.0);
}

PetscErrorCode PSB_csurf::compute(IceModelVec* &output) {
  PetscErrorCode ierr;

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

  ierr = result->mask_by(*thickness, convert(grid.config.get("fill_value"), "m/year", "m/s")); CHKERRQ(ierr); // mask out ice-free areas

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
  vars[0].set("_FillValue", convert(grid.config.get("fill_value"), "m/year", "m/s"));

  vars[1].set("valid_min", convert(-1e6, "m/year", "m/second"));
  vars[1].set("valid_max", convert(1e6, "m/year", "m/second"));
  vars[1].set("_FillValue", convert(grid.config.get("fill_value"), "m/year", "m/s"));
}

PetscErrorCode PSB_velsurf::compute(IceModelVec* &output) {
  PetscErrorCode ierr;
  IceModelVec2V *result;
  IceModelVec3 *u3, *v3, *w3;
  IceModelVec2S *thickness, tmp;
  PetscScalar fill_value = convert(grid.config.get("fill_value"), "m/year", "m/s");

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
  IceModelVec2S *bed, *uplift, *thickness;
  PetscScalar *u, *v, *w, *res;

  result = new IceModelVec3;
  ierr = result->create(grid, "wvel", false); CHKERRQ(ierr);
  ierr = result->set_metadata(vars[0], 0); CHKERRQ(ierr);

  bed = dynamic_cast<IceModelVec2S*>(variables.get("bedrock_altitude"));
  if (bed == NULL) SETERRQ(grid.com, 1, "bedrock_altitude is not available");

  uplift = dynamic_cast<IceModelVec2S*>(variables.get("tendency_of_bedrock_altitude"));
  if (uplift == NULL) SETERRQ(grid.com, 1, "tendency_of_bedrock_altitude is not available");

  thickness = dynamic_cast<IceModelVec2S*>(variables.get("land_ice_thickness"));
  if (thickness == NULL) SETERRQ(grid.com, 1, "land_ice_thickness is not available");

  ierr = model->get_3d_velocity(u3, v3, w3); CHKERRQ(ierr);

  ierr = thickness->begin_access(); CHKERRQ(ierr);
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

      int ks = grid.kBelowHeight((*thickness)(i,j));

      // in the ice:
      for (int k = 0; k <= ks ; k++) {
	res[k] = w[k] + (*uplift)(i,j) + u[k] * bed->diff_x_p(i,j) + v[k] * bed->diff_y_p(i,j);
      }
      // above the ice:
      for (int k = ks+1; k < grid.Mz ; k++) {
        res[k] = 0.0;
      }
    }
  }

  ierr = result->end_access(); CHKERRQ(ierr);
  ierr = uplift->end_access(); CHKERRQ(ierr);
  ierr = w3->end_access(); CHKERRQ(ierr);
  ierr = v3->end_access(); CHKERRQ(ierr);
  ierr = u3->end_access(); CHKERRQ(ierr);
  ierr = bed->end_access(); CHKERRQ(ierr);
  ierr = thickness->end_access(); CHKERRQ(ierr);

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
  vars[0].set("_FillValue", convert(grid.config.get("fill_value"), "m/year", "m/s"));
}

PetscErrorCode PSB_wvelsurf::compute(IceModelVec* &output) {
  PetscErrorCode ierr;
  IceModelVec *tmp;
  IceModelVec3 *w3;
  IceModelVec2S *result, *thickness;
  PetscScalar fill_value = convert(grid.config.get("fill_value"), "m/year", "m/s");

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
  vars[0].set("_FillValue", convert(grid.config.get("fill_value"), "m/year", "m/s"));
}

PetscErrorCode PSB_wvelbase::compute(IceModelVec* &output) {
  PetscErrorCode ierr;
  IceModelVec *tmp;
  IceModelVec3 *w3;
  IceModelVec2S *result;
  PetscScalar fill_value = convert(grid.config.get("fill_value"), "m/year", "m/s");

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
  vars[0].set("_FillValue", convert(grid.config.get("fill_value"), "m/year", "m/s"));

  vars[1].set("valid_min", convert(-1e6, "m/year", "m/second"));
  vars[1].set("valid_max", convert(1e6, "m/year", "m/second"));
  vars[1].set("_FillValue", convert(grid.config.get("fill_value"), "m/year", "m/s"));
}

PetscErrorCode PSB_velbase::compute(IceModelVec* &output) {
  PetscErrorCode ierr;
  IceModelVec2V *result;
  IceModelVec3 *u3, *v3, *w3;
  IceModelVec2S tmp;            // will be de-allocated automatically
  PetscScalar fill_value = convert(grid.config.get("fill_value"), "m/year", "m/s");

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

PSB_strain_rates::PSB_strain_rates(PISMStressBalance *m, IceGrid &g, PISMVars &my_vars)
  : PISMDiag<PISMStressBalance>(m, g, my_vars) {
  dof = 2;
  vars.resize(dof);

  // set metadata:
  vars[0].init_2d("eigen1", grid);
  vars[1].init_2d("eigen2", grid);

  set_attrs("first eigenvalue of the horizontal, vertically-integrated strain rate tensor",
            "", "s-1", "s-1", 0);
  set_attrs("second eigenvalue of the horizontal, vertically-integrated strain rate tensor",
            "", "s-1", "s-1", 1);
}

PetscErrorCode PSB_strain_rates::compute(IceModelVec* &output) {
  PetscErrorCode ierr;
  IceModelVec2 *result;
  IceModelVec *velbar;
  IceModelVec2Int *mask;
  PSB_velbar diag(model, grid, variables);

  result = new IceModelVec2;
  ierr = result->create(grid, "strain_rates", false, 1, 2); CHKERRQ(ierr);
  ierr = result->set_metadata(vars[0], 0); CHKERRQ(ierr);
  ierr = result->set_metadata(vars[1], 1); CHKERRQ(ierr);

  mask = dynamic_cast<IceModelVec2Int*>(variables.get("mask"));
  if (mask == NULL) SETERRQ(grid.com, 1, "mask is not available");

  ierr = diag.compute(velbar); CHKERRQ(ierr);
  IceModelVec2V *v_tmp = dynamic_cast<IceModelVec2V*>(velbar);
  if (v_tmp == NULL) SETERRQ(grid.com, 1, "velbar is expected to be an IceModelVec2V");

  IceModelVec2V velbar_with_ghosts;
  ierr = velbar_with_ghosts.create(grid, "velbar", true); CHKERRQ(ierr);

  // copy_from communicates ghosts
  ierr = velbar_with_ghosts.copy_from(*v_tmp); CHKERRQ(ierr);

  ierr = model->compute_2D_principal_strain_rates(velbar_with_ghosts, *mask, *result); CHKERRQ(ierr);

  delete velbar;
  output = result;
  return 0;
}

PSB_deviatoric_stresses::PSB_deviatoric_stresses(PISMStressBalance *m, IceGrid &g, PISMVars &my_vars)
  : PISMDiag<PISMStressBalance>(m, g, my_vars) {
  dof = 3;
  vars.resize(dof);

  // set metadata:
  vars[0].init_2d("sigma_xx", grid);
  vars[1].init_2d("sigma_yy", grid);
  vars[2].init_2d("sigma_xy", grid);

  set_attrs("deviatoric stress in X direction", "", "Pa", "Pa", 0);
  set_attrs("deviatoric stress in Y direction", "", "Pa", "Pa", 1);
  set_attrs("deviatoric shear stress", "", "Pa", "Pa", 2);

}

PetscErrorCode PSB_deviatoric_stresses::compute(IceModelVec* &output) {
  PetscErrorCode ierr;
  IceModelVec2 *result;
  IceModelVec *velbar;
  IceModelVec2Int *mask;
  PSB_velbar diag(model, grid, variables);

  result = new IceModelVec2;
  ierr = result->create(grid, "strain_rates", false, 1, 3); CHKERRQ(ierr);
  ierr = result->set_metadata(vars[0], 0); CHKERRQ(ierr);
  ierr = result->set_metadata(vars[1], 1); CHKERRQ(ierr);
  ierr = result->set_metadata(vars[2], 2); CHKERRQ(ierr);

  mask = dynamic_cast<IceModelVec2Int*>(variables.get("mask"));
  if (mask == NULL) SETERRQ(grid.com, 1, "mask is not available");

  ierr = diag.compute(velbar); CHKERRQ(ierr);
  IceModelVec2V *v_tmp = dynamic_cast<IceModelVec2V*>(velbar);
  if (v_tmp == NULL) SETERRQ(grid.com, 1, "velbar is expected to be an IceModelVec2V");

  IceModelVec2V velbar_with_ghosts;
  ierr = velbar_with_ghosts.create(grid, "velbar", true); CHKERRQ(ierr);

  // copy_from communicates ghosts
  ierr = velbar_with_ghosts.copy_from(*v_tmp); CHKERRQ(ierr);

  ierr = model->compute_2D_stresses(velbar_with_ghosts, *mask, *result); CHKERRQ(ierr);

  output = result;

  return 0;
}

