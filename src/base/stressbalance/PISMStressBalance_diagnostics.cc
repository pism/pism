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

#include "PISMStressBalance_diagnostics.hh"
#include "Mask.hh"
#include "ShallowStressBalance.hh"
#include "SSB_Modifier.hh"
#include "PISMVars.hh"
#include "PISMConfig.hh"

#include "error_handling.hh"

namespace pism {

void StressBalance::get_diagnostics(std::map<std::string, Diagnostic*> &dict,
                                        std::map<std::string, TSDiagnostic*> &ts_dict) {

  dict["bfrict"]   = new PSB_bfrict(this, grid, *m_variables);

  dict["velbar_mag"]     = new PSB_velbar_mag(this,     grid, *m_variables);
  dict["flux_mag"]     = new PSB_flux_mag(this,     grid, *m_variables);
  dict["velbase_mag"]    = new PSB_velbase_mag(this,    grid, *m_variables);
  dict["velsurf_mag"]    = new PSB_velsurf_mag(this,    grid, *m_variables);

  dict["uvel"]     = new PSB_uvel(this, grid, *m_variables);
  dict["vvel"]     = new PSB_vvel(this, grid, *m_variables);

  dict["strainheat"] = new PSB_strainheat(this, grid, *m_variables);

  dict["velbar"]   = new PSB_velbar(this,   grid, *m_variables);
  dict["velbase"]  = new PSB_velbase(this,  grid, *m_variables);
  dict["velsurf"]  = new PSB_velsurf(this,  grid, *m_variables);

  dict["wvel"]     = new PSB_wvel(this,     grid, *m_variables);
  dict["wvelbase"] = new PSB_wvelbase(this, grid, *m_variables);
  dict["wvelsurf"] = new PSB_wvelsurf(this, grid, *m_variables);
  dict["wvel_rel"] = new PSB_wvel_rel(this, grid, *m_variables);
  dict["strain_rates"] = new PSB_strain_rates(this, grid, *m_variables);
  dict["deviatoric_stresses"] = new PSB_deviatoric_stresses(this, grid, *m_variables);

  dict["pressure"] = new PSB_pressure(this, grid, *m_variables);
  dict["tauxz"] = new PSB_tauxz(this, grid, *m_variables);
  dict["tauyz"] = new PSB_tauyz(this, grid, *m_variables);

  m_stress_balance->get_diagnostics(dict, ts_dict);
  m_modifier->get_diagnostics(dict, ts_dict);
}

PSB_velbar::PSB_velbar(StressBalance *m, IceGrid &g, Vars &my_vars)
  : Diag<StressBalance>(m, g, my_vars) {

  dof = 2;

  // set metadata:
  vars.push_back(NCSpatialVariable(grid.get_unit_system(), "ubar", grid));
  vars.push_back(NCSpatialVariable(grid.get_unit_system(), "vbar", grid));

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
  IceModelVec2S *thickness = variables.get_2d_scalar("land_ice_thickness");
  IceModelVec2V *result;
  double *u_ij, *v_ij;
  double icefree_thickness = grid.config.get("mask_icefree_thickness_standard");

  result = new IceModelVec2V;
  ierr = result->create(grid, "bar", WITHOUT_GHOSTS); CHKERRQ(ierr);
  result->metadata() = vars[0];
  result->metadata(1) = vars[1];

  ierr = model->get_3d_velocity(u3, v3, w3); CHKERRQ(ierr);

  IceModelVec::AccessList list;
  list.add(*u3);
  list.add(*v3);
  list.add(*thickness);
  list.add(*result);

  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    double u_sum = 0, v_sum = 0,
      thk = (*thickness)(i,j);
    int ks = grid.kBelowHeight(thk);

    // an "ice-free" cell:
    if (thk < icefree_thickness) {
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


  output = result;
  return 0;
}

PSB_velbar_mag::PSB_velbar_mag(StressBalance *m, IceGrid &g, Vars &my_vars)
  : Diag<StressBalance>(m, g, my_vars) {

  // set metadata:
  vars.push_back(NCSpatialVariable(grid.get_unit_system(), "velbar_mag", grid));

  set_attrs("magnitude of vertically-integrated horizontal velocity of ice", "",
            "m s-1", "m year-1", 0);
  vars[0].set_double("_FillValue", grid.convert(grid.config.get("fill_value"),
                                         "m/year", "m/s"));
  vars[0].set_double("valid_min", 0.0);
}

PetscErrorCode PSB_velbar_mag::compute(IceModelVec* &output) {
  PetscErrorCode ierr;
  IceModelVec *tmp;
  IceModelVec2V *velbar_vec;
  IceModelVec2S *thickness = variables.get_2d_scalar("land_ice_thickness"),
    *result = NULL;

  result = new IceModelVec2S;
  ierr = result->create(grid, "velbar_mag", WITHOUT_GHOSTS);
  result->metadata() = vars[0];

  // compute vertically-averaged horizontal velocity:
  PSB_velbar velbar(model, grid, variables);
  ierr = velbar.compute(tmp); CHKERRQ(ierr);

  velbar_vec = dynamic_cast<IceModelVec2V*>(tmp);
  if (velbar_vec == NULL) {
    throw RuntimeError("dynamic cast failure");
  }

  // compute its magnitude:
  ierr = velbar_vec->magnitude(*result); CHKERRQ(ierr);

  // mask out ice-free areas:
  ierr = result->mask_by(*thickness, grid.config.get("fill_value", "m/year", "m/s")); CHKERRQ(ierr);

  delete tmp;
  output = result;
  return 0;
}

PSB_flux_mag::PSB_flux_mag(StressBalance *m, IceGrid &g, Vars &my_vars)
  : Diag<StressBalance>(m, g, my_vars) {

  // set metadata:
  vars.push_back(NCSpatialVariable(grid.get_unit_system(), "flux_mag", grid));

  set_attrs("magnitude of vertically-integrated horizontal flux of ice", "",
            "m2 s-1", "m2 year-1", 0);
  vars[0].set_double("_FillValue", grid.config.get("fill_value", "m2/year", "m2/s"));
  vars[0].set_double("valid_min", 0.0);
}

PetscErrorCode PSB_flux_mag::compute(IceModelVec* &output) {
  PetscErrorCode ierr;
  IceModelVec2S *thickness = variables.get_2d_scalar("land_ice_thickness"),
    *result = NULL;
  IceModelVec *tmp = NULL;

  // Compute the vertically-average horizontal ice velocity:
  PSB_velbar_mag velbar_mag(model, grid, variables);
  ierr = velbar_mag.compute(tmp); CHKERRQ(ierr);
  // NB: the call above allocates memory

  result = dynamic_cast<IceModelVec2S*>(tmp);
  if (result == NULL) {
    throw RuntimeError("dynamic_cast failure");
  }

  IceModelVec::AccessList list;
  list.add(*thickness);
  list.add(*result);

  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    (*result)(i,j) *= (*thickness)(i,j);
  }


  ierr = result->mask_by(*thickness, grid.config.get("fill_value", "m/year", "m/s")); CHKERRQ(ierr);

  result->metadata() = vars[0];

  output = result;
  return 0;
}

PSB_velbase_mag::PSB_velbase_mag(StressBalance *m, IceGrid &g, Vars &my_vars)
 : Diag<StressBalance>(m, g, my_vars) {

  // set metadata:
  vars.push_back(NCSpatialVariable(grid.get_unit_system(), "velbase_mag", grid));

  set_attrs("magnitude of horizontal velocity of ice at base of ice", "",
            "m s-1", "m year-1", 0);
  vars[0].set_double("_FillValue", grid.config.get("fill_value", "m/year", "m/s"));
  vars[0].set_double("valid_min", 0.0);
}

PetscErrorCode PSB_velbase_mag::compute(IceModelVec* &output) {
  PetscErrorCode ierr;
  IceModelVec3 *u3, *v3, *w3;
  IceModelVec2S tmp, *result, *thickness;

  ierr = tmp.create(grid, "tmp", WITHOUT_GHOSTS); CHKERRQ(ierr);

  result = new IceModelVec2S;
  ierr = result->create(grid, "velbase_mag", WITHOUT_GHOSTS); CHKERRQ(ierr);
  result->metadata() = vars[0];

  ierr = model->get_3d_velocity(u3, v3, w3); CHKERRQ(ierr);

  thickness = variables.get_2d_scalar("land_ice_thickness");

  ierr = u3->getHorSlice(*result, 0.0); CHKERRQ(ierr); // result = u_{z=0}
  ierr = v3->getHorSlice(tmp, 0.0); CHKERRQ(ierr);    // tmp = v_{z=0}

  ierr = result->set_to_magnitude(*result, tmp); CHKERRQ(ierr);

  // mask out ice-free areas
  ierr = result->mask_by(*thickness, grid.config.get("fill_value", "m/year", "m/s")); CHKERRQ(ierr);

  output = result;
  return 0;
}

PSB_velsurf_mag::PSB_velsurf_mag(StressBalance *m, IceGrid &g, Vars &my_vars)
  : Diag<StressBalance>(m, g, my_vars) {
  // set metadata:
  vars.push_back(NCSpatialVariable(grid.get_unit_system(), "velsurf_mag", grid));

  set_attrs("magnitude of horizontal velocity of ice at ice surface", "",
            "m s-1", "m year-1", 0);
  vars[0].set_double("_FillValue", grid.config.get("fill_value", "m/year", "m/s"));
  vars[0].set_double("valid_min",  0.0);
}

PetscErrorCode PSB_velsurf_mag::compute(IceModelVec* &output) {
  PetscErrorCode ierr;

  IceModelVec3 *u3, *v3, *w3;
  IceModelVec2S tmp, *result, *thickness;

  ierr = tmp.create(grid, "tmp", WITHOUT_GHOSTS); CHKERRQ(ierr);

  result = new IceModelVec2S;
  ierr = result->create(grid, "velsurf_mag", WITHOUT_GHOSTS); CHKERRQ(ierr);
  result->metadata() = vars[0];

  ierr = model->get_3d_velocity(u3, v3, w3); CHKERRQ(ierr);

  thickness = variables.get_2d_scalar("land_ice_thickness");

  ierr = u3->getSurfaceValues(*result, *thickness); CHKERRQ(ierr);
  ierr = v3->getSurfaceValues(tmp, *thickness); CHKERRQ(ierr);

  ierr = result->set_to_magnitude(*result, tmp); CHKERRQ(ierr);

  // mask out ice-free areas
  ierr = result->mask_by(*thickness, grid.config.get("fill_value", "m/year", "m/s")); CHKERRQ(ierr);

  output = result;
  return 0;
}


PSB_velsurf::PSB_velsurf(StressBalance *m, IceGrid &g, Vars &my_vars)
  : Diag<StressBalance>(m, g, my_vars) {

  // set metadata:
  dof = 2;

  vars.push_back(NCSpatialVariable(grid.get_unit_system(), "uvelsurf", grid));
  vars.push_back(NCSpatialVariable(grid.get_unit_system(), "vvelsurf", grid));

  set_attrs("x-component of the horizontal velocity of ice at ice surface", "",
            "m s-1", "m year-1", 0);
  set_attrs("y-component of the horizontal velocity of ice at ice surface", "",
            "m s-1", "m year-1", 1);

  vars[0].set_double("valid_min", grid.convert(-1e6, "m/year", "m/second"));
  vars[0].set_double("valid_max", grid.convert(1e6, "m/year", "m/second"));
  vars[0].set_double("_FillValue", grid.config.get("fill_value", "m/year", "m/s"));

  vars[1].set_double("valid_min", grid.convert(-1e6, "m/year", "m/second"));
  vars[1].set_double("valid_max", grid.convert(1e6, "m/year", "m/second"));
  vars[1].set_double("_FillValue", grid.config.get("fill_value", "m/year", "m/s"));
}

PetscErrorCode PSB_velsurf::compute(IceModelVec* &output) {
  PetscErrorCode ierr;
  IceModelVec2V *result;
  IceModelVec3 *u3, *v3, *w3;
  IceModelVec2S *thickness, tmp;
  double fill_value = grid.config.get("fill_value", "m/year", "m/s");

  result = new IceModelVec2V;
  ierr = result->create(grid, "surf", WITHOUT_GHOSTS); CHKERRQ(ierr);
  result->metadata(0) = vars[0];
  result->metadata(1) = vars[1];

  ierr = tmp.create(grid, "tmp", WITHOUT_GHOSTS); CHKERRQ(ierr);

  ierr = model->get_3d_velocity(u3, v3, w3); CHKERRQ(ierr);

  thickness = variables.get_2d_scalar("land_ice_thickness");

  ierr = u3->getSurfaceValues(tmp, *thickness); CHKERRQ(ierr);
  ierr = result->set_component(0, tmp); CHKERRQ(ierr);

  ierr = v3->getSurfaceValues(tmp, *thickness); CHKERRQ(ierr);
  ierr = result->set_component(1, tmp); CHKERRQ(ierr);

  IceModelVec2Int *mask = variables.get_2d_mask("mask");

  MaskQuery M(*mask);

  IceModelVec::AccessList list;
  list.add(*mask);
  list.add(*result);

  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (M.ice_free(i, j)) {
      (*result)(i, j).u = fill_value;
      (*result)(i, j).v = fill_value;
    }
  }


  output = result;
  return 0;
}

PSB_wvel::PSB_wvel(StressBalance *m, IceGrid &g, Vars &my_vars)
  : Diag<StressBalance>(m, g, my_vars) {

  // set metadata:
  vars.push_back(NCSpatialVariable(grid.get_unit_system(), "wvel", grid, g.zlevels));

  set_attrs("vertical velocity of ice, relative to geoid", "",
            "m s-1", "m year-1", 0);
  vars[0].set_double("valid_min", grid.convert(-1e6, "m/year", "m/second"));
  vars[0].set_double("valid_max", grid.convert(1e6, "m/year", "m/second"));
}

PetscErrorCode PSB_wvel::compute(IceModelVec* &output) {
  PetscErrorCode ierr;
  IceModelVec3 *result, *u3, *v3, *w3;
  IceModelVec2S *bed, *uplift, *thickness;
  IceModelVec2Int *mask;
  double *u, *v, *w, *res;

  result = new IceModelVec3;
  ierr = result->create(grid, "wvel", WITHOUT_GHOSTS); CHKERRQ(ierr);
  result->metadata() = vars[0];

  mask      = variables.get_2d_mask("mask");
  bed       = variables.get_2d_scalar("bedrock_altitude");
  uplift    = variables.get_2d_scalar("tendency_of_bedrock_altitude");
  thickness = variables.get_2d_scalar("land_ice_thickness");

  ierr = model->get_3d_velocity(u3, v3, w3); CHKERRQ(ierr);

  IceModelVec::AccessList list;
  list.add(*thickness);
  list.add(*mask);
  list.add(*bed);
  list.add(*u3);
  list.add(*v3);
  list.add(*w3);
  list.add(*uplift);
  list.add(*result);

  MaskQuery M(*mask);

  const double ice_density = grid.config.get("ice_density"),
    sea_water_density = grid.config.get("sea_water_density"),
    R = ice_density / sea_water_density;

  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    ierr = u3->getInternalColumn(i, j, &u); CHKERRQ(ierr);
    ierr = v3->getInternalColumn(i, j, &v); CHKERRQ(ierr);
    ierr = w3->getInternalColumn(i, j, &w); CHKERRQ(ierr);
    ierr = result->getInternalColumn(i, j, &res); CHKERRQ(ierr);

    int ks = grid.kBelowHeight((*thickness)(i,j));

    // in the ice:
    if (M.grounded(i,j)) {
      for (int k = 0; k <= ks ; k++) {
        res[k] = w[k] + (*uplift)(i,j) + u[k] * bed->diff_x_p(i,j) + v[k] * bed->diff_y_p(i,j);
      }

    } else {                  // floating
      const double
        z_sl = R * (*thickness)(i,j),
        w_sl = w3->getValZ(i, j, z_sl);

      for (int k = 0; k <= ks ; k++) {
        res[k] = w[k] - w_sl;
      }

    }

    // above the ice:
    for (unsigned int k = ks+1; k < grid.Mz ; k++) {
      res[k] = 0.0;
    }

  }

  output = result;
  return 0;
}

PSB_wvelsurf::PSB_wvelsurf(StressBalance *m, IceGrid &g, Vars &my_vars)
  : Diag<StressBalance>(m, g, my_vars) {

  // set metadata:
  vars.push_back(NCSpatialVariable(grid.get_unit_system(), "wvelsurf", grid));

  set_attrs("vertical velocity of ice at ice surface, relative to the geoid", "",
            "m s-1", "m year-1", 0);
  vars[0].set_double("valid_min", grid.convert(-1e6, "m/year", "m/second"));
  vars[0].set_double("valid_max", grid.convert(1e6, "m/year", "m/second"));
  vars[0].set_double("_FillValue", grid.config.get("fill_value", "m/year", "m/s"));
}

PetscErrorCode PSB_wvelsurf::compute(IceModelVec* &output) {
  PetscErrorCode ierr;
  IceModelVec *tmp;
  IceModelVec3 *w3;
  IceModelVec2S *result, *thickness;
  double fill_value = grid.config.get("fill_value", "m/year", "m/s");

  result = new IceModelVec2S;
  ierr = result->create(grid, "wvelsurf", WITHOUT_GHOSTS); CHKERRQ(ierr);
  result->metadata() = vars[0];

  PSB_wvel wvel(model, grid, variables);

  ierr = wvel.compute(tmp); CHKERRQ(ierr);

  w3 = dynamic_cast<IceModelVec3*>(tmp);
  if (tmp == NULL) {
    throw RuntimeError("dynamic_cast failure");
  }

  thickness = variables.get_2d_scalar("land_ice_thickness");

  ierr = w3->getSurfaceValues(*result, *thickness); CHKERRQ(ierr);

  IceModelVec2Int *mask = variables.get_2d_mask("mask");

  MaskQuery M(*mask);

  IceModelVec::AccessList list;
  list.add(*mask);
  list.add(*result);

  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (M.ice_free(i, j)) {
      (*result)(i, j) = fill_value;
    }
  }


  delete tmp;
  output = result;
  return 0;
}

PSB_wvelbase::PSB_wvelbase(StressBalance *m, IceGrid &g, Vars &my_vars)
  : Diag<StressBalance>(m, g, my_vars) {

  // set metadata:
  vars.push_back(NCSpatialVariable(grid.get_unit_system(), "wvelbase", grid));

  set_attrs("vertical velocity of ice at the base of ice, relative to the geoid", "",
            "m s-1", "m year-1", 0);
  vars[0].set_double("valid_min", grid.convert(-1e6, "m/year", "m/second"));
  vars[0].set_double("valid_max", grid.convert(1e6, "m/year", "m/second"));
  vars[0].set_double("_FillValue", grid.config.get("fill_value", "m/year", "m/s"));
}

PetscErrorCode PSB_wvelbase::compute(IceModelVec* &output) {
  PetscErrorCode ierr;
  IceModelVec *tmp;
  IceModelVec3 *w3;
  IceModelVec2S *result;
  double fill_value = grid.config.get("fill_value", "m/year", "m/s");

  result = new IceModelVec2S;
  ierr = result->create(grid, "wvelbase", WITHOUT_GHOSTS); CHKERRQ(ierr);
  result->metadata() = vars[0];

  PSB_wvel wvel(model, grid, variables);

  ierr = wvel.compute(tmp); CHKERRQ(ierr);

  w3 = dynamic_cast<IceModelVec3*>(tmp);
  if (tmp == NULL) {
    throw RuntimeError("dynamic_cast failure");
  }

  ierr = w3->getHorSlice(*result, 0.0); CHKERRQ(ierr);

  IceModelVec2Int *mask = variables.get_2d_mask("mask");

  MaskQuery M(*mask);

  IceModelVec::AccessList list;
  list.add(*mask);
  list.add(*result);

  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (M.ice_free(i, j)) {
      (*result)(i, j) = fill_value;
    }
  }


  delete tmp;
  output = result;
  return 0;
}

PSB_velbase::PSB_velbase(StressBalance *m, IceGrid &g, Vars &my_vars)
  : Diag<StressBalance>(m, g, my_vars) {

  // set metadata:
  dof = 2;

  vars.push_back(NCSpatialVariable(grid.get_unit_system(), "uvelbase", grid));
  vars.push_back(NCSpatialVariable(grid.get_unit_system(), "vvelbase", grid));

  set_attrs("x-component of the horizontal velocity of ice at the base of ice", "",
            "m s-1", "m year-1", 0);
  set_attrs("y-component of the horizontal velocity of ice at the base of ice", "",
            "m s-1", "m year-1", 1);

  vars[0].set_double("valid_min", grid.convert(-1e6, "m/year", "m/second"));
  vars[0].set_double("valid_max", grid.convert(1e6, "m/year", "m/second"));
  vars[0].set_double("_FillValue", grid.config.get("fill_value", "m/year", "m/s"));

  vars[1].set_double("valid_min", grid.convert(-1e6, "m/year", "m/second"));
  vars[1].set_double("valid_max", grid.convert(1e6, "m/year", "m/second"));
  vars[1].set_double("_FillValue", grid.config.get("fill_value", "m/year", "m/s"));
}

PetscErrorCode PSB_velbase::compute(IceModelVec* &output) {
  PetscErrorCode ierr;
  IceModelVec2V *result;
  IceModelVec3 *u3, *v3, *w3;
  IceModelVec2S tmp;            // will be de-allocated automatically
  double fill_value = grid.config.get("fill_value", "m/year", "m/s");

  result = new IceModelVec2V;
  ierr = result->create(grid, "base", WITHOUT_GHOSTS); CHKERRQ(ierr);
  result->metadata() = vars[0];
  result->metadata(1) = vars[1];

  ierr = tmp.create(grid, "tmp", WITHOUT_GHOSTS); CHKERRQ(ierr);

  ierr = model->get_3d_velocity(u3, v3, w3); CHKERRQ(ierr);

  ierr = u3->getHorSlice(tmp, 0.0); CHKERRQ(ierr);
  ierr = result->set_component(0, tmp); CHKERRQ(ierr);

  ierr = v3->getHorSlice(tmp, 0.0); CHKERRQ(ierr);
  ierr = result->set_component(1, tmp); CHKERRQ(ierr);

  IceModelVec2Int *mask = variables.get_2d_mask("mask");

  MaskQuery M(*mask);

  IceModelVec::AccessList list;
  list.add(*mask);
  list.add(*result);

  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (M.ice_free(i, j)) {
      (*result)(i, j).u = fill_value;
      (*result)(i, j).v = fill_value;
    }
  }


  output = result;
  return 0;
}


PSB_bfrict::PSB_bfrict(StressBalance *m, IceGrid &g, Vars &my_vars)
  : Diag<StressBalance>(m, g, my_vars) {

  // set metadata:
  vars.push_back(NCSpatialVariable(grid.get_unit_system(), "bfrict", grid));

  set_attrs("basal frictional heating", "",
            "W m-2", "W m-2", 0);
}

PetscErrorCode PSB_bfrict::compute(IceModelVec* &output) {
  PetscErrorCode ierr;

  IceModelVec2S *result = new IceModelVec2S;
  ierr = result->create(grid, "bfrict", WITHOUT_GHOSTS); CHKERRQ(ierr);
  result->metadata() = vars[0];

  IceModelVec2S *bfrict;
  ierr = model->get_basal_frictional_heating(bfrict); CHKERRQ(ierr);

  ierr = bfrict->copy_to(*result); CHKERRQ(ierr);

  output = result;
  return 0;
}


PSB_uvel::PSB_uvel(StressBalance *m, IceGrid &g, Vars &my_vars)
  : Diag<StressBalance>(m, g, my_vars) {

  // set metadata:
  vars.push_back(NCSpatialVariable(grid.get_unit_system(), "uvel", grid, g.zlevels));

  set_attrs("horizontal velocity of ice in the X direction", "land_ice_x_velocity",
            "m s-1", "m year-1", 0);
}

PetscErrorCode PSB_uvel::compute(IceModelVec* &output) {
  PetscErrorCode ierr;

  IceModelVec3 *result = new IceModelVec3;
  ierr = result->create(grid, "uvel", WITHOUT_GHOSTS); CHKERRQ(ierr);
  result->metadata() = vars[0];

  IceModelVec2S *thickness = variables.get_2d_scalar("land_ice_thickness");

  IceModelVec3 *u3, *v3, *w3;
  ierr = model->get_3d_velocity(u3, v3, w3); CHKERRQ(ierr);

  IceModelVec::AccessList list;
  list.add(*u3);
  list.add(*result);
  list.add(*thickness);

  double *u_ij, *u_out_ij;
  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    int ks = grid.kBelowHeight((*thickness)(i,j));

    ierr = u3->getInternalColumn(i,j,&u_ij); CHKERRQ(ierr);
    ierr = result->getInternalColumn(i,j,&u_out_ij); CHKERRQ(ierr);

    // in the ice:
    for (int k = 0; k <= ks ; k++) {
      u_out_ij[k] = u_ij[k];
    }
    // above the ice:
    for (unsigned int k = ks+1; k < grid.Mz ; k++) {
      u_out_ij[k] = 0.0;
    }
  }


  output = result;
  return 0;
}

PSB_vvel::PSB_vvel(StressBalance *m, IceGrid &g, Vars &my_vars)
  : Diag<StressBalance>(m, g, my_vars) {

  // set metadata:
  vars.push_back(NCSpatialVariable(grid.get_unit_system(), "vvel", grid, g.zlevels));

  set_attrs("horizontal velocity of ice in the Y direction", "land_ice_y_velocity",
            "m s-1", "m year-1", 0);
}

PetscErrorCode PSB_vvel::compute(IceModelVec* &output) {
  PetscErrorCode ierr;

  IceModelVec3 *result = new IceModelVec3;
  ierr = result->create(grid, "vvel", WITHOUT_GHOSTS); CHKERRQ(ierr);
  result->metadata() = vars[0];

  IceModelVec2S *thickness = variables.get_2d_scalar("land_ice_thickness");

  IceModelVec3 *u3, *v3, *w3;
  ierr = model->get_3d_velocity(u3, v3, w3); CHKERRQ(ierr);

  IceModelVec::AccessList list;
  list.add(*v3);
  list.add(*result);
  list.add(*thickness);

  double *v_ij, *v_out_ij;
  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    int ks = grid.kBelowHeight((*thickness)(i,j));

    ierr = v3->getInternalColumn(i,j,&v_ij); CHKERRQ(ierr);
    ierr = result->getInternalColumn(i,j,&v_out_ij); CHKERRQ(ierr);

    // in the ice:
    for (int k = 0; k <= ks ; k++) {
      v_out_ij[k] = v_ij[k];
    }
    // above the ice:
    for (unsigned int k = ks+1; k < grid.Mz ; k++) {
      v_out_ij[k] = 0.0;
    }
  }


  output = result;
  return 0;
}

PSB_wvel_rel::PSB_wvel_rel(StressBalance *m, IceGrid &g, Vars &my_vars)
  : Diag<StressBalance>(m, g, my_vars) {

  // set metadata:
  vars.push_back(NCSpatialVariable(grid.get_unit_system(), "wvel_rel", grid, g.zlevels));

  set_attrs("vertical velocity of ice, relative to base of ice directly below", "",
            "m s-1", "m year-1", 0);
}

PetscErrorCode PSB_wvel_rel::compute(IceModelVec* &output) {
  PetscErrorCode ierr;

  IceModelVec3 *result = new IceModelVec3;
  ierr = result->create(grid, "wvel_rel", WITHOUT_GHOSTS); CHKERRQ(ierr);
  result->metadata() = vars[0];

  IceModelVec2S *thickness = variables.get_2d_scalar("land_ice_thickness");

  IceModelVec3 *u3, *v3, *w3;
  ierr = model->get_3d_velocity(u3, v3, w3); CHKERRQ(ierr);

  IceModelVec::AccessList list;
  list.add(*w3);
  list.add(*result);
  list.add(*thickness);

  double *w_ij, *w_out_ij;
  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    int ks = grid.kBelowHeight((*thickness)(i,j));

    ierr = w3->getInternalColumn(i,j,&w_ij); CHKERRQ(ierr);
    ierr = result->getInternalColumn(i,j,&w_out_ij); CHKERRQ(ierr);

    // in the ice:
    for (int k = 0; k <= ks ; k++) {
      w_out_ij[k] = w_ij[k];
    }
    // above the ice:
    for (unsigned int k = ks+1; k < grid.Mz ; k++) {
      w_out_ij[k] = 0.0;
    }
  }


  output = result;
  return 0;
}


PSB_strainheat::PSB_strainheat(StressBalance *m, IceGrid &g, Vars &my_vars)
  : Diag<StressBalance>(m, g, my_vars) {

  // set metadata:
  vars.push_back(NCSpatialVariable(grid.get_unit_system(), "strainheat", grid, grid.zlevels));

  set_attrs("rate of strain heating in ice (dissipation heating)", "",
            "W m-3", "mW m-3", 0);
}

PetscErrorCode PSB_strainheat::compute(IceModelVec* &output) {
  PetscErrorCode ierr;

  IceModelVec3 *result = new IceModelVec3;
  ierr = result->create(grid, "strainheat", WITHOUT_GHOSTS); CHKERRQ(ierr);
  result->metadata() = vars[0];
  result->write_in_glaciological_units = true;

  IceModelVec3 *tmp;
  ierr = model->get_volumetric_strain_heating(tmp); CHKERRQ(ierr);

  ierr = tmp->copy_to(*result); CHKERRQ(ierr);

  output = result;
  return 0;
}

PSB_strain_rates::PSB_strain_rates(StressBalance *m, IceGrid &g, Vars &my_vars)
  : Diag<StressBalance>(m, g, my_vars) {
  dof = 2;

  // set metadata:
  vars.push_back(NCSpatialVariable(grid.get_unit_system(), "eigen1", grid));
  vars.push_back(NCSpatialVariable(grid.get_unit_system(), "eigen2", grid));

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
  ierr = result->create(grid, "strain_rates", WITHOUT_GHOSTS, 1, 2); CHKERRQ(ierr);
  result->metadata() = vars[0];
  result->metadata(1) = vars[1];

  mask = variables.get_2d_mask("mask");

  ierr = diag.compute(velbar); CHKERRQ(ierr);
  IceModelVec2V *v_tmp = dynamic_cast<IceModelVec2V*>(velbar);
  if (v_tmp == NULL) {
    throw RuntimeError("velbar is expected to be an IceModelVec2V");
  }

  IceModelVec2V velbar_with_ghosts;
  ierr = velbar_with_ghosts.create(grid, "velbar", WITH_GHOSTS); CHKERRQ(ierr);

  // copy_from communicates ghosts
  ierr = velbar_with_ghosts.copy_from(*v_tmp); CHKERRQ(ierr);

  ierr = model->compute_2D_principal_strain_rates(velbar_with_ghosts, *mask, *result); CHKERRQ(ierr);

  delete velbar;
  output = result;
  return 0;
}

PSB_deviatoric_stresses::PSB_deviatoric_stresses(StressBalance *m, IceGrid &g, Vars &my_vars)
  : Diag<StressBalance>(m, g, my_vars) {
  dof = 3;

  // set metadata:
  vars.push_back(NCSpatialVariable(grid.get_unit_system(), "sigma_xx", grid));
  vars.push_back(NCSpatialVariable(grid.get_unit_system(), "sigma_yy", grid));
  vars.push_back(NCSpatialVariable(grid.get_unit_system(), "sigma_xy", grid));

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
  ierr = result->create(grid, "deviatoric_stresses", WITHOUT_GHOSTS, 1, 3); CHKERRQ(ierr);
  result->metadata() = vars[0];
  result->metadata(1) = vars[1];
  result->metadata(2) = vars[2];

  mask = variables.get_2d_mask("mask");

  ierr = diag.compute(velbar); CHKERRQ(ierr);
  IceModelVec2V *v_tmp = dynamic_cast<IceModelVec2V*>(velbar);
  if (v_tmp == NULL) {
    throw RuntimeError("velbar is expected to be an IceModelVec2V");
  }

  IceModelVec2V velbar_with_ghosts;
  ierr = velbar_with_ghosts.create(grid, "velbar", WITH_GHOSTS); CHKERRQ(ierr);

  // copy_from communicates ghosts
  ierr = velbar_with_ghosts.copy_from(*v_tmp); CHKERRQ(ierr);

  ierr = model->compute_2D_stresses(velbar_with_ghosts, *mask, *result); CHKERRQ(ierr);

  output = result;

  return 0;
}

PSB_pressure::PSB_pressure(StressBalance *m, IceGrid &g, Vars &my_vars)
  : Diag<StressBalance>(m, g, my_vars) {

  // set metadata:
  vars.push_back(NCSpatialVariable(grid.get_unit_system(), "pressure", grid, grid.zlevels));

  set_attrs("pressure in ice (hydrostatic)", "",
            "Pa", "Pa", 0);
}

PetscErrorCode PSB_pressure::compute(IceModelVec* &output) {
  PetscErrorCode ierr;

  IceModelVec3 *result = new IceModelVec3;
  ierr = result->create(grid, "pressure", WITHOUT_GHOSTS); CHKERRQ(ierr);
  result->metadata() = vars[0];

  IceModelVec2S *thickness = variables.get_2d_scalar("land_ice_thickness");

  IceModelVec::AccessList list;
  list.add(*result);
  list.add(*thickness);

  double *P_out_ij;
  const double rg = grid.config.get("ice_density") * grid.config.get("standard_gravity");

  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    unsigned int ks = grid.kBelowHeight((*thickness)(i,j));
    ierr = result->getInternalColumn(i,j,&P_out_ij); CHKERRQ(ierr);
    const double H = (*thickness)(i,j);
    // within the ice:
    for (unsigned int k = 0; k <= ks; ++k) {
      P_out_ij[k] = rg * (H - grid.zlevels[k]);  // FIXME: add atmospheric pressure?
    }
    // above the ice:
    for (unsigned int k = ks + 1; k < grid.Mz; ++k) {
      P_out_ij[k] = 0.0;  // FIXME: use atmospheric pressure?
    }
  }


  output = result;

  return 0;
}


PSB_tauxz::PSB_tauxz(StressBalance *m, IceGrid &g, Vars &my_vars)
  : Diag<StressBalance>(m, g, my_vars) {

  // set metadata:
  vars.push_back(NCSpatialVariable(grid.get_unit_system(), "tauxz", grid, grid.zlevels));

  set_attrs("shear stress xz component (in shallow ice approximation SIA)", "",
            "Pa", "Pa", 0);
}


/*!
 * The SIA-applicable shear stress component tauxz computed here is not used
 * by the model.  This implementation intentionally does not use the
 * eta-transformation or special cases at ice margins.
 * CODE DUPLICATION WITH PSB_tauyz
 */
PetscErrorCode PSB_tauxz::compute(IceModelVec* &output) {
  PetscErrorCode ierr;

  IceModelVec3 *result = new IceModelVec3;
  ierr = result->create(grid, "tauxz", WITHOUT_GHOSTS); CHKERRQ(ierr);
  result->metadata() = vars[0];

  IceModelVec2S *thickness, *surface;

  thickness = variables.get_2d_scalar("land_ice_thickness");
  surface = variables.get_2d_scalar("surface_altitude");

  IceModelVec::AccessList list;
  list.add(*result);
  list.add(*surface);
  list.add(*thickness);

  double *tauxz_out_ij;
  const double rg = grid.config.get("ice_density") * grid.config.get("standard_gravity");

  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();


    unsigned int ks = grid.kBelowHeight((*thickness)(i,j));
    ierr = result->getInternalColumn(i,j,&tauxz_out_ij); CHKERRQ(ierr);
    const double H    = (*thickness)(i,j),
      dhdx = surface->diff_x_p(i,j);
    // within the ice:
    for (unsigned int k = 0; k <= ks; ++k) {
      tauxz_out_ij[k] = - rg * (H - grid.zlevels[k]) * dhdx;
    }
    // above the ice:
    for (unsigned int k = ks + 1; k < grid.Mz; ++k) {
      tauxz_out_ij[k] = 0.0;
    }

  }


  output = result;
  return 0;
}


PSB_tauyz::PSB_tauyz(StressBalance *m, IceGrid &g, Vars &my_vars)
  : Diag<StressBalance>(m, g, my_vars) {

  // set metadata:
  vars.push_back(NCSpatialVariable(grid.get_unit_system(), "tauyz", grid, grid.zlevels));

  set_attrs("shear stress yz component (in shallow ice approximation SIA)", "",
            "Pa", "Pa", 0);
}


/*!
 * The SIA-applicable shear stress component tauyz computed here is not used
 * by the model.  This implementation intentionally does not use the
 * eta-transformation or special cases at ice margins.
 * CODE DUPLICATION WITH PSB_tauxz
 */
PetscErrorCode PSB_tauyz::compute(IceModelVec* &output) {
  PetscErrorCode ierr;

  IceModelVec3 *result = new IceModelVec3;
  ierr = result->create(grid, "tauyz", WITHOUT_GHOSTS); CHKERRQ(ierr);
  result->metadata() = vars[0];

  IceModelVec2S *thickness, *surface;

  thickness = variables.get_2d_scalar("land_ice_thickness");
  surface = variables.get_2d_scalar("surface_altitude");

  IceModelVec::AccessList list;
  list.add(*result);
  list.add(*surface);
  list.add(*thickness);

  double *tauyz_out_ij;
  const double rg = grid.config.get("ice_density") * grid.config.get("standard_gravity");

  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();


    unsigned int ks = grid.kBelowHeight((*thickness)(i,j));
    ierr = result->getInternalColumn(i,j,&tauyz_out_ij); CHKERRQ(ierr);
    const double H    = (*thickness)(i,j),
      dhdy = surface->diff_y_p(i,j);
    // within the ice:
    for (unsigned int k = 0; k <= ks; ++k) {
      tauyz_out_ij[k] = - rg * (H - grid.zlevels[k]) * dhdy;
    }
    // above the ice:
    for (unsigned int k = ks + 1; k < grid.Mz; ++k) {
      tauyz_out_ij[k] = 0.0;
    }

  }


  output = result;
  return 0;
}


} // end of namespace pism
