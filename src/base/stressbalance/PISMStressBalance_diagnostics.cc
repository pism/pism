// Copyright (C) 2010, 2011, 2012, 2013, 2014, 2015 Constantine Khroulev
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

  dict["bfrict"]   = new PSB_bfrict(this);

  dict["velbar_mag"]     = new PSB_velbar_mag(this);
  dict["flux"]           = new PSB_flux(this);
  dict["flux_mag"]       = new PSB_flux_mag(this);
  dict["velbase_mag"]    = new PSB_velbase_mag(this);
  dict["velsurf_mag"]    = new PSB_velsurf_mag(this);

  dict["uvel"]     = new PSB_uvel(this);
  dict["vvel"]     = new PSB_vvel(this);

  dict["strainheat"] = new PSB_strainheat(this);

  dict["velbar"]   = new PSB_velbar(this);
  dict["velbase"]  = new PSB_velbase(this);
  dict["velsurf"]  = new PSB_velsurf(this);

  dict["wvel"]     = new PSB_wvel(this);
  dict["wvelbase"] = new PSB_wvelbase(this);
  dict["wvelsurf"] = new PSB_wvelsurf(this);
  dict["wvel_rel"] = new PSB_wvel_rel(this);
  dict["strain_rates"] = new PSB_strain_rates(this);
  dict["deviatoric_stresses"] = new PSB_deviatoric_stresses(this);

  dict["pressure"] = new PSB_pressure(this);
  dict["tauxz"] = new PSB_tauxz(this);
  dict["tauyz"] = new PSB_tauyz(this);

  m_stress_balance->get_diagnostics(dict, ts_dict);
  m_modifier->get_diagnostics(dict, ts_dict);
}

PSB_velbar::PSB_velbar(StressBalance *m)
  : Diag<StressBalance>(m) {

  m_dof = 2;

  // set metadata:
  m_vars.push_back(NCSpatialVariable(m_grid.config.get_unit_system(), "ubar", m_grid));
  m_vars.push_back(NCSpatialVariable(m_grid.config.get_unit_system(), "vbar", m_grid));

  set_attrs("vertical mean of horizontal ice velocity in the X direction",
            "land_ice_vertical_mean_x_velocity",
            "m s-1", "m year-1", 0);
  set_attrs("vertical mean of horizontal ice velocity in the Y direction",
            "land_ice_vertical_mean_y_velocity",
            "m s-1", "m year-1", 1);
}

void PSB_velbar::compute(IceModelVec* &output) {
  IceModelVec2V *result;

  // get the thickness
  const IceModelVec2S *thickness = m_grid.variables().get_2d_scalar("land_ice_thickness");

  // Compute the vertically-integrated horizontal ice flux:
  PSB_flux flux(model);
  IceModelVec *tmp;
  flux.compute(tmp);
  result = dynamic_cast<IceModelVec2V*>(tmp);
  if (result == NULL) {
    throw RuntimeError("dynamic cast failure");
  }
  // NB: the call above allocates memory

  // Override metadata set by the flux computation
  result->metadata(0) = m_vars[0];
  result->metadata(1) = m_vars[1];

  IceModelVec::AccessList list;
  list.add(*thickness);
  list.add(*result);

  for (Points p(m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();
    double thk = (*thickness)(i,j);

    // Ice flux is masked already, but we need to check for division
    // by zero anyway.
    if (thk > 0.0) {
      (*result)(i,j) /= thk;
    } else {
      (*result)(i,j).u = 0.0;
      (*result)(i,j).v = 0.0;
    }
  }

  output = result;
}

PSB_velbar_mag::PSB_velbar_mag(StressBalance *m)
  : Diag<StressBalance>(m) {

  // set metadata:
  m_vars.push_back(NCSpatialVariable(m_grid.config.get_unit_system(), "velbar_mag", m_grid));

  set_attrs("magnitude of vertically-integrated horizontal velocity of ice", "",
            "m s-1", "m year-1", 0);
  m_vars[0].set_double("_FillValue", m_grid.convert(m_grid.config.get("fill_value"),
                                         "m/year", "m/s"));
  m_vars[0].set_double("valid_min", 0.0);
}

void PSB_velbar_mag::compute(IceModelVec* &output) {
  IceModelVec *tmp;
  IceModelVec2V *velbar_vec;

  IceModelVec2S *result = new IceModelVec2S;
  result->create(m_grid, "velbar_mag", WITHOUT_GHOSTS);
  result->metadata() = m_vars[0];

  // compute vertically-averaged horizontal velocity:
  PSB_velbar velbar(model);
  velbar.compute(tmp);

  velbar_vec = dynamic_cast<IceModelVec2V*>(tmp);
  if (velbar_vec == NULL) {
    throw RuntimeError("dynamic cast failure");
  }

  // compute its magnitude:
  velbar_vec->magnitude(*result);

  const IceModelVec2S *thickness = m_grid.variables().get_2d_scalar("land_ice_thickness");

  // mask out ice-free areas:
  result->mask_by(*thickness, m_grid.config.get("fill_value", "m/year", "m/s"));

  delete tmp;
  output = result;
}


PSB_flux::PSB_flux(StressBalance *m)
  : Diag<StressBalance>(m) {

  m_dof = 2;

  // set metadata:
  m_vars.push_back(NCSpatialVariable(m_grid.config.get_unit_system(), "uflux", m_grid));
  m_vars.push_back(NCSpatialVariable(m_grid.config.get_unit_system(), "vflux", m_grid));

  set_attrs("Vertically integrated horizontal flux of ice in the X direction",
            "",                 // no CF standard name
            "m2 s-1", "m2 year-1", 0);
  set_attrs("Vertically integrated horizontal flux of ice in the Y direction",
            "",                 // no CF standard name
            "m2 s-1", "m2 year-1", 1);
}

void PSB_flux::compute(IceModelVec* &output) {
  IceModelVec3 *u3, *v3, *w3;
  IceModelVec2V *result;
  double *u_ij, *v_ij;
  double icefree_thickness = m_grid.config.get("mask_icefree_thickness_standard");

  result = new IceModelVec2V;
  result->create(m_grid, "flux", WITHOUT_GHOSTS);
  result->metadata(0) = m_vars[0];
  result->metadata(1) = m_vars[1];

  // get the thickness
  const IceModelVec2S *thickness = m_grid.variables().get_2d_scalar("land_ice_thickness");

  model->get_3d_velocity(u3, v3, w3);

  IceModelVec::AccessList list;
  list.add(*u3);
  list.add(*v3);
  list.add(*thickness);
  list.add(*result);

  for (Points p(m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    double u_sum = 0, v_sum = 0,
      thk = (*thickness)(i,j);
    int ks = m_grid.kBelowHeight(thk);

    // an "ice-free" cell:
    if (thk < icefree_thickness) {
      (*result)(i,j).u = 0;
      (*result)(i,j).v = 0;
      continue;
    }

    // an ice-filled cell:
    u3->getInternalColumn(i, j, &u_ij);
    v3->getInternalColumn(i, j, &v_ij);

    if (thk <= m_grid.z(1)) {
      (*result)(i,j).u = u_ij[0];
      (*result)(i,j).v = v_ij[0];
      continue;
    }

    for (int k = 1; k <= ks; ++k) {
      u_sum += (m_grid.z(k) - m_grid.z(k-1)) * (u_ij[k] + u_ij[k-1]);
      v_sum += (m_grid.z(k) - m_grid.z(k-1)) * (v_ij[k] + v_ij[k-1]);
    }

    // Finish the trapezoidal rule integration (multiply by 1/2).
    (*result)(i,j).u = 0.5 * u_sum;
    (*result)(i,j).v = 0.5 * v_sum;

    // The top surface of the ice is not aligned with the grid, so
    // we have at most dz meters of ice above grid.z(ks).
    // Assume that its velocity is (u_ij[ks], v_ij[ks]) and add its
    // contribution.
    (*result)(i,j).u += u_ij[ks] * (thk - m_grid.z(ks));
    (*result)(i,j).v += v_ij[ks] * (thk - m_grid.z(ks));
  }

  output = result;
}


PSB_flux_mag::PSB_flux_mag(StressBalance *m)
  : Diag<StressBalance>(m) {

  // set metadata:
  m_vars.push_back(NCSpatialVariable(m_grid.config.get_unit_system(), "flux_mag", m_grid));

  set_attrs("magnitude of vertically-integrated horizontal flux of ice", "",
            "m2 s-1", "m2 year-1", 0);
  m_vars[0].set_double("_FillValue", m_grid.config.get("fill_value", "m2/year", "m2/s"));
  m_vars[0].set_double("valid_min", 0.0);
}

void PSB_flux_mag::compute(IceModelVec* &output) {
  const IceModelVec2S *thickness = m_grid.variables().get_2d_scalar("land_ice_thickness");
  IceModelVec *tmp = NULL;

  // Compute the vertically-averaged horizontal ice velocity:
  PSB_velbar_mag velbar_mag(model);
  velbar_mag.compute(tmp);
  // NB: the call above allocates memory

  IceModelVec2S *result = dynamic_cast<IceModelVec2S*>(tmp);
  if (result == NULL) {
    throw RuntimeError("dynamic_cast failure");
  }

  IceModelVec::AccessList list;
  list.add(*thickness);
  list.add(*result);

  for (Points p(m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    (*result)(i,j) *= (*thickness)(i,j);
  }


  result->mask_by(*thickness, m_grid.config.get("fill_value", "m/year", "m/s"));

  result->metadata() = m_vars[0];

  output = result;
}

PSB_velbase_mag::PSB_velbase_mag(StressBalance *m)
  : Diag<StressBalance>(m) {

  // set metadata:
  m_vars.push_back(NCSpatialVariable(m_grid.config.get_unit_system(), "velbase_mag", m_grid));

  set_attrs("magnitude of horizontal velocity of ice at base of ice", "",
            "m s-1", "m year-1", 0);
  m_vars[0].set_double("_FillValue", m_grid.config.get("fill_value", "m/year", "m/s"));
  m_vars[0].set_double("valid_min", 0.0);
}

void PSB_velbase_mag::compute(IceModelVec* &output) {
  IceModelVec3 *u3, *v3, *w3;
  IceModelVec2S tmp, *result;

  tmp.create(m_grid, "tmp", WITHOUT_GHOSTS);

  result = new IceModelVec2S;
  result->create(m_grid, "velbase_mag", WITHOUT_GHOSTS);
  result->metadata() = m_vars[0];

  model->get_3d_velocity(u3, v3, w3);

  const IceModelVec2S *thickness = m_grid.variables().get_2d_scalar("land_ice_thickness");

  u3->getHorSlice(*result, 0.0); // result = u_{z=0}
  v3->getHorSlice(tmp, 0.0);    // tmp = v_{z=0}

  result->set_to_magnitude(*result, tmp);

  // mask out ice-free areas
  result->mask_by(*thickness, m_grid.config.get("fill_value", "m/year", "m/s"));

  output = result;
}

PSB_velsurf_mag::PSB_velsurf_mag(StressBalance *m)
  : Diag<StressBalance>(m) {
  // set metadata:
  m_vars.push_back(NCSpatialVariable(m_grid.config.get_unit_system(), "velsurf_mag", m_grid));

  set_attrs("magnitude of horizontal velocity of ice at ice surface", "",
            "m s-1", "m year-1", 0);
  m_vars[0].set_double("_FillValue", m_grid.config.get("fill_value", "m/year", "m/s"));
  m_vars[0].set_double("valid_min",  0.0);
}

void PSB_velsurf_mag::compute(IceModelVec* &output) {

  IceModelVec3 *u3, *v3, *w3;
  IceModelVec2S tmp, *result;

  tmp.create(m_grid, "tmp", WITHOUT_GHOSTS);

  result = new IceModelVec2S;
  result->create(m_grid, "velsurf_mag", WITHOUT_GHOSTS);
  result->metadata() = m_vars[0];

  model->get_3d_velocity(u3, v3, w3);

  const IceModelVec2S *thickness = m_grid.variables().get_2d_scalar("land_ice_thickness");

  u3->getSurfaceValues(*result, *thickness);
  v3->getSurfaceValues(tmp, *thickness);

  result->set_to_magnitude(*result, tmp);

  // mask out ice-free areas
  result->mask_by(*thickness, m_grid.config.get("fill_value", "m/year", "m/s"));

  output = result;
}


PSB_velsurf::PSB_velsurf(StressBalance *m)
  : Diag<StressBalance>(m) {

  // set metadata:
  m_dof = 2;

  m_vars.push_back(NCSpatialVariable(m_grid.config.get_unit_system(), "uvelsurf", m_grid));
  m_vars.push_back(NCSpatialVariable(m_grid.config.get_unit_system(), "vvelsurf", m_grid));

  set_attrs("x-component of the horizontal velocity of ice at ice surface", "",
            "m s-1", "m year-1", 0);
  set_attrs("y-component of the horizontal velocity of ice at ice surface", "",
            "m s-1", "m year-1", 1);

  m_vars[0].set_double("valid_min", m_grid.convert(-1e6, "m/year", "m/second"));
  m_vars[0].set_double("valid_max", m_grid.convert(1e6, "m/year", "m/second"));
  m_vars[0].set_double("_FillValue", m_grid.config.get("fill_value", "m/year", "m/s"));

  m_vars[1].set_double("valid_min", m_grid.convert(-1e6, "m/year", "m/second"));
  m_vars[1].set_double("valid_max", m_grid.convert(1e6, "m/year", "m/second"));
  m_vars[1].set_double("_FillValue", m_grid.config.get("fill_value", "m/year", "m/s"));
}

void PSB_velsurf::compute(IceModelVec* &output) {
  IceModelVec2V *result;
  IceModelVec3 *u3, *v3, *w3;
  IceModelVec2S tmp;
  double fill_value = m_grid.config.get("fill_value", "m/year", "m/s");

  result = new IceModelVec2V;
  result->create(m_grid, "surf", WITHOUT_GHOSTS);
  result->metadata(0) = m_vars[0];
  result->metadata(1) = m_vars[1];

  tmp.create(m_grid, "tmp", WITHOUT_GHOSTS);

  model->get_3d_velocity(u3, v3, w3);

  const IceModelVec2S *thickness = m_grid.variables().get_2d_scalar("land_ice_thickness");

  u3->getSurfaceValues(tmp, *thickness);
  result->set_component(0, tmp);

  v3->getSurfaceValues(tmp, *thickness);
  result->set_component(1, tmp);

  const IceModelVec2Int *mask = m_grid.variables().get_2d_mask("mask");

  MaskQuery M(*mask);

  IceModelVec::AccessList list;
  list.add(*mask);
  list.add(*result);

  for (Points p(m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (M.ice_free(i, j)) {
      (*result)(i, j).u = fill_value;
      (*result)(i, j).v = fill_value;
    }
  }


  output = result;
}

PSB_wvel::PSB_wvel(StressBalance *m)
  : Diag<StressBalance>(m) {

  // set metadata:
  m_vars.push_back(NCSpatialVariable(m_grid.config.get_unit_system(), "wvel", m_grid,
                                     m_grid.z()));

  set_attrs("vertical velocity of ice, relative to geoid", "",
            "m s-1", "m year-1", 0);
  m_vars[0].set_double("valid_min", m_grid.convert(-1e6, "m/year", "m/second"));
  m_vars[0].set_double("valid_max", m_grid.convert(1e6, "m/year", "m/second"));
}

void PSB_wvel::compute(IceModelVec* &output) {
  IceModelVec3 *result, *u3, *v3, *w3;
  double *u, *v, *w, *res;

  result = new IceModelVec3;
  result->create(m_grid, "wvel", WITHOUT_GHOSTS);
  result->metadata() = m_vars[0];

  const IceModelVec2S *bed, *uplift;
  bed       = m_grid.variables().get_2d_scalar("bedrock_altitude");
  uplift    = m_grid.variables().get_2d_scalar("tendency_of_bedrock_altitude");

  const IceModelVec2S   *thickness = m_grid.variables().get_2d_scalar("land_ice_thickness");
  const IceModelVec2Int *mask      = m_grid.variables().get_2d_mask("mask");

  model->get_3d_velocity(u3, v3, w3);

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

  const double ice_density = m_grid.config.get("ice_density"),
    sea_water_density = m_grid.config.get("sea_water_density"),
    R = ice_density / sea_water_density;

  for (Points p(m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    u3->getInternalColumn(i, j, &u);
    v3->getInternalColumn(i, j, &v);
    w3->getInternalColumn(i, j, &w);
    result->getInternalColumn(i, j, &res);

    int ks = m_grid.kBelowHeight((*thickness)(i,j));

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
    for (unsigned int k = ks+1; k < m_grid.Mz() ; k++) {
      res[k] = 0.0;
    }

  }

  output = result;
}

PSB_wvelsurf::PSB_wvelsurf(StressBalance *m)
  : Diag<StressBalance>(m) {

  // set metadata:
  m_vars.push_back(NCSpatialVariable(m_grid.config.get_unit_system(), "wvelsurf", m_grid));

  set_attrs("vertical velocity of ice at ice surface, relative to the geoid", "",
            "m s-1", "m year-1", 0);
  m_vars[0].set_double("valid_min", m_grid.convert(-1e6, "m/year", "m/second"));
  m_vars[0].set_double("valid_max", m_grid.convert(1e6, "m/year", "m/second"));
  m_vars[0].set_double("_FillValue", m_grid.config.get("fill_value", "m/year", "m/s"));
}

void PSB_wvelsurf::compute(IceModelVec* &output) {
  IceModelVec *tmp;
  IceModelVec3 *w3;
  IceModelVec2S *result;
  double fill_value = m_grid.config.get("fill_value", "m/year", "m/s");

  result = new IceModelVec2S;
  result->create(m_grid, "wvelsurf", WITHOUT_GHOSTS);
  result->metadata() = m_vars[0];

  PSB_wvel wvel(model);

  wvel.compute(tmp);

  w3 = dynamic_cast<IceModelVec3*>(tmp);
  if (tmp == NULL) {
    throw RuntimeError("dynamic_cast failure");
  }

  const IceModelVec2S *thickness = m_grid.variables().get_2d_scalar("land_ice_thickness");

  w3->getSurfaceValues(*result, *thickness);

  const IceModelVec2Int *mask = m_grid.variables().get_2d_mask("mask");

  MaskQuery M(*mask);

  IceModelVec::AccessList list;
  list.add(*mask);
  list.add(*result);

  for (Points p(m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (M.ice_free(i, j)) {
      (*result)(i, j) = fill_value;
    }
  }


  delete tmp;
  output = result;
}

PSB_wvelbase::PSB_wvelbase(StressBalance *m)
  : Diag<StressBalance>(m) {

  // set metadata:
  m_vars.push_back(NCSpatialVariable(m_grid.config.get_unit_system(), "wvelbase", m_grid));

  set_attrs("vertical velocity of ice at the base of ice, relative to the geoid", "",
            "m s-1", "m year-1", 0);
  m_vars[0].set_double("valid_min", m_grid.convert(-1e6, "m/year", "m/second"));
  m_vars[0].set_double("valid_max", m_grid.convert(1e6, "m/year", "m/second"));
  m_vars[0].set_double("_FillValue", m_grid.config.get("fill_value", "m/year", "m/s"));
}

void PSB_wvelbase::compute(IceModelVec* &output) {
  IceModelVec *tmp;
  IceModelVec3 *w3;
  IceModelVec2S *result;
  double fill_value = m_grid.config.get("fill_value", "m/year", "m/s");

  result = new IceModelVec2S;
  result->create(m_grid, "wvelbase", WITHOUT_GHOSTS);
  result->metadata() = m_vars[0];

  PSB_wvel wvel(model);

  wvel.compute(tmp);

  w3 = dynamic_cast<IceModelVec3*>(tmp);
  if (tmp == NULL) {
    throw RuntimeError("dynamic_cast failure");
  }

  w3->getHorSlice(*result, 0.0);

  const IceModelVec2Int *mask = m_grid.variables().get_2d_mask("mask");

  MaskQuery M(*mask);

  IceModelVec::AccessList list;
  list.add(*mask);
  list.add(*result);

  for (Points p(m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (M.ice_free(i, j)) {
      (*result)(i, j) = fill_value;
    }
  }


  delete tmp;
  output = result;
}

PSB_velbase::PSB_velbase(StressBalance *m)
  : Diag<StressBalance>(m) {

  // set metadata:
  m_dof = 2;

  m_vars.push_back(NCSpatialVariable(m_grid.config.get_unit_system(), "uvelbase", m_grid));
  m_vars.push_back(NCSpatialVariable(m_grid.config.get_unit_system(), "vvelbase", m_grid));

  set_attrs("x-component of the horizontal velocity of ice at the base of ice", "",
            "m s-1", "m year-1", 0);
  set_attrs("y-component of the horizontal velocity of ice at the base of ice", "",
            "m s-1", "m year-1", 1);

  m_vars[0].set_double("valid_min", m_grid.convert(-1e6, "m/year", "m/second"));
  m_vars[0].set_double("valid_max", m_grid.convert(1e6, "m/year", "m/second"));
  m_vars[0].set_double("_FillValue", m_grid.config.get("fill_value", "m/year", "m/s"));

  m_vars[1].set_double("valid_min", m_grid.convert(-1e6, "m/year", "m/second"));
  m_vars[1].set_double("valid_max", m_grid.convert(1e6, "m/year", "m/second"));
  m_vars[1].set_double("_FillValue", m_grid.config.get("fill_value", "m/year", "m/s"));
}

void PSB_velbase::compute(IceModelVec* &output) {
  IceModelVec2V *result;
  IceModelVec3 *u3, *v3, *w3;
  IceModelVec2S tmp;            // will be de-allocated automatically
  double fill_value = m_grid.config.get("fill_value", "m/year", "m/s");

  result = new IceModelVec2V;
  result->create(m_grid, "base", WITHOUT_GHOSTS);
  result->metadata() = m_vars[0];
  result->metadata(1) = m_vars[1];

  tmp.create(m_grid, "tmp", WITHOUT_GHOSTS);

  model->get_3d_velocity(u3, v3, w3);

  u3->getHorSlice(tmp, 0.0);
  result->set_component(0, tmp);

  v3->getHorSlice(tmp, 0.0);
  result->set_component(1, tmp);

  const IceModelVec2Int *mask = m_grid.variables().get_2d_mask("mask");

  MaskQuery M(*mask);

  IceModelVec::AccessList list;
  list.add(*mask);
  list.add(*result);

  for (Points p(m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (M.ice_free(i, j)) {
      (*result)(i, j).u = fill_value;
      (*result)(i, j).v = fill_value;
    }
  }


  output = result;
}


PSB_bfrict::PSB_bfrict(StressBalance *m)
  : Diag<StressBalance>(m) {

  // set metadata:
  m_vars.push_back(NCSpatialVariable(m_grid.config.get_unit_system(), "bfrict", m_grid));

  set_attrs("basal frictional heating", "",
            "W m-2", "W m-2", 0);
}

void PSB_bfrict::compute(IceModelVec* &output) {

  IceModelVec2S *result = new IceModelVec2S;
  result->create(m_grid, "bfrict", WITHOUT_GHOSTS);
  result->metadata() = m_vars[0];

  IceModelVec2S *bfrict;
  model->get_basal_frictional_heating(bfrict);

  bfrict->copy_to(*result);

  output = result;
}


PSB_uvel::PSB_uvel(StressBalance *m)
  : Diag<StressBalance>(m) {

  // set metadata:
  m_vars.push_back(NCSpatialVariable(m_grid.config.get_unit_system(), "uvel", m_grid, m_grid.z()));

  set_attrs("horizontal velocity of ice in the X direction", "land_ice_x_velocity",
            "m s-1", "m year-1", 0);
}

void PSB_uvel::compute(IceModelVec* &output) {

  IceModelVec3 *result = new IceModelVec3;
  result->create(m_grid, "uvel", WITHOUT_GHOSTS);
  result->metadata() = m_vars[0];

  const IceModelVec2S *thickness = m_grid.variables().get_2d_scalar("land_ice_thickness");

  IceModelVec3 *u3, *v3, *w3;
  model->get_3d_velocity(u3, v3, w3);

  IceModelVec::AccessList list;
  list.add(*u3);
  list.add(*result);
  list.add(*thickness);

  double *u_ij, *u_out_ij;
  for (Points p(m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    int ks = m_grid.kBelowHeight((*thickness)(i,j));

    u3->getInternalColumn(i,j,&u_ij);
    result->getInternalColumn(i,j,&u_out_ij);

    // in the ice:
    for (int k = 0; k <= ks ; k++) {
      u_out_ij[k] = u_ij[k];
    }
    // above the ice:
    for (unsigned int k = ks+1; k < m_grid.Mz() ; k++) {
      u_out_ij[k] = 0.0;
    }
  }


  output = result;
}

PSB_vvel::PSB_vvel(StressBalance *m)
  : Diag<StressBalance>(m) {

  // set metadata:
  m_vars.push_back(NCSpatialVariable(m_grid.config.get_unit_system(), "vvel", m_grid, m_grid.z()));

  set_attrs("horizontal velocity of ice in the Y direction", "land_ice_y_velocity",
            "m s-1", "m year-1", 0);
}

void PSB_vvel::compute(IceModelVec* &output) {

  IceModelVec3 *result = new IceModelVec3;
  result->create(m_grid, "vvel", WITHOUT_GHOSTS);
  result->metadata() = m_vars[0];

  const IceModelVec2S *thickness = m_grid.variables().get_2d_scalar("land_ice_thickness");

  IceModelVec3 *u3, *v3, *w3;
  model->get_3d_velocity(u3, v3, w3);

  IceModelVec::AccessList list;
  list.add(*v3);
  list.add(*result);
  list.add(*thickness);

  double *v_ij, *v_out_ij;
  for (Points p(m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    int ks = m_grid.kBelowHeight((*thickness)(i,j));

    v3->getInternalColumn(i,j,&v_ij);
    result->getInternalColumn(i,j,&v_out_ij);

    // in the ice:
    for (int k = 0; k <= ks ; k++) {
      v_out_ij[k] = v_ij[k];
    }
    // above the ice:
    for (unsigned int k = ks+1; k < m_grid.Mz() ; k++) {
      v_out_ij[k] = 0.0;
    }
  }


  output = result;
}

PSB_wvel_rel::PSB_wvel_rel(StressBalance *m)
  : Diag<StressBalance>(m) {

  // set metadata:
  m_vars.push_back(NCSpatialVariable(m_grid.config.get_unit_system(), "wvel_rel", m_grid, m_grid.z()));

  set_attrs("vertical velocity of ice, relative to base of ice directly below", "",
            "m s-1", "m year-1", 0);
}

void PSB_wvel_rel::compute(IceModelVec* &output) {

  IceModelVec3 *result = new IceModelVec3;
  result->create(m_grid, "wvel_rel", WITHOUT_GHOSTS);
  result->metadata() = m_vars[0];

  const IceModelVec2S *thickness = m_grid.variables().get_2d_scalar("land_ice_thickness");

  IceModelVec3 *u3, *v3, *w3;
  model->get_3d_velocity(u3, v3, w3);

  IceModelVec::AccessList list;
  list.add(*w3);
  list.add(*result);
  list.add(*thickness);

  double *w_ij, *w_out_ij;
  for (Points p(m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    int ks = m_grid.kBelowHeight((*thickness)(i,j));

    w3->getInternalColumn(i,j,&w_ij);
    result->getInternalColumn(i,j,&w_out_ij);

    // in the ice:
    for (int k = 0; k <= ks ; k++) {
      w_out_ij[k] = w_ij[k];
    }
    // above the ice:
    for (unsigned int k = ks+1; k < m_grid.Mz() ; k++) {
      w_out_ij[k] = 0.0;
    }
  }


  output = result;
}


PSB_strainheat::PSB_strainheat(StressBalance *m)
  : Diag<StressBalance>(m) {

  // set metadata:
  m_vars.push_back(NCSpatialVariable(m_grid.config.get_unit_system(), "strainheat", m_grid, m_grid.z()));

  set_attrs("rate of strain heating in ice (dissipation heating)", "",
            "W m-3", "mW m-3", 0);
}

void PSB_strainheat::compute(IceModelVec* &output) {

  IceModelVec3 *result = new IceModelVec3;
  result->create(m_grid, "strainheat", WITHOUT_GHOSTS);
  result->metadata() = m_vars[0];
  result->write_in_glaciological_units = true;

  IceModelVec3 *tmp;
  model->get_volumetric_strain_heating(tmp);

  tmp->copy_to(*result);

  output = result;
}

PSB_strain_rates::PSB_strain_rates(StressBalance *m)
  : Diag<StressBalance>(m) {
  m_dof = 2;

  // set metadata:
  m_vars.push_back(NCSpatialVariable(m_grid.config.get_unit_system(), "eigen1", m_grid));
  m_vars.push_back(NCSpatialVariable(m_grid.config.get_unit_system(), "eigen2", m_grid));

  set_attrs("first eigenvalue of the horizontal, vertically-integrated strain rate tensor",
            "", "s-1", "s-1", 0);
  set_attrs("second eigenvalue of the horizontal, vertically-integrated strain rate tensor",
            "", "s-1", "s-1", 1);
}

void PSB_strain_rates::compute(IceModelVec* &output) {
  IceModelVec2 *result;
  IceModelVec *velbar;
  PSB_velbar diag(model);

  result = new IceModelVec2;
  result->create(m_grid, "strain_rates", WITHOUT_GHOSTS, 1, 2);
  result->metadata() = m_vars[0];
  result->metadata(1) = m_vars[1];

  const IceModelVec2Int *mask = m_grid.variables().get_2d_mask("mask");

  diag.compute(velbar);
  IceModelVec2V *v_tmp = dynamic_cast<IceModelVec2V*>(velbar);
  if (v_tmp == NULL) {
    throw RuntimeError("velbar is expected to be an IceModelVec2V");
  }

  IceModelVec2V velbar_with_ghosts;
  velbar_with_ghosts.create(m_grid, "velbar", WITH_GHOSTS);

  // copy_from communicates ghosts
  velbar_with_ghosts.copy_from(*v_tmp);

  model->compute_2D_principal_strain_rates(velbar_with_ghosts, *mask, *result);

  delete velbar;
  output = result;
}

PSB_deviatoric_stresses::PSB_deviatoric_stresses(StressBalance *m)
  : Diag<StressBalance>(m) {
  m_dof = 3;

  // set metadata:
  m_vars.push_back(NCSpatialVariable(m_grid.config.get_unit_system(), "sigma_xx", m_grid));
  m_vars.push_back(NCSpatialVariable(m_grid.config.get_unit_system(), "sigma_yy", m_grid));
  m_vars.push_back(NCSpatialVariable(m_grid.config.get_unit_system(), "sigma_xy", m_grid));

  set_attrs("deviatoric stress in X direction", "", "Pa", "Pa", 0);
  set_attrs("deviatoric stress in Y direction", "", "Pa", "Pa", 1);
  set_attrs("deviatoric shear stress", "", "Pa", "Pa", 2);

}

void PSB_deviatoric_stresses::compute(IceModelVec* &output) {
  IceModelVec2 *result;
  IceModelVec *velbar;
  PSB_velbar diag(model);

  result = new IceModelVec2;
  result->create(m_grid, "deviatoric_stresses", WITHOUT_GHOSTS, 1, 3);
  result->metadata() = m_vars[0];
  result->metadata(1) = m_vars[1];
  result->metadata(2) = m_vars[2];

  const IceModelVec2Int *mask = m_grid.variables().get_2d_mask("mask");

  diag.compute(velbar);
  IceModelVec2V *v_tmp = dynamic_cast<IceModelVec2V*>(velbar);
  if (v_tmp == NULL) {
    throw RuntimeError("velbar is expected to be an IceModelVec2V");
  }

  IceModelVec2V velbar_with_ghosts;
  velbar_with_ghosts.create(m_grid, "velbar", WITH_GHOSTS);

  // copy_from communicates ghosts
  velbar_with_ghosts.copy_from(*v_tmp);

  model->compute_2D_stresses(velbar_with_ghosts, *mask, *result);

  output = result;
}

PSB_pressure::PSB_pressure(StressBalance *m)
  : Diag<StressBalance>(m) {

  // set metadata:
  m_vars.push_back(NCSpatialVariable(m_grid.config.get_unit_system(), "pressure", m_grid, m_grid.z()));

  set_attrs("pressure in ice (hydrostatic)", "",
            "Pa", "Pa", 0);
}

void PSB_pressure::compute(IceModelVec* &output) {

  IceModelVec3 *result = new IceModelVec3;
  result->create(m_grid, "pressure", WITHOUT_GHOSTS);
  result->metadata() = m_vars[0];

  const IceModelVec2S *thickness = m_grid.variables().get_2d_scalar("land_ice_thickness");

  IceModelVec::AccessList list;
  list.add(*result);
  list.add(*thickness);

  double *P_out_ij;
  const double rg = m_grid.config.get("ice_density") * m_grid.config.get("standard_gravity");

  for (Points p(m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    unsigned int ks = m_grid.kBelowHeight((*thickness)(i,j));
    result->getInternalColumn(i,j,&P_out_ij);
    const double H = (*thickness)(i,j);
    // within the ice:
    for (unsigned int k = 0; k <= ks; ++k) {
      P_out_ij[k] = rg * (H - m_grid.z(k));  // FIXME: add atmospheric pressure?
    }
    // above the ice:
    for (unsigned int k = ks + 1; k < m_grid.Mz(); ++k) {
      P_out_ij[k] = 0.0;  // FIXME: use atmospheric pressure?
    }
  }


  output = result;
}


PSB_tauxz::PSB_tauxz(StressBalance *m)
  : Diag<StressBalance>(m) {

  // set metadata:
  m_vars.push_back(NCSpatialVariable(m_grid.config.get_unit_system(), "tauxz", m_grid, m_grid.z()));

  set_attrs("shear stress xz component (in shallow ice approximation SIA)", "",
            "Pa", "Pa", 0);
}


/*!
 * The SIA-applicable shear stress component tauxz computed here is not used
 * by the model.  This implementation intentionally does not use the
 * eta-transformation or special cases at ice margins.
 * CODE DUPLICATION WITH PSB_tauyz
 */
void PSB_tauxz::compute(IceModelVec* &output) {

  IceModelVec3 *result = new IceModelVec3;
  result->create(m_grid, "tauxz", WITHOUT_GHOSTS);
  result->metadata() = m_vars[0];

  const IceModelVec2S *thickness, *surface;

  thickness = m_grid.variables().get_2d_scalar("land_ice_thickness");
  surface = m_grid.variables().get_2d_scalar("surface_altitude");

  IceModelVec::AccessList list;
  list.add(*result);
  list.add(*surface);
  list.add(*thickness);

  double *tauxz_out_ij;
  const double rg = m_grid.config.get("ice_density") * m_grid.config.get("standard_gravity");

  for (Points p(m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();


    unsigned int ks = m_grid.kBelowHeight((*thickness)(i,j));
    result->getInternalColumn(i,j,&tauxz_out_ij);
    const double H    = (*thickness)(i,j),
      dhdx = surface->diff_x_p(i,j);
    // within the ice:
    for (unsigned int k = 0; k <= ks; ++k) {
      tauxz_out_ij[k] = - rg * (H - m_grid.z(k)) * dhdx;
    }
    // above the ice:
    for (unsigned int k = ks + 1; k < m_grid.Mz(); ++k) {
      tauxz_out_ij[k] = 0.0;
    }

  }


  output = result;
}


PSB_tauyz::PSB_tauyz(StressBalance *m)
  : Diag<StressBalance>(m) {

  // set metadata:
  m_vars.push_back(NCSpatialVariable(m_grid.config.get_unit_system(), "tauyz", m_grid, m_grid.z()));

  set_attrs("shear stress yz component (in shallow ice approximation SIA)", "",
            "Pa", "Pa", 0);
}


/*!
 * The SIA-applicable shear stress component tauyz computed here is not used
 * by the model.  This implementation intentionally does not use the
 * eta-transformation or special cases at ice margins.
 * CODE DUPLICATION WITH PSB_tauxz
 */
void PSB_tauyz::compute(IceModelVec* &output) {

  IceModelVec3 *result = new IceModelVec3;
  result->create(m_grid, "tauyz", WITHOUT_GHOSTS);
  result->metadata() = m_vars[0];

  const IceModelVec2S *thickness = m_grid.variables().get_2d_scalar("land_ice_thickness");
  const IceModelVec2S *surface = m_grid.variables().get_2d_scalar("surface_altitude");

  IceModelVec::AccessList list;
  list.add(*result);
  list.add(*surface);
  list.add(*thickness);

  double *tauyz_out_ij;
  const double rg = m_grid.config.get("ice_density") * m_grid.config.get("standard_gravity");

  for (Points p(m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();


    unsigned int ks = m_grid.kBelowHeight((*thickness)(i,j));
    result->getInternalColumn(i,j,&tauyz_out_ij);
    const double H    = (*thickness)(i,j),
      dhdy = surface->diff_y_p(i,j);
    // within the ice:
    for (unsigned int k = 0; k <= ks; ++k) {
      tauyz_out_ij[k] = - rg * (H - m_grid.z(k)) * dhdy;
    }
    // above the ice:
    for (unsigned int k = ks + 1; k < m_grid.Mz(); ++k) {
      tauyz_out_ij[k] = 0.0;
    }

  }


  output = result;
}


} // end of namespace pism
