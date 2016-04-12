// Copyright (C) 2010, 2011, 2012, 2013, 2014, 2015, 2016 Constantine Khroulev
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
#include "SSB_Modifier.hh"
#include "ShallowStressBalance.hh"
#include "base/util/Mask.hh"
#include "base/util/PISMConfigInterface.hh"
#include "base/util/PISMVars.hh"
#include "base/util/error_handling.hh"
#include "base/util/IceModelVec2CellType.hh"

namespace pism {
namespace stressbalance {

using units::convert;

void StressBalance::get_diagnostics_impl(std::map<std::string, Diagnostic::Ptr> &dict,
                                        std::map<std::string, TSDiagnostic::Ptr> &ts_dict) {

  dict["bfrict"]   = Diagnostic::Ptr(new PSB_bfrict(this));

  dict["velbar_mag"]     = Diagnostic::Ptr(new PSB_velbar_mag(this));
  dict["flux"]           = Diagnostic::Ptr(new PSB_flux(this));
  dict["flux_mag"]       = Diagnostic::Ptr(new PSB_flux_mag(this));
  dict["velbase_mag"]    = Diagnostic::Ptr(new PSB_velbase_mag(this));
  dict["velsurf_mag"]    = Diagnostic::Ptr(new PSB_velsurf_mag(this));

  dict["uvel"]     = Diagnostic::Ptr(new PSB_uvel(this));
  dict["vvel"]     = Diagnostic::Ptr(new PSB_vvel(this));

  dict["strainheat"] = Diagnostic::Ptr(new PSB_strainheat(this));

  dict["velbar"]   = Diagnostic::Ptr(new PSB_velbar(this));
  dict["velbase"]  = Diagnostic::Ptr(new PSB_velbase(this));
  dict["velsurf"]  = Diagnostic::Ptr(new PSB_velsurf(this));

  dict["wvel"]     = Diagnostic::Ptr(new PSB_wvel(this));
  dict["wvelbase"] = Diagnostic::Ptr(new PSB_wvelbase(this));
  dict["wvelsurf"] = Diagnostic::Ptr(new PSB_wvelsurf(this));
  dict["wvel_rel"] = Diagnostic::Ptr(new PSB_wvel_rel(this));
  dict["strain_rates"] = Diagnostic::Ptr(new PSB_strain_rates(this));
  dict["deviatoric_stresses"] = Diagnostic::Ptr(new PSB_deviatoric_stresses(this));

  dict["pressure"] = Diagnostic::Ptr(new PSB_pressure(this));
  dict["tauxz"] = Diagnostic::Ptr(new PSB_tauxz(this));
  dict["tauyz"] = Diagnostic::Ptr(new PSB_tauyz(this));

  m_shallow_stress_balance->get_diagnostics(dict, ts_dict);
  m_modifier->get_diagnostics(dict, ts_dict);
}

PSB_velbar::PSB_velbar(StressBalance *m)
  : Diag<StressBalance>(m) {

  m_dof = 2;

  // set metadata:
  m_vars.push_back(SpatialVariableMetadata(m_sys, "ubar"));
  m_vars.push_back(SpatialVariableMetadata(m_sys, "vbar"));

  set_attrs("vertical mean of horizontal ice velocity in the X direction",
            "land_ice_vertical_mean_x_velocity",
            "m s-1", "m year-1", 0);
  set_attrs("vertical mean of horizontal ice velocity in the Y direction",
            "land_ice_vertical_mean_y_velocity",
            "m s-1", "m year-1", 1);
}

IceModelVec::Ptr PSB_velbar::compute_impl() {
  // get the thickness
  const IceModelVec2S* thickness = m_grid->variables().get_2d_scalar("land_ice_thickness");

  // Compute the vertically-integrated horizontal ice flux:
  IceModelVec2V::Ptr result = IceModelVec2V::ToVector(PSB_flux(model).compute());

  // Override metadata set by the flux computation
  result->metadata(0) = m_vars[0];
  result->metadata(1) = m_vars[1];

  IceModelVec::AccessList list;
  list.add(*thickness);
  list.add(*result);

  for (Points p(*m_grid); p; p.next()) {
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

  return result;
}

PSB_velbar_mag::PSB_velbar_mag(StressBalance *m)
  : Diag<StressBalance>(m) {

  // set metadata:
  m_vars.push_back(SpatialVariableMetadata(m_sys, "velbar_mag"));

  set_attrs("magnitude of vertically-integrated horizontal velocity of ice", "",
            "m s-1", "m year-1", 0);
  m_vars[0].set_double("_FillValue", convert(m_sys, m_config->get_double("fill_value"),
                                         "m year-1", "m second-1"));
  m_vars[0].set_double("valid_min", 0.0);
}

IceModelVec::Ptr PSB_velbar_mag::compute_impl() {

  IceModelVec2S::Ptr result(new IceModelVec2S);
  result->create(m_grid, "velbar_mag", WITHOUT_GHOSTS);
  result->metadata(0) = m_vars[0];

  // compute vertically-averaged horizontal velocity:
  IceModelVec2V::Ptr velbar = IceModelVec2V::ToVector(PSB_velbar(model).compute());

  // compute its magnitude:
  result->set_to_magnitude(*velbar);

  const IceModelVec2S *thickness = m_grid->variables().get_2d_scalar("land_ice_thickness");

  // mask out ice-free areas:
  double fill_value = convert(m_sys, m_config->get_double("fill_value"), "m year-1", "m second-1");
  result->mask_by(*thickness, fill_value);

  return result;
}


PSB_flux::PSB_flux(StressBalance *m)
  : Diag<StressBalance>(m) {

  m_dof = 2;

  // set metadata:
  m_vars.push_back(SpatialVariableMetadata(m_sys, "uflux"));
  m_vars.push_back(SpatialVariableMetadata(m_sys, "vflux"));

  set_attrs("Vertically integrated horizontal flux of ice in the X direction",
            "",                 // no CF standard name
            "m2 s-1", "m2 year-1", 0);
  set_attrs("Vertically integrated horizontal flux of ice in the Y direction",
            "",                 // no CF standard name
            "m2 s-1", "m2 year-1", 1);
}

IceModelVec::Ptr PSB_flux::compute_impl() {
  double icefree_thickness = m_config->get_double("mask_icefree_thickness_standard");

  IceModelVec2V::Ptr result(new IceModelVec2V);
  result->create(m_grid, "flux", WITHOUT_GHOSTS);
  result->metadata(0) = m_vars[0];
  result->metadata(1) = m_vars[1];

  // get the thickness
  const IceModelVec2S *thickness = m_grid->variables().get_2d_scalar("land_ice_thickness");

  const IceModelVec3
    &u3 = model->velocity_u(),
    &v3 = model->velocity_v();

  IceModelVec::AccessList list;
  list.add(u3);
  list.add(v3);
  list.add(*thickness);
  list.add(*result);

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      double u_sum = 0, v_sum = 0,
        thk = (*thickness)(i,j);
      int ks = m_grid->kBelowHeight(thk);

      // an "ice-free" cell:
      if (thk < icefree_thickness) {
        (*result)(i,j).u = 0;
        (*result)(i,j).v = 0;
        continue;
      }

      // an ice-filled cell:
      const double *u_ij = NULL, *v_ij = NULL;
      u_ij = u3.get_column(i, j);
      v_ij = v3.get_column(i, j);

      if (thk <= m_grid->z(1)) {
        (*result)(i,j).u = u_ij[0];
        (*result)(i,j).v = v_ij[0];
        continue;
      }

      for (int k = 1; k <= ks; ++k) {
        u_sum += (m_grid->z(k) - m_grid->z(k-1)) * (u_ij[k] + u_ij[k-1]);
        v_sum += (m_grid->z(k) - m_grid->z(k-1)) * (v_ij[k] + v_ij[k-1]);
      }

      // Finish the trapezoidal rule integration (multiply by 1/2).
      (*result)(i,j).u = 0.5 * u_sum;
      (*result)(i,j).v = 0.5 * v_sum;

      // The top surface of the ice is not aligned with the grid, so
      // we have at most dz meters of ice above grid.z(ks).
      // Assume that its velocity is (u_ij[ks], v_ij[ks]) and add its
      // contribution.
      (*result)(i,j).u += u_ij[ks] * (thk - m_grid->z(ks));
      (*result)(i,j).v += v_ij[ks] * (thk - m_grid->z(ks));
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  return result;
}


PSB_flux_mag::PSB_flux_mag(StressBalance *m)
  : Diag<StressBalance>(m) {

  // set metadata:
  m_vars.push_back(SpatialVariableMetadata(m_sys, "flux_mag"));

  set_attrs("magnitude of vertically-integrated horizontal flux of ice", "",
            "m2 s-1", "m2 year-1", 0);

  double fill_value = convert(m_sys, m_config->get_double("fill_value"),
                              "m2 year-1", "m2 second-1");
  m_vars[0].set_double("_FillValue", fill_value);
  m_vars[0].set_double("valid_min", 0.0);
}

IceModelVec::Ptr PSB_flux_mag::compute_impl() {
  const IceModelVec2S *thickness = m_grid->variables().get_2d_scalar("land_ice_thickness");

  // Compute the vertically-averaged horizontal ice velocity:
  IceModelVec2S::Ptr result = IceModelVec2S::To2DScalar(PSB_velbar_mag(model).compute());

  IceModelVec::AccessList list;
  list.add(*thickness);
  list.add(*result);

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    (*result)(i,j) *= (*thickness)(i,j);
  }

  double fill_value = convert(m_sys, m_config->get_double("fill_value"), "m year-1", "m second-1");
  result->mask_by(*thickness, fill_value);

  result->metadata() = m_vars[0];

  return result;
}

PSB_velbase_mag::PSB_velbase_mag(StressBalance *m)
  : Diag<StressBalance>(m) {

  // set metadata:
  m_vars.push_back(SpatialVariableMetadata(m_sys, "velbase_mag"));

  set_attrs("magnitude of horizontal velocity of ice at base of ice", "",
            "m s-1", "m year-1", 0);

  double fill_value = convert(m_sys, m_config->get_double("fill_value"), "m year-1", "m second-1");

  m_vars[0].set_double("_FillValue", fill_value);
  m_vars[0].set_double("valid_min", 0.0);
}

IceModelVec::Ptr PSB_velbase_mag::compute_impl() {
  // FIXME: compute this using PSB_velbase.

  IceModelVec2S tmp;
  tmp.create(m_grid, "tmp", WITHOUT_GHOSTS);

  IceModelVec2S::Ptr result(new IceModelVec2S);
  result->create(m_grid, "velbase_mag", WITHOUT_GHOSTS);
  result->metadata(0) = m_vars[0];

  const IceModelVec3
    &u3 = model->velocity_u(),
    &v3 = model->velocity_v();

  const IceModelVec2S *thickness = m_grid->variables().get_2d_scalar("land_ice_thickness");

  u3.getHorSlice(*result, 0.0); // result = u_{z=0}
  v3.getHorSlice(tmp, 0.0);    // tmp = v_{z=0}

  result->set_to_magnitude(*result, tmp);

  // mask out ice-free areas
  double fill_value = convert(m_sys, m_config->get_double("fill_value"), "m year-1", "m second-1");
  result->mask_by(*thickness, fill_value);

  return result;
}

PSB_velsurf_mag::PSB_velsurf_mag(StressBalance *m)
  : Diag<StressBalance>(m) {
  // set metadata:
  m_vars.push_back(SpatialVariableMetadata(m_sys, "velsurf_mag"));

  set_attrs("magnitude of horizontal velocity of ice at ice surface", "",
            "m s-1", "m year-1", 0);

  double fill_value = convert(m_sys, m_config->get_double("fill_value"), "m year-1", "m second-1");

  m_vars[0].set_double("_FillValue", fill_value);
  m_vars[0].set_double("valid_min",  0.0);
}

IceModelVec::Ptr PSB_velsurf_mag::compute_impl() {

  // FIXME: Compute this using PSB_velsurf.

  IceModelVec2S tmp;
  tmp.create(m_grid, "tmp", WITHOUT_GHOSTS);

  IceModelVec2S::Ptr result(new IceModelVec2S);
  result->create(m_grid, "velsurf_mag", WITHOUT_GHOSTS);
  result->metadata(0) = m_vars[0];

  const IceModelVec3
    &u3 = model->velocity_u(),
    &v3 = model->velocity_v();

  const IceModelVec2S *thickness = m_grid->variables().get_2d_scalar("land_ice_thickness");

  u3.getSurfaceValues(*result, *thickness);
  v3.getSurfaceValues(tmp, *thickness);

  result->set_to_magnitude(*result, tmp);

  // mask out ice-free areas
  double fill_value = convert(m_sys, m_config->get_double("fill_value"), "m year-1", "m second-1");
  result->mask_by(*thickness, fill_value);

  return result;
}


PSB_velsurf::PSB_velsurf(StressBalance *m)
  : Diag<StressBalance>(m) {

  // set metadata:
  m_dof = 2;

  m_vars.push_back(SpatialVariableMetadata(m_sys, "uvelsurf"));
  m_vars.push_back(SpatialVariableMetadata(m_sys, "vvelsurf"));

  set_attrs("x-component of the horizontal velocity of ice at ice surface", "",
            "m s-1", "m year-1", 0);
  set_attrs("y-component of the horizontal velocity of ice at ice surface", "",
            "m s-1", "m year-1", 1);

  double fill_value = convert(m_sys, m_config->get_double("fill_value"), "m year-1", "m second-1");

  m_vars[0].set_double("valid_min", convert(m_sys, -1e6, "m year-1", "m second-1"));
  m_vars[0].set_double("valid_max", convert(m_sys, 1e6, "m year-1", "m second-1"));
  m_vars[0].set_double("_FillValue", fill_value);

  m_vars[1].set_double("valid_min", convert(m_sys, -1e6, "m year-1", "m second-1"));
  m_vars[1].set_double("valid_max", convert(m_sys, 1e6, "m year-1", "m second-1"));
  m_vars[1].set_double("_FillValue", fill_value);
}

IceModelVec::Ptr PSB_velsurf::compute_impl() {
  double fill_value = convert(m_sys, m_config->get_double("fill_value"), "m year-1", "m second-1");

  IceModelVec2V::Ptr result(new IceModelVec2V);
  result->create(m_grid, "surf", WITHOUT_GHOSTS);
  result->metadata(0) = m_vars[0];
  result->metadata(1) = m_vars[1];

  IceModelVec2S tmp;
  tmp.create(m_grid, "tmp", WITHOUT_GHOSTS);

  const IceModelVec3
    &u3 = model->velocity_u(),
    &v3 = model->velocity_v();

  const IceModelVec2S *thickness = m_grid->variables().get_2d_scalar("land_ice_thickness");

  u3.getSurfaceValues(tmp, *thickness);
  result->set_component(0, tmp);

  v3.getSurfaceValues(tmp, *thickness);
  result->set_component(1, tmp);

  const IceModelVec2CellType &mask = *m_grid->variables().get_2d_cell_type("mask");

  IceModelVec::AccessList list;
  list.add(mask);
  list.add(*result);

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (mask.ice_free(i, j)) {
      (*result)(i, j).u = fill_value;
      (*result)(i, j).v = fill_value;
    }
  }

  return result;
}

PSB_wvel::PSB_wvel(StressBalance *m)
  : Diag<StressBalance>(m) {

  // set metadata:
  m_vars.push_back(SpatialVariableMetadata(m_sys, "wvel",
                                           m_grid->z()));

  set_attrs("vertical velocity of ice, relative to geoid", "",
            "m s-1", "m year-1", 0);
  m_vars[0].set_double("valid_min", convert(m_sys, -1e6, "m year-1", "m second-1"));
  m_vars[0].set_double("valid_max", convert(m_sys, 1e6, "m year-1", "m second-1"));
}

IceModelVec::Ptr PSB_wvel::compute(bool zero_above_ice) {
  IceModelVec3::Ptr result3(new IceModelVec3);
  result3->create(m_grid, "wvel", WITHOUT_GHOSTS);
  result3->metadata() = m_vars[0];

  const IceModelVec2S *bed, *uplift;
  bed    = m_grid->variables().get_2d_scalar("bedrock_altitude");
  uplift = m_grid->variables().get_2d_scalar("tendency_of_bedrock_altitude");

  const IceModelVec2S        &thickness = *m_grid->variables().get_2d_scalar("land_ice_thickness");
  const IceModelVec2CellType &mask      = *m_grid->variables().get_2d_cell_type("mask");

  const IceModelVec3
    &u3 = model->velocity_u(),
    &v3 = model->velocity_v(),
    &w3 = model->velocity_w();

  IceModelVec::AccessList list;
  list.add(thickness);
  list.add(mask);
  list.add(*bed);
  list.add(u3);
  list.add(v3);
  list.add(w3);
  list.add(*uplift);
  list.add(*result3);

  const double ice_density = m_config->get_double("ice_density"),
    sea_water_density = m_config->get_double("sea_water_density"),
    R = ice_density / sea_water_density;

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      const double
        *u = u3.get_column(i, j),
        *v = v3.get_column(i, j),
        *w = w3.get_column(i, j);
      double *result = result3->get_column(i, j);

      int ks = m_grid->kBelowHeight(thickness(i,j));

      // in the ice:
      if (mask.grounded(i,j)) {
        const double
          bed_dx = bed->diff_x_p(i,j),
          bed_dy = bed->diff_y_p(i,j),
          uplift_ij = (*uplift)(i,j);
        for (int k = 0; k <= ks ; k++) {
          result[k] = w[k] + uplift_ij + u[k] * bed_dx + v[k] * bed_dy;
        }

      } else {                  // floating
        const double
          z_sl = R * thickness(i,j),
          w_sl = w3.getValZ(i, j, z_sl);

        for (int k = 0; k <= ks ; k++) {
          result[k] = w[k] - w_sl;
        }

      }

      // above the ice:
      if (zero_above_ice) {
        // set to zeros
        for (unsigned int k = ks+1; k < m_grid->Mz() ; k++) {
          result[k] = 0.0;
        }
      } else {
        // extrapolate using the topmost value
        for (unsigned int k = ks+1; k < m_grid->Mz() ; k++) {
          result[k] = result[ks];
        }
      }

    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  return result3;
}

IceModelVec::Ptr PSB_wvel::compute_impl() {
  return this->compute(true);   // fill wvel above the ice with zeros
}

PSB_wvelsurf::PSB_wvelsurf(StressBalance *m)
  : Diag<StressBalance>(m) {

  // set metadata:
  m_vars.push_back(SpatialVariableMetadata(m_sys, "wvelsurf"));

  set_attrs("vertical velocity of ice at ice surface, relative to the geoid", "",
            "m s-1", "m year-1", 0);
  m_vars[0].set_double("valid_min", convert(m_sys, -1e6, "m year-1", "m second-1"));
  m_vars[0].set_double("valid_max", convert(m_sys, 1e6, "m year-1", "m second-1"));

  double fill_value = convert(m_sys, m_config->get_double("fill_value"), "m year-1", "m second-1");
  m_vars[0].set_double("_FillValue", fill_value);
}

IceModelVec::Ptr PSB_wvelsurf::compute_impl() {
  double fill_value = convert(m_sys, m_config->get_double("fill_value"), "m year-1", "m second-1");

  IceModelVec2S::Ptr result(new IceModelVec2S);
  result->create(m_grid, "wvelsurf", WITHOUT_GHOSTS);
  result->metadata() = m_vars[0];

  // here "false" means "don't fill w3 above the ice surface with zeros"
  IceModelVec3::Ptr w3 = IceModelVec3::To3DScalar(PSB_wvel(model).compute(false));

  const IceModelVec2S *thickness = m_grid->variables().get_2d_scalar("land_ice_thickness");

  w3->getSurfaceValues(*result, *thickness);

  const IceModelVec2CellType &mask = *m_grid->variables().get_2d_cell_type("mask");

  IceModelVec::AccessList list;
  list.add(mask);
  list.add(*result);

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (mask.ice_free(i, j)) {
      (*result)(i, j) = fill_value;
    }
  }

  return result;
}

PSB_wvelbase::PSB_wvelbase(StressBalance *m)
  : Diag<StressBalance>(m) {

  // set metadata:
  m_vars.push_back(SpatialVariableMetadata(m_sys, "wvelbase"));

  set_attrs("vertical velocity of ice at the base of ice, relative to the geoid", "",
            "m s-1", "m year-1", 0);
  m_vars[0].set_double("valid_min", convert(m_sys, -1e6, "m year-1", "m second-1"));
  m_vars[0].set_double("valid_max", convert(m_sys, 1e6, "m year-1", "m second-1"));

  double fill_value = convert(m_sys, m_config->get_double("fill_value"), "m year-1", "m second-1");
  m_vars[0].set_double("_FillValue", fill_value);
}

IceModelVec::Ptr PSB_wvelbase::compute_impl() {
  double fill_value = convert(m_sys, m_config->get_double("fill_value"), "m year-1", "m second-1");

  IceModelVec2S::Ptr result(new IceModelVec2S);
  result->create(m_grid, "wvelbase", WITHOUT_GHOSTS);
  result->metadata() = m_vars[0];

  // here "false" means "don't fill w3 above the ice surface with zeros"
  IceModelVec3::Ptr w3 = IceModelVec3::To3DScalar(PSB_wvel(model).compute(false));

  w3->getHorSlice(*result, 0.0);

  const IceModelVec2CellType &mask = *m_grid->variables().get_2d_cell_type("mask");

  IceModelVec::AccessList list;
  list.add(mask);
  list.add(*result);

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (mask.ice_free(i, j)) {
      (*result)(i, j) = fill_value;
    }
  }

  return result;
}

PSB_velbase::PSB_velbase(StressBalance *m)
  : Diag<StressBalance>(m) {

  // set metadata:
  m_dof = 2;

  m_vars.push_back(SpatialVariableMetadata(m_sys, "uvelbase"));
  m_vars.push_back(SpatialVariableMetadata(m_sys, "vvelbase"));

  set_attrs("x-component of the horizontal velocity of ice at the base of ice", "",
            "m s-1", "m year-1", 0);
  set_attrs("y-component of the horizontal velocity of ice at the base of ice", "",
            "m s-1", "m year-1", 1);

  double fill_value = convert(m_sys, m_config->get_double("fill_value"), "m year-1", "m second-1");

  m_vars[0].set_double("valid_min", convert(m_sys, -1e6, "m year-1", "m second-1"));
  m_vars[0].set_double("valid_max", convert(m_sys, 1e6, "m year-1", "m second-1"));
  m_vars[0].set_double("_FillValue", fill_value);

  m_vars[1].set_double("valid_min", convert(m_sys, -1e6, "m year-1", "m second-1"));
  m_vars[1].set_double("valid_max", convert(m_sys, 1e6, "m year-1", "m second-1"));
  m_vars[1].set_double("_FillValue", fill_value);
}

IceModelVec::Ptr PSB_velbase::compute_impl() {
  double fill_value = convert(m_sys, m_config->get_double("fill_value"), "m year-1", "m second-1");

  IceModelVec2V::Ptr result(new IceModelVec2V);
  result->create(m_grid, "base", WITHOUT_GHOSTS);
  result->metadata(0) = m_vars[0];
  result->metadata(1) = m_vars[1];

  IceModelVec2S tmp;            // will be de-allocated automatically
  tmp.create(m_grid, "tmp", WITHOUT_GHOSTS);

  const IceModelVec3
    &u3 = model->velocity_u(),
    &v3 = model->velocity_v();

  u3.getHorSlice(tmp, 0.0);
  result->set_component(0, tmp);

  v3.getHorSlice(tmp, 0.0);
  result->set_component(1, tmp);

  const IceModelVec2CellType &mask = *m_grid->variables().get_2d_cell_type("mask");

  IceModelVec::AccessList list;
  list.add(mask);
  list.add(*result);

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (mask.ice_free(i, j)) {
      (*result)(i, j).u = fill_value;
      (*result)(i, j).v = fill_value;
    }
  }

  return result;
}


PSB_bfrict::PSB_bfrict(StressBalance *m)
  : Diag<StressBalance>(m) {

  // set metadata:
  m_vars.push_back(SpatialVariableMetadata(m_sys, "bfrict"));

  set_attrs("basal frictional heating", "",
            "W m-2", "W m-2", 0);
}

IceModelVec::Ptr PSB_bfrict::compute_impl() {

  IceModelVec2S::Ptr result(new IceModelVec2S);
  result->create(m_grid, "bfrict", WITHOUT_GHOSTS);
  result->metadata() = m_vars[0];

  result->copy_from(model->basal_frictional_heating());

  return result;
}


PSB_uvel::PSB_uvel(StressBalance *m)
  : Diag<StressBalance>(m) {

  // set metadata:
  m_vars.push_back(SpatialVariableMetadata(m_sys, "uvel", m_grid->z()));

  set_attrs("horizontal velocity of ice in the X direction", "land_ice_x_velocity",
            "m s-1", "m year-1", 0);
}

IceModelVec::Ptr PSB_uvel::compute_impl() {

  IceModelVec3::Ptr result(new IceModelVec3);
  result->create(m_grid, "uvel", WITHOUT_GHOSTS);
  result->metadata() = m_vars[0];

  const IceModelVec2S *thickness = m_grid->variables().get_2d_scalar("land_ice_thickness");

  const IceModelVec3
    &u3 = model->velocity_u();

  IceModelVec::AccessList list;
  list.add(u3);
  list.add(*result);
  list.add(*thickness);

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      int ks = m_grid->kBelowHeight((*thickness)(i,j));

      const double *u_ij = u3.get_column(i,j);
      double *u_out_ij = result->get_column(i,j);

      // in the ice:
      for (int k = 0; k <= ks ; k++) {
        u_out_ij[k] = u_ij[k];
      }
      // above the ice:
      for (unsigned int k = ks+1; k < m_grid->Mz() ; k++) {
        u_out_ij[k] = 0.0;
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  return result;
}

PSB_vvel::PSB_vvel(StressBalance *m)
  : Diag<StressBalance>(m) {

  // set metadata:
  m_vars.push_back(SpatialVariableMetadata(m_sys, "vvel", m_grid->z()));

  set_attrs("horizontal velocity of ice in the Y direction", "land_ice_y_velocity",
            "m s-1", "m year-1", 0);
}

IceModelVec::Ptr PSB_vvel::compute_impl() {

  IceModelVec3::Ptr result(new IceModelVec3);
  result->create(m_grid, "vvel", WITHOUT_GHOSTS);
  result->metadata() = m_vars[0];

  const IceModelVec2S *thickness = m_grid->variables().get_2d_scalar("land_ice_thickness");

  const IceModelVec3
    &v3 = model->velocity_v();

  IceModelVec::AccessList list;
  list.add(v3);
  list.add(*result);
  list.add(*thickness);

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      int ks = m_grid->kBelowHeight((*thickness)(i,j));

      const double *v_ij = v3.get_column(i,j);
      double *v_out_ij = result->get_column(i,j);

      // in the ice:
      for (int k = 0; k <= ks ; k++) {
        v_out_ij[k] = v_ij[k];
      }
      // above the ice:
      for (unsigned int k = ks+1; k < m_grid->Mz() ; k++) {
        v_out_ij[k] = 0.0;
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  return result;
}

PSB_wvel_rel::PSB_wvel_rel(StressBalance *m)
  : Diag<StressBalance>(m) {

  // set metadata:
  m_vars.push_back(SpatialVariableMetadata(m_sys, "wvel_rel", m_grid->z()));

  set_attrs("vertical velocity of ice, relative to base of ice directly below", "",
            "m s-1", "m year-1", 0);
}

IceModelVec::Ptr PSB_wvel_rel::compute_impl() {

  IceModelVec3::Ptr result(new IceModelVec3);
  result->create(m_grid, "wvel_rel", WITHOUT_GHOSTS);
  result->metadata() = m_vars[0];

  const IceModelVec2S *thickness = m_grid->variables().get_2d_scalar("land_ice_thickness");

  const IceModelVec3
    &w3 = model->velocity_w();

  IceModelVec::AccessList list;
  list.add(w3);
  list.add(*result);
  list.add(*thickness);

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      int ks = m_grid->kBelowHeight((*thickness)(i,j));

      const double *w_ij = w3.get_column(i,j);
      double *w_out_ij = result->get_column(i,j);

      // in the ice:
      for (int k = 0; k <= ks ; k++) {
        w_out_ij[k] = w_ij[k];
      }
      // above the ice:
      for (unsigned int k = ks+1; k < m_grid->Mz() ; k++) {
        w_out_ij[k] = 0.0;
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();


  return result;
}


PSB_strainheat::PSB_strainheat(StressBalance *m)
  : Diag<StressBalance>(m) {

  // set metadata:
  m_vars.push_back(SpatialVariableMetadata(m_sys, "strainheat", m_grid->z()));

  set_attrs("rate of strain heating in ice (dissipation heating)", "",
            "W m-3", "mW m-3", 0);
}

IceModelVec::Ptr PSB_strainheat::compute_impl() {

  IceModelVec3::Ptr result(new IceModelVec3);
  result->create(m_grid, "strainheat", WITHOUT_GHOSTS);
  result->metadata() = m_vars[0];
  result->write_in_glaciological_units = true;

  result->copy_from(model->volumetric_strain_heating());

  return result;
}

PSB_strain_rates::PSB_strain_rates(StressBalance *m)
  : Diag<StressBalance>(m) {
  m_dof = 2;

  // set metadata:
  m_vars.push_back(SpatialVariableMetadata(m_sys, "eigen1"));
  m_vars.push_back(SpatialVariableMetadata(m_sys, "eigen2"));

  set_attrs("first eigenvalue of the horizontal, vertically-integrated strain rate tensor",
            "", "s-1", "s-1", 0);
  set_attrs("second eigenvalue of the horizontal, vertically-integrated strain rate tensor",
            "", "s-1", "s-1", 1);
}

IceModelVec::Ptr PSB_strain_rates::compute_impl() {
  IceModelVec2V::Ptr velbar = IceModelVec2V::ToVector(PSB_velbar(model).compute());

  IceModelVec2::Ptr result(new IceModelVec2);
  result->create(m_grid, "strain_rates", WITHOUT_GHOSTS, 1, 2);
  result->metadata(0) = m_vars[0];
  result->metadata(1) = m_vars[1];

  const IceModelVec2CellType &mask = *m_grid->variables().get_2d_cell_type("mask");

  IceModelVec2V velbar_with_ghosts;
  velbar_with_ghosts.create(m_grid, "velbar", WITH_GHOSTS);

  // copy_from communicates ghosts
  velbar_with_ghosts.copy_from(*velbar);

  model->compute_2D_principal_strain_rates(velbar_with_ghosts, mask, *result);

  return result;
}

PSB_deviatoric_stresses::PSB_deviatoric_stresses(StressBalance *m)
  : Diag<StressBalance>(m) {
  m_dof = 3;

  // set metadata:
  m_vars.push_back(SpatialVariableMetadata(m_sys, "sigma_xx"));
  m_vars.push_back(SpatialVariableMetadata(m_sys, "sigma_yy"));
  m_vars.push_back(SpatialVariableMetadata(m_sys, "sigma_xy"));

  set_attrs("deviatoric stress in X direction", "", "Pa", "Pa", 0);
  set_attrs("deviatoric stress in Y direction", "", "Pa", "Pa", 1);
  set_attrs("deviatoric shear stress", "", "Pa", "Pa", 2);

}

IceModelVec::Ptr PSB_deviatoric_stresses::compute_impl() {

  IceModelVec2::Ptr velbar = IceModelVec2V::ToVector(PSB_velbar(model).compute());

  IceModelVec2::Ptr result(new IceModelVec2);
  result->create(m_grid, "deviatoric_stresses", WITHOUT_GHOSTS, 1, 3);
  result->metadata(0) = m_vars[0];
  result->metadata(1) = m_vars[1];
  result->metadata(2) = m_vars[2];

  const IceModelVec2CellType &mask = *m_grid->variables().get_2d_cell_type("mask");

  IceModelVec2V velbar_with_ghosts;
  velbar_with_ghosts.create(m_grid, "velbar", WITH_GHOSTS);

  // copy_from communicates ghosts
  velbar_with_ghosts.copy_from(*velbar);

  model->compute_2D_stresses(velbar_with_ghosts, mask, *result);

  return result;
}

PSB_pressure::PSB_pressure(StressBalance *m)
  : Diag<StressBalance>(m) {

  // set metadata:
  m_vars.push_back(SpatialVariableMetadata(m_sys, "pressure",
                                           m_grid->z()));

  set_attrs("pressure in ice (hydrostatic)", "", "Pa", "Pa", 0);
}

IceModelVec::Ptr PSB_pressure::compute_impl() {

  IceModelVec3::Ptr result(new IceModelVec3);
  result->create(m_grid, "pressure", WITHOUT_GHOSTS);
  result->metadata(0) = m_vars[0];

  const IceModelVec2S *thickness = m_grid->variables().get_2d_scalar("land_ice_thickness");

  IceModelVec::AccessList list;
  list.add(*result);
  list.add(*thickness);

  const double rg = m_config->get_double("ice_density") * m_config->get_double("standard_gravity");

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      unsigned int ks = m_grid->kBelowHeight((*thickness)(i,j));
      double *P_out_ij = result->get_column(i,j);
      const double H = (*thickness)(i,j);
      // within the ice:
      for (unsigned int k = 0; k <= ks; ++k) {
        P_out_ij[k] = rg * (H - m_grid->z(k));  // FIXME: add atmospheric pressure?
      }
      // above the ice:
      for (unsigned int k = ks + 1; k < m_grid->Mz(); ++k) {
        P_out_ij[k] = 0.0;  // FIXME: use atmospheric pressure?
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  return result;
}


PSB_tauxz::PSB_tauxz(StressBalance *m)
  : Diag<StressBalance>(m) {

  // set metadata:
  m_vars.push_back(SpatialVariableMetadata(m_sys, "tauxz", m_grid->z()));

  set_attrs("shear stress xz component (in shallow ice approximation SIA)", "",
            "Pa", "Pa", 0);
}


/*!
 * The SIA-applicable shear stress component tauxz computed here is not used
 * by the model.  This implementation intentionally does not use the
 * eta-transformation or special cases at ice margins.
 * CODE DUPLICATION WITH PSB_tauyz
 */
IceModelVec::Ptr PSB_tauxz::compute_impl() {

  IceModelVec3::Ptr result(new IceModelVec3);
  result->create(m_grid, "tauxz", WITHOUT_GHOSTS);
  result->metadata() = m_vars[0];

  const IceModelVec2S *thickness, *surface;

  thickness = m_grid->variables().get_2d_scalar("land_ice_thickness");
  surface   = m_grid->variables().get_2d_scalar("surface_altitude");

  IceModelVec::AccessList list;
  list.add(*result);
  list.add(*surface);
  list.add(*thickness);

  const double rg = m_config->get_double("ice_density") * m_config->get_double("standard_gravity");

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      unsigned int ks = m_grid->kBelowHeight((*thickness)(i,j));
      double *tauxz_out_ij = result->get_column(i, j);
      const double
        H    = (*thickness)(i,j),
        dhdx = surface->diff_x_p(i,j);

      // within the ice:
      for (unsigned int k = 0; k <= ks; ++k) {
        tauxz_out_ij[k] = - rg * (H - m_grid->z(k)) * dhdx;
      }
      // above the ice:
      for (unsigned int k = ks + 1; k < m_grid->Mz(); ++k) {
        tauxz_out_ij[k] = 0.0;
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  return result;
}


PSB_tauyz::PSB_tauyz(StressBalance *m)
  : Diag<StressBalance>(m) {

  // set metadata:
  m_vars.push_back(SpatialVariableMetadata(m_sys, "tauyz",
                                           m_grid->z()));

  set_attrs("shear stress yz component (in shallow ice approximation SIA)", "",
            "Pa", "Pa", 0);
}


/*!
 * The SIA-applicable shear stress component tauyz computed here is not used
 * by the model.  This implementation intentionally does not use the
 * eta-transformation or special cases at ice margins.
 * CODE DUPLICATION WITH PSB_tauxz
 */
IceModelVec::Ptr PSB_tauyz::compute_impl() {

  IceModelVec3::Ptr result(new IceModelVec3);
  result->create(m_grid, "tauyz", WITHOUT_GHOSTS);
  result->metadata(0) = m_vars[0];

  const IceModelVec2S *thickness = m_grid->variables().get_2d_scalar("land_ice_thickness");
  const IceModelVec2S *surface   = m_grid->variables().get_2d_scalar("surface_altitude");

  IceModelVec::AccessList list;
  list.add(*result);
  list.add(*surface);
  list.add(*thickness);

  const double rg = m_config->get_double("ice_density") * m_config->get_double("standard_gravity");

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      unsigned int ks = m_grid->kBelowHeight((*thickness)(i,j));
      double *tauyz_out_ij = result->get_column(i, j);
      const double
        H    = (*thickness)(i,j),
        dhdy = surface->diff_y_p(i,j);

      // within the ice:
      for (unsigned int k = 0; k <= ks; ++k) {
        tauyz_out_ij[k] = - rg * (H - m_grid->z(k)) * dhdy;
      }
      // above the ice:
      for (unsigned int k = ks + 1; k < m_grid->Mz(); ++k) {
        tauyz_out_ij[k] = 0.0;
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  return result;
}


} // end of namespace stressbalance
} // end of namespace pism
