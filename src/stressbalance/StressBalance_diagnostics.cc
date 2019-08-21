// Copyright (C) 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018 Constantine Khroulev
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

#include "StressBalance_diagnostics.hh"
#include "SSB_Modifier.hh"
#include "ShallowStressBalance.hh"
#include "pism/util/Mask.hh"
#include "pism/util/ConfigInterface.hh"
#include "pism/util/Vars.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/util/IceModelVec2CellType.hh"
#include "pism/rheology/FlowLaw.hh"

namespace pism {
namespace stressbalance {

using units::convert;

DiagnosticList StressBalance::diagnostics_impl() const {
  DiagnosticList result = {
    {"bfrict",              Diagnostic::Ptr(new PSB_bfrict(this))},
    {"velbar_mag",          Diagnostic::Ptr(new PSB_velbar_mag(this))},
    {"flux",                Diagnostic::Ptr(new PSB_flux(this))},
    {"flux_mag",            Diagnostic::Ptr(new PSB_flux_mag(this))},
    {"velbase_mag",         Diagnostic::Ptr(new PSB_velbase_mag(this))},
    {"velsurf_mag",         Diagnostic::Ptr(new PSB_velsurf_mag(this))},
    {"uvel",                Diagnostic::Ptr(new PSB_uvel(this))},
    {"vvel",                Diagnostic::Ptr(new PSB_vvel(this))},
    {"strainheat",          Diagnostic::Ptr(new PSB_strainheat(this))},
    {"velbar",              Diagnostic::Ptr(new PSB_velbar(this))},
    {"velbase",             Diagnostic::Ptr(new PSB_velbase(this))},
    {"velsurf",             Diagnostic::Ptr(new PSB_velsurf(this))},
    {"wvel",                Diagnostic::Ptr(new PSB_wvel(this))},
    {"wvelbase",            Diagnostic::Ptr(new PSB_wvelbase(this))},
    {"wvelsurf",            Diagnostic::Ptr(new PSB_wvelsurf(this))},
    {"wvel_rel",            Diagnostic::Ptr(new PSB_wvel_rel(this))},
    {"strain_rates",        Diagnostic::Ptr(new PSB_strain_rates(this))},
    {"vonmises_stress",     Diagnostic::Ptr(new PSB_vonmises_stress(this))},
    {"deviatoric_stresses", Diagnostic::Ptr(new PSB_deviatoric_stresses(this))},
    {"pressure",            Diagnostic::Ptr(new PSB_pressure(this))},
    {"tauxz",               Diagnostic::Ptr(new PSB_tauxz(this))},
    {"tauyz",               Diagnostic::Ptr(new PSB_tauyz(this))}
  };

  // add diagnostics from the shallow stress balance and the "modifier"
  result = pism::combine(result, m_shallow_stress_balance->diagnostics());
  result = pism::combine(result, m_modifier->diagnostics());

  return result;
}

TSDiagnosticList StressBalance::ts_diagnostics_impl() const {
  return pism::combine(m_shallow_stress_balance->ts_diagnostics(),
                       m_modifier->ts_diagnostics());
}

PSB_velbar::PSB_velbar(const StressBalance *m)
  : Diag<StressBalance>(m) {

  // set metadata:
  m_vars = {SpatialVariableMetadata(m_sys, "ubar"),
            SpatialVariableMetadata(m_sys, "vbar")};

  set_attrs("vertical mean of horizontal ice velocity in the X direction",
            "land_ice_vertical_mean_x_velocity",
            "m s-1", "m year-1", 0);
  set_attrs("vertical mean of horizontal ice velocity in the Y direction",
            "land_ice_vertical_mean_y_velocity",
            "m s-1", "m year-1", 1);
}

IceModelVec::Ptr PSB_velbar::compute_impl() const {
  // get the thickness
  const IceModelVec2S* thickness = m_grid->variables().get_2d_scalar("land_ice_thickness");

  // Compute the vertically-integrated horizontal ice flux:
  IceModelVec2V::Ptr result = IceModelVec2V::ToVector(PSB_flux(model).compute());

  // Override metadata set by the flux computation
  result->metadata(0) = m_vars[0];
  result->metadata(1) = m_vars[1];

  IceModelVec::AccessList list{thickness, result.get()};

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

PSB_velbar_mag::PSB_velbar_mag(const StressBalance *m)
  : Diag<StressBalance>(m) {

  // set metadata:
  m_vars = {SpatialVariableMetadata(m_sys, "velbar_mag")};

  set_attrs("magnitude of vertically-integrated horizontal velocity of ice", "",
            "m second-1", "m year-1", 0);
  m_vars[0].set_double("_FillValue", convert(m_sys, m_fill_value,
                                             "m year-1", "m second-1"));
  m_vars[0].set_double("valid_min", 0.0);
}

IceModelVec::Ptr PSB_velbar_mag::compute_impl() const {

  IceModelVec2S::Ptr result(new IceModelVec2S(m_grid, "velbar_mag", WITHOUT_GHOSTS));
  result->metadata(0) = m_vars[0];

  // compute vertically-averaged horizontal velocity:
  IceModelVec2V::Ptr velbar = IceModelVec2V::ToVector(PSB_velbar(model).compute());

  // compute its magnitude:
  result->set_to_magnitude(*velbar);

  const IceModelVec2S *thickness = m_grid->variables().get_2d_scalar("land_ice_thickness");

  // mask out ice-free areas:
  double fill_value = convert(m_sys, m_fill_value, "m year-1", "m second-1");
  result->mask_by(*thickness, fill_value);

  return result;
}


PSB_flux::PSB_flux(const StressBalance *m)
  : Diag<StressBalance>(m) {

  // set metadata:
  m_vars = {SpatialVariableMetadata(m_sys, "uflux"),
            SpatialVariableMetadata(m_sys, "vflux")};

  set_attrs("Vertically integrated horizontal flux of ice in the X direction",
            "",                 // no CF standard name
            "m2 s-1", "m2 year-1", 0);
  set_attrs("Vertically integrated horizontal flux of ice in the Y direction",
            "",                 // no CF standard name
            "m2 s-1", "m2 year-1", 1);
}

IceModelVec::Ptr PSB_flux::compute_impl() const {
  double H_threshold = m_config->get_double("geometry.ice_free_thickness_standard");

  IceModelVec2V::Ptr result(new IceModelVec2V);
  result->create(m_grid, "flux", WITHOUT_GHOSTS);
  result->metadata(0) = m_vars[0];
  result->metadata(1) = m_vars[1];

  // get the thickness
  const IceModelVec2S *thickness = m_grid->variables().get_2d_scalar("land_ice_thickness");

  const IceModelVec3
    &u3 = model->velocity_u(),
    &v3 = model->velocity_v();

  IceModelVec::AccessList list{&u3, &v3, thickness, result.get()};

  auto &z = m_grid->z();

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      double H = (*thickness)(i,j);

      // an "ice-free" cell:
      if (H < H_threshold) {
        (*result)(i, j) = 0.0;
        continue;
      }

      // an icy cell:
      {
        auto u = u3.get_column(i, j);
        auto v = v3.get_column(i, j);

        Vector2 Q(0.0, 0.0);

        // ks is "k just below the surface"
        int ks = m_grid->kBelowHeight(H);

        if (ks > 0) {
          Vector2 v0(u[0], v[0]);

          for (int k = 1; k <= ks; ++k) {
            Vector2 v1(u[k], v[k]);

            // trapezoid rule
            Q += (z[k] - z[k - 1]) * 0.5 * (v0 + v1);

            v0 = v1;
          }
        }

        // rectangle method to integrate over the last level
        Q += (H - z[ks]) * Vector2(u[ks], v[ks]);

        (*result)(i, j) = Q;
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  return result;
}


PSB_flux_mag::PSB_flux_mag(const StressBalance *m)
  : Diag<StressBalance>(m) {

  // set metadata:
  m_vars = {SpatialVariableMetadata(m_sys, "flux_mag")};

  set_attrs("magnitude of vertically-integrated horizontal flux of ice", "",
            "m2 s-1", "m2 year-1", 0);

  double fill_value = convert(m_sys, m_fill_value,
                              "m2 year-1", "m2 second-1");
  m_vars[0].set_double("_FillValue", fill_value);
  m_vars[0].set_double("valid_min", 0.0);
}

IceModelVec::Ptr PSB_flux_mag::compute_impl() const {
  const IceModelVec2S *thickness = m_grid->variables().get_2d_scalar("land_ice_thickness");

  // Compute the vertically-averaged horizontal ice velocity:
  IceModelVec2S::Ptr result = IceModelVec2S::To2DScalar(PSB_velbar_mag(model).compute());

  IceModelVec::AccessList list{thickness, result.get()};

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    (*result)(i,j) *= (*thickness)(i,j);
  }

  double fill_value = convert(m_sys, m_fill_value, "m year-1", "m second-1");
  result->mask_by(*thickness, fill_value);

  result->metadata() = m_vars[0];

  return result;
}

PSB_velbase_mag::PSB_velbase_mag(const StressBalance *m)
  : Diag<StressBalance>(m) {

  // set metadata:
  m_vars = {SpatialVariableMetadata(m_sys, "velbase_mag")};

  set_attrs("magnitude of horizontal velocity of ice at base of ice", "",
            "m s-1", "m year-1", 0);

  double fill_value = convert(m_sys, m_fill_value, "m year-1", "m second-1");

  m_vars[0].set_double("_FillValue", fill_value);
  m_vars[0].set_double("valid_min", 0.0);
}

IceModelVec::Ptr PSB_velbase_mag::compute_impl() const {
  // FIXME: compute this using PSB_velbase.

  IceModelVec2S tmp(m_grid, "tmp", WITHOUT_GHOSTS);

  IceModelVec2S::Ptr result(new IceModelVec2S(m_grid, "velbase_mag", WITHOUT_GHOSTS));
  result->metadata(0) = m_vars[0];

  const IceModelVec3
    &u3 = model->velocity_u(),
    &v3 = model->velocity_v();

  const IceModelVec2S *thickness = m_grid->variables().get_2d_scalar("land_ice_thickness");

  u3.getHorSlice(*result, 0.0); // result = u_{z=0}
  v3.getHorSlice(tmp, 0.0);    // tmp = v_{z=0}

  result->set_to_magnitude(*result, tmp);

  // mask out ice-free areas
  double fill_value = convert(m_sys, m_fill_value, "m year-1", "m second-1");
  result->mask_by(*thickness, fill_value);

  return result;
}

PSB_velsurf_mag::PSB_velsurf_mag(const StressBalance *m)
  : Diag<StressBalance>(m) {
  // set metadata:
  m_vars = {SpatialVariableMetadata(m_sys, "velsurf_mag")};

  set_attrs("magnitude of horizontal velocity of ice at ice surface", "",
            "m s-1", "m year-1", 0);

  double fill_value = convert(m_sys, m_fill_value, "m year-1", "m second-1");

  m_vars[0].set_double("_FillValue", fill_value);
  m_vars[0].set_double("valid_min",  0.0);
}

IceModelVec::Ptr PSB_velsurf_mag::compute_impl() const {

  // FIXME: Compute this using PSB_velsurf.

  IceModelVec2S tmp;
  tmp.create(m_grid, "tmp", WITHOUT_GHOSTS);

  IceModelVec2S::Ptr result(new IceModelVec2S(m_grid, "velsurf_mag", WITHOUT_GHOSTS));
  result->metadata(0) = m_vars[0];

  const IceModelVec3
    &u3 = model->velocity_u(),
    &v3 = model->velocity_v();

  const IceModelVec2S *thickness = m_grid->variables().get_2d_scalar("land_ice_thickness");

  u3.getSurfaceValues(*result, *thickness);
  v3.getSurfaceValues(tmp, *thickness);

  result->set_to_magnitude(*result, tmp);

  // mask out ice-free areas
  double fill_value = convert(m_sys, m_fill_value, "m year-1", "m second-1");
  result->mask_by(*thickness, fill_value);

  return result;
}


PSB_velsurf::PSB_velsurf(const StressBalance *m)
  : Diag<StressBalance>(m) {

  // set metadata:
  m_vars = {SpatialVariableMetadata(m_sys, "uvelsurf"),
            SpatialVariableMetadata(m_sys, "vvelsurf")};

  set_attrs("x-component of the horizontal velocity of ice at ice surface",
            "land_ice_surface_x_velocity", // InitMIP "standard" name
            "m s-1", "m year-1", 0);
  set_attrs("y-component of the horizontal velocity of ice at ice surface",
            "land_ice_surface_y_velocity", // InitMIP "standard" name
            "m s-1", "m year-1", 1);

  double fill_value = convert(m_sys, m_fill_value, "m year-1", "m second-1");

  m_vars[0].set_double("valid_min", convert(m_sys, -1e6, "m year-1", "m second-1"));
  m_vars[0].set_double("valid_max", convert(m_sys, 1e6, "m year-1", "m second-1"));
  m_vars[0].set_double("_FillValue", fill_value);

  m_vars[1].set_double("valid_min", convert(m_sys, -1e6, "m year-1", "m second-1"));
  m_vars[1].set_double("valid_max", convert(m_sys, 1e6, "m year-1", "m second-1"));
  m_vars[1].set_double("_FillValue", fill_value);
}

IceModelVec::Ptr PSB_velsurf::compute_impl() const {
  double fill_value = convert(m_sys, m_fill_value, "m year-1", "m second-1");

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

  IceModelVec::AccessList list{&mask, result.get()};

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (mask.ice_free(i, j)) {
      (*result)(i, j).u = fill_value;
      (*result)(i, j).v = fill_value;
    }
  }

  return result;
}

PSB_wvel::PSB_wvel(const StressBalance *m)
  : Diag<StressBalance>(m) {

  // set metadata:
  m_vars = {SpatialVariableMetadata(m_sys, "wvel", m_grid->z())};

  set_attrs("vertical velocity of ice, relative to geoid", "",
            "m s-1", "m year-1", 0);
  m_vars[0].set_double("valid_min", convert(m_sys, -1e6, "m year-1", "m second-1"));
  m_vars[0].set_double("valid_max", convert(m_sys, 1e6, "m year-1", "m second-1"));
}

IceModelVec::Ptr PSB_wvel::compute(bool zero_above_ice) const {
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

  IceModelVec::AccessList list{&thickness, &mask, bed, &u3, &v3, &w3, uplift, result3.get()};

  const double ice_density = m_config->get_double("constants.ice.density"),
    sea_water_density = m_config->get_double("constants.sea_water.density"),
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

IceModelVec::Ptr PSB_wvel::compute_impl() const {
  return this->compute(true);   // fill wvel above the ice with zeros
}

PSB_wvelsurf::PSB_wvelsurf(const StressBalance *m)
  : Diag<StressBalance>(m) {

  // set metadata:
  m_vars = {SpatialVariableMetadata(m_sys, "wvelsurf")};

  set_attrs("vertical velocity of ice at ice surface, relative to the geoid",
            "land_ice_surface_upward_velocity", // InitMIP "standard" name
            "m s-1", "m year-1", 0);
  m_vars[0].set_double("valid_min", convert(m_sys, -1e6, "m year-1", "m second-1"));
  m_vars[0].set_double("valid_max", convert(m_sys, 1e6, "m year-1", "m second-1"));

  double fill_value = convert(m_sys, m_fill_value, "m year-1", "m second-1");
  m_vars[0].set_double("_FillValue", fill_value);
}

IceModelVec::Ptr PSB_wvelsurf::compute_impl() const {
  double fill_value = convert(m_sys, m_fill_value, "m year-1", "m second-1");

  IceModelVec2S::Ptr result(new IceModelVec2S(m_grid, "wvelsurf", WITHOUT_GHOSTS));
  result->metadata() = m_vars[0];

  // here "false" means "don't fill w3 above the ice surface with zeros"
  IceModelVec3::Ptr w3 = IceModelVec3::To3DScalar(PSB_wvel(model).compute(false));

  const IceModelVec2S *thickness = m_grid->variables().get_2d_scalar("land_ice_thickness");

  w3->getSurfaceValues(*result, *thickness);

  const IceModelVec2CellType &mask = *m_grid->variables().get_2d_cell_type("mask");

  IceModelVec::AccessList list{&mask, result.get()};

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (mask.ice_free(i, j)) {
      (*result)(i, j) = fill_value;
    }
  }

  return result;
}

PSB_wvelbase::PSB_wvelbase(const StressBalance *m)
  : Diag<StressBalance>(m) {

  // set metadata:
  m_vars = {SpatialVariableMetadata(m_sys, "wvelbase")};

  set_attrs("vertical velocity of ice at the base of ice, relative to the geoid",
            "land_ice_basal_upward_velocity", // InitMIP "standard" name
            "m s-1", "m year-1", 0);
  m_vars[0].set_double("valid_min", convert(m_sys, -1e6, "m year-1", "m second-1"));
  m_vars[0].set_double("valid_max", convert(m_sys, 1e6, "m year-1", "m second-1"));

  double fill_value = convert(m_sys, m_fill_value, "m year-1", "m second-1");
  m_vars[0].set_double("_FillValue", fill_value);
}

IceModelVec::Ptr PSB_wvelbase::compute_impl() const {
  double fill_value = convert(m_sys, m_fill_value, "m year-1", "m second-1");

  IceModelVec2S::Ptr result(new IceModelVec2S(m_grid, "wvelbase", WITHOUT_GHOSTS));
  result->metadata() = m_vars[0];

  // here "false" means "don't fill w3 above the ice surface with zeros"
  IceModelVec3::Ptr w3 = IceModelVec3::To3DScalar(PSB_wvel(model).compute(false));

  w3->getHorSlice(*result, 0.0);

  const IceModelVec2CellType &mask = *m_grid->variables().get_2d_cell_type("mask");

  IceModelVec::AccessList list{&mask, result.get()};

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (mask.ice_free(i, j)) {
      (*result)(i, j) = fill_value;
    }
  }

  return result;
}

PSB_velbase::PSB_velbase(const StressBalance *m)
  : Diag<StressBalance>(m) {

  // set metadata:
  m_vars = {SpatialVariableMetadata(m_sys, "uvelbase"),
            SpatialVariableMetadata(m_sys, "vvelbase")};

  set_attrs("x-component of the horizontal velocity of ice at the base of ice",
            "land_ice_basal_x_velocity", // InitMIP "standard" name
            "m s-1", "m year-1", 0);
  set_attrs("y-component of the horizontal velocity of ice at the base of ice",
            "land_ice_basal_y_velocity", // InitMIP "standard" name
            "m s-1", "m year-1", 1);

  double fill_value = convert(m_sys, m_fill_value, "m year-1", "m second-1");

  m_vars[0].set_double("valid_min", convert(m_sys, -1e6, "m year-1", "m second-1"));
  m_vars[0].set_double("valid_max", convert(m_sys, 1e6, "m year-1", "m second-1"));
  m_vars[0].set_double("_FillValue", fill_value);

  m_vars[1].set_double("valid_min", convert(m_sys, -1e6, "m year-1", "m second-1"));
  m_vars[1].set_double("valid_max", convert(m_sys, 1e6, "m year-1", "m second-1"));
  m_vars[1].set_double("_FillValue", fill_value);
}

IceModelVec::Ptr PSB_velbase::compute_impl() const {
  double fill_value = convert(m_sys, m_fill_value, "m year-1", "m second-1");

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

  IceModelVec::AccessList list{&mask, result.get()};

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (mask.ice_free(i, j)) {
      (*result)(i, j).u = fill_value;
      (*result)(i, j).v = fill_value;
    }
  }

  return result;
}


PSB_bfrict::PSB_bfrict(const StressBalance *m)
  : Diag<StressBalance>(m) {

  // set metadata:
  m_vars = {SpatialVariableMetadata(m_sys, "bfrict")};

  set_attrs("basal frictional heating", "",
            "W m-2", "W m-2", 0);
}

IceModelVec::Ptr PSB_bfrict::compute_impl() const {

  IceModelVec2S::Ptr result(new IceModelVec2S(m_grid, "bfrict", WITHOUT_GHOSTS));
  result->metadata() = m_vars[0];

  result->copy_from(model->basal_frictional_heating());

  return result;
}


PSB_uvel::PSB_uvel(const StressBalance *m)
  : Diag<StressBalance>(m) {

  // set metadata:
  m_vars = {SpatialVariableMetadata(m_sys, "uvel", m_grid->z())};

  set_attrs("horizontal velocity of ice in the X direction", "land_ice_x_velocity",
            "m s-1", "m year-1", 0);
}

/*!
 * Copy F to result and set it to zero above the surface of the ice.
 */
static void zero_above_ice(const IceModelVec3 &F, const IceModelVec2S &H,
                           IceModelVec3 &result) {

  IceModelVec::AccessList list{&F, &H, &result};

  IceGrid::ConstPtr grid = result.grid();

  auto Mz = grid->Mz();

  ParallelSection loop(grid->com);
  try {
    for (Points p(*grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      int ks = grid->kBelowHeight(H(i,j));

      const double *F_ij = F.get_column(i,j);
      double *F_out_ij = result.get_column(i,j);

      // in the ice:
      for (int k = 0; k <= ks ; k++) {
        F_out_ij[k] = F_ij[k];
      }
      // above the ice:
      for (unsigned int k = ks+1; k < Mz ; k++) {
        F_out_ij[k] = 0.0;
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();
}

IceModelVec::Ptr PSB_uvel::compute_impl() const {

  IceModelVec3::Ptr result(new IceModelVec3);
  result->create(m_grid, "uvel", WITHOUT_GHOSTS);
  result->metadata() = m_vars[0];

  zero_above_ice(model->velocity_u(),
                 *m_grid->variables().get_2d_scalar("land_ice_thickness"),
                 *result);

  return result;
}

PSB_vvel::PSB_vvel(const StressBalance *m)
  : Diag<StressBalance>(m) {

  // set metadata:
  m_vars = {SpatialVariableMetadata(m_sys, "vvel", m_grid->z())};

  set_attrs("horizontal velocity of ice in the Y direction", "land_ice_y_velocity",
            "m s-1", "m year-1", 0);
}

IceModelVec::Ptr PSB_vvel::compute_impl() const {

  IceModelVec3::Ptr result(new IceModelVec3);
  result->create(m_grid, "vvel", WITHOUT_GHOSTS);
  result->metadata() = m_vars[0];

  zero_above_ice(model->velocity_v(),
                 *m_grid->variables().get_2d_scalar("land_ice_thickness"),
                 *result);

  return result;
}

PSB_wvel_rel::PSB_wvel_rel(const StressBalance *m)
  : Diag<StressBalance>(m) {

  // set metadata:
  m_vars = {SpatialVariableMetadata(m_sys, "wvel_rel", m_grid->z())};

  set_attrs("vertical velocity of ice, relative to base of ice directly below", "",
            "m s-1", "m year-1", 0);
}

IceModelVec::Ptr PSB_wvel_rel::compute_impl() const {

  IceModelVec3::Ptr result(new IceModelVec3);
  result->create(m_grid, "wvel_rel", WITHOUT_GHOSTS);
  result->metadata() = m_vars[0];

  zero_above_ice(model->velocity_w(),
                 *m_grid->variables().get_2d_scalar("land_ice_thickness"),
                 *result);

  return result;
}


PSB_strainheat::PSB_strainheat(const StressBalance *m)
  : Diag<StressBalance>(m) {

  // set metadata:
  m_vars = {SpatialVariableMetadata(m_sys, "strainheat", m_grid->z())};

  set_attrs("rate of strain heating in ice (dissipation heating)", "",
            "W m-3", "mW m-3", 0);
}

IceModelVec::Ptr PSB_strainheat::compute_impl() const {
  IceModelVec3::Ptr result(new IceModelVec3(m_grid, "strainheat", WITHOUT_GHOSTS));
  result->metadata() = m_vars[0];

  result->copy_from(model->volumetric_strain_heating());

  return result;
}

PSB_strain_rates::PSB_strain_rates(const StressBalance *m)
  : Diag<StressBalance>(m) {
  // set metadata:
  m_vars = {SpatialVariableMetadata(m_sys, "eigen1"),
            SpatialVariableMetadata(m_sys, "eigen2")};

  set_attrs("first eigenvalue of the horizontal, vertically-integrated strain rate tensor",
            "", "s-1", "s-1", 0);
  set_attrs("second eigenvalue of the horizontal, vertically-integrated strain rate tensor",
            "", "s-1", "s-1", 1);
}

IceModelVec::Ptr PSB_strain_rates::compute_impl() const {
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

  compute_2D_principal_strain_rates(velbar_with_ghosts, mask, *result);

  return result;
}

PSB_deviatoric_stresses::PSB_deviatoric_stresses(const StressBalance *m)
  : Diag<StressBalance>(m) {
  // set metadata:
  m_vars = {SpatialVariableMetadata(m_sys, "sigma_xx"),
            SpatialVariableMetadata(m_sys, "sigma_yy"),
            SpatialVariableMetadata(m_sys, "sigma_xy")};

  set_attrs("deviatoric stress in X direction", "", "Pa", "Pa", 0);
  set_attrs("deviatoric stress in Y direction", "", "Pa", "Pa", 1);
  set_attrs("deviatoric shear stress", "", "Pa", "Pa", 2);

}

IceModelVec::Ptr PSB_deviatoric_stresses::compute_impl() const {

  IceModelVec2::Ptr result(new IceModelVec2);
  result->create(m_grid, "deviatoric_stresses", WITHOUT_GHOSTS, 1, 3);
  result->metadata(0) = m_vars[0];
  result->metadata(1) = m_vars[1];
  result->metadata(2) = m_vars[2];

  const IceModelVec2CellType &cell_type = *m_grid->variables().get_2d_cell_type("mask");
  const IceModelVec3         *enthalpy  = m_grid->variables().get_3d_scalar("enthalpy");
  const IceModelVec2S        *thickness = m_grid->variables().get_2d_scalar("land_ice_thickness");

  IceModelVec2S hardness(m_grid, "hardness", WITHOUT_GHOSTS);
  IceModelVec2V velocity(m_grid, "velocity", WITH_GHOSTS);

  averaged_hardness_vec(*model->shallow()->flow_law(), *thickness, *enthalpy,
                        hardness);

  // copy_from updates ghosts
  velocity.copy_from(*IceModelVec2V::ToVector(PSB_velbar(model).compute()));

  model->compute_2D_stresses(velocity, hardness, cell_type, *result);

  return result;
}

PSB_pressure::PSB_pressure(const StressBalance *m)
  : Diag<StressBalance>(m) {

  // set metadata:
  m_vars = {SpatialVariableMetadata(m_sys, "pressure", m_grid->z())};

  set_attrs("pressure in ice (hydrostatic)", "", "Pa", "Pa", 0);
}

IceModelVec::Ptr PSB_pressure::compute_impl() const {

  IceModelVec3::Ptr result(new IceModelVec3);
  result->create(m_grid, "pressure", WITHOUT_GHOSTS);
  result->metadata(0) = m_vars[0];

  const IceModelVec2S *thickness = m_grid->variables().get_2d_scalar("land_ice_thickness");

  IceModelVec::AccessList list{thickness, result.get()};

  const double rg = m_config->get_double("constants.ice.density") * m_config->get_double("constants.standard_gravity");

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


PSB_tauxz::PSB_tauxz(const StressBalance *m)
  : Diag<StressBalance>(m) {

  // set metadata:
  m_vars = {SpatialVariableMetadata(m_sys, "tauxz", m_grid->z())};

  set_attrs("shear stress xz component (in shallow ice approximation SIA)", "",
            "Pa", "Pa", 0);
}


/*!
 * The SIA-applicable shear stress component tauxz computed here is not used
 * by the model.  This implementation intentionally does not use the
 * eta-transformation or special cases at ice margins.
 * CODE DUPLICATION WITH PSB_tauyz
 */
IceModelVec::Ptr PSB_tauxz::compute_impl() const {

  IceModelVec3::Ptr result(new IceModelVec3);
  result->create(m_grid, "tauxz", WITHOUT_GHOSTS);
  result->metadata() = m_vars[0];

  const IceModelVec2S *thickness, *surface;

  thickness = m_grid->variables().get_2d_scalar("land_ice_thickness");
  surface   = m_grid->variables().get_2d_scalar("surface_altitude");

  IceModelVec::AccessList list{surface, thickness, result.get()};

  const double rg = m_config->get_double("constants.ice.density") * m_config->get_double("constants.standard_gravity");

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


PSB_tauyz::PSB_tauyz(const StressBalance *m)
  : Diag<StressBalance>(m) {

  // set metadata:
  m_vars = {SpatialVariableMetadata(m_sys, "tauyz", m_grid->z())};

  set_attrs("shear stress yz component (in shallow ice approximation SIA)", "",
            "Pa", "Pa", 0);
}


/*!
 * The SIA-applicable shear stress component tauyz computed here is not used
 * by the model.  This implementation intentionally does not use the
 * eta-transformation or special cases at ice margins.
 * CODE DUPLICATION WITH PSB_tauxz
 */
IceModelVec::Ptr PSB_tauyz::compute_impl() const {

  IceModelVec3::Ptr result(new IceModelVec3);
  result->create(m_grid, "tauyz", WITHOUT_GHOSTS);
  result->metadata(0) = m_vars[0];

  const IceModelVec2S *thickness = m_grid->variables().get_2d_scalar("land_ice_thickness");
  const IceModelVec2S *surface   = m_grid->variables().get_2d_scalar("surface_altitude");

  IceModelVec::AccessList list{surface, thickness, result.get()};

  const double rg = m_config->get_double("constants.ice.density") * m_config->get_double("constants.standard_gravity");

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

PSB_vonmises_stress::PSB_vonmises_stress(const StressBalance *m)
  : Diag<StressBalance>(m) {

  /* set metadata: */
  m_vars = {SpatialVariableMetadata(m_sys, "vonmises_stress")};

  set_attrs("tensile von Mises stress",
            "",                 // no standard name
            "Pascal", "Pascal", 0);
}

IceModelVec::Ptr PSB_vonmises_stress::compute_impl() const {

  using std::max;

  IceModelVec2S::Ptr result(new IceModelVec2S(m_grid, "vonmises_stress", WITHOUT_GHOSTS));
  result->metadata(0) = m_vars[0];

  IceModelVec2S &vonmises_stress = *result;

  IceModelVec2V::Ptr velbar = IceModelVec2V::ToVector(PSB_velbar(model).compute());
  IceModelVec2V &velocity = *velbar;

  IceModelVec2::Ptr eigen12 = IceModelVec2::To2D(PSB_strain_rates(model).compute());
  IceModelVec2 &strain_rates = *eigen12;

  const IceModelVec2S &ice_thickness = *m_grid->variables().get_2d_scalar("land_ice_thickness");
  const IceModelVec3 *enthalpy = m_grid->variables().get_3d_scalar("enthalpy");
  const IceModelVec2CellType &mask = *m_grid->variables().get_2d_cell_type("mask");

  const rheology::FlowLaw &flow_law = *model->shallow()->flow_law();

  const double *z = &m_grid->z()[0];
  const double ssa_n = flow_law.exponent();

  IceModelVec::AccessList list{&vonmises_stress, &velocity, &strain_rates, &ice_thickness,
      enthalpy, &mask};

  for (Points pt(*m_grid); pt; pt.next()) {
    const int i = pt.i(), j = pt.j();

    // Find partially filled or empty grid boxes on the icefree ocean, which
    // have floating ice neighbors after the mass continuity step
    if (mask.icy(i, j)) {

      const double       H = ice_thickness(i, j);
      const unsigned int k = m_grid->kBelowHeight(H);

      const double
        *enthalpy_column   = enthalpy->get_column(i, j),
        hardness           = averaged_hardness(flow_law, H, k, z, enthalpy_column),
        eigen1             = strain_rates(i, j, 0),
        eigen2             = strain_rates(i, j, 1);

      // [\ref Morlighem2016] equation 6
      const double effective_tensile_strain_rate = sqrt(0.5 * (PetscSqr(max(0.0, eigen1)) +
                                                               PetscSqr(max(0.0, eigen2))));
      // [\ref Morlighem2016] equation 7
      vonmises_stress(i, j) = sqrt(3.0) * hardness * pow(effective_tensile_strain_rate,
                                                         1.0 / ssa_n);

    } else { // end of "if (mask.icy(i, j))"
      vonmises_stress(i, j) = 0.0;
    }
  }   // end of loop over grid points

  return result;
}

} // end of namespace stressbalance
} // end of namespace pism
