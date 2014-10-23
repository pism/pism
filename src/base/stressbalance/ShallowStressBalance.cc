// Copyright (C) 2010, 2011, 2012, 2013, 2014 Constantine Khroulev and Ed Bueler
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

#include "ShallowStressBalance.hh"
#include "Mask.hh"
#include "PISMVars.hh"
#include "flowlaws.hh"
#include "basal_resistance.hh"
#include "pism_options.hh"
#include <cassert>

#include "error_handling.hh"

namespace pism {

using mask::ice_free;

ShallowStressBalance::ShallowStressBalance(IceGrid &g, EnthalpyConverter &e, const Config &conf)
  : Component(g, conf), basal_sliding_law(NULL), flow_law(NULL), EC(e) {

  m_vel_bc = NULL;
  bc_locations = NULL;
  variables = NULL;
  sea_level = 0;

  allocate();
}

ShallowStressBalance::~ShallowStressBalance() {
  delete basal_sliding_law;
}

//! \brief Allocate a shallow stress balance object.
PetscErrorCode ShallowStressBalance::allocate() {
  PetscErrorCode ierr;
  const unsigned int WIDE_STENCIL = config.get("grid_max_stencil_width");

  if (config.get_flag("do_pseudo_plastic_till") == true)
    basal_sliding_law = new IceBasalResistancePseudoPlasticLaw(config);
  else
    basal_sliding_law = new IceBasalResistancePlasticLaw(config);

  ierr = m_velocity.create(grid, "bar", WITH_GHOSTS, WIDE_STENCIL); CHKERRQ(ierr); // components ubar, vbar
  ierr = m_velocity.set_attrs("model_state",
                              "thickness-advective ice velocity (x-component)", 
                              "m s-1", "", 0); CHKERRQ(ierr);
  ierr = m_velocity.set_attrs("model_state",
                              "thickness-advective ice velocity (y-component)",
                              "m s-1", "", 1); CHKERRQ(ierr);
  ierr = m_velocity.set_glaciological_units("m year-1"); CHKERRQ(ierr);
  m_velocity.write_in_glaciological_units = true;

  ierr = basal_frictional_heating.create(grid, "bfrict", WITHOUT_GHOSTS); CHKERRQ(ierr);
  ierr = basal_frictional_heating.set_attrs("diagnostic",
                                            "basal frictional heating",
                                            "W m-2", ""); CHKERRQ(ierr);
  ierr = basal_frictional_heating.set_glaciological_units("mW m-2"); CHKERRQ(ierr);
  basal_frictional_heating.write_in_glaciological_units = true;

  return 0;
}

void ShallowStressBalance::get_diagnostics(std::map<std::string, Diagnostic*> &dict,
                                           std::map<std::string, TSDiagnostic*> &/*ts_dict*/) {
  dict["beta"]     = new SSB_beta(this, grid, *variables);
  dict["taub"]     = new SSB_taub(this, grid, *variables);
  dict["taub_mag"] = new SSB_taub_mag(this, grid, *variables);
  dict["taud"]     = new SSB_taud(this, grid, *variables);
  dict["taud_mag"] = new SSB_taud_mag(this, grid, *variables);
}


ZeroSliding::ZeroSliding(IceGrid &g, EnthalpyConverter &e, const Config &conf)
  : ShallowStressBalance(g, e, conf) {

  // Use the SIA flow law.
  IceFlowLawFactory ice_factory(grid.com, "sia_", config, &EC);
  ice_factory.setType(config.get_string("sia_flow_law"));

  ice_factory.setFromOptions();
  ice_factory.create(&flow_law);
}

ZeroSliding::~ZeroSliding() {
  delete flow_law;
}

void ZeroSliding::add_vars_to_output(const std::string &/*keyword*/, std::set<std::string> &/*result*/)
{
  // empty
}

PetscErrorCode ZeroSliding::define_variables(const std::set<std::string> &/*vars*/, const PIO &/*nc*/,
                                             IO_Type /*nctype*/) {
  return 0;
}

PetscErrorCode ZeroSliding::write_variables(const std::set<std::string> &/*vars*/, const PIO &/*nc*/) {
  return 0;
}


//! \brief Update the trivial shallow stress balance object.
PetscErrorCode ZeroSliding::update(bool fast, IceModelVec2S &melange_back_pressure) {
  PetscErrorCode ierr;
  if (fast)
    return 0;

  (void) melange_back_pressure;

  ierr = m_velocity.set(0.0); CHKERRQ(ierr);

  ierr = basal_frictional_heating.set(0.0); CHKERRQ(ierr);

  return 0;
}

//! \brief Compute the basal frictional heating.
/*!
  Ice shelves have zero basal friction heating.

  \param[in] velocity *basal* sliding velocity
  \param[in] tauc basal yield stress
  \param[in] mask (used to determine if floating or grounded)
  \param[out] result
 */
PetscErrorCode ShallowStressBalance::compute_basal_frictional_heating(IceModelVec2V &velocity,
                                                                      IceModelVec2S &tauc,
                                                                      IceModelVec2Int &mask,
                                                                      IceModelVec2S &result) {
  MaskQuery m(mask);

  IceModelVec::AccessList list;
  list.add(velocity);
  list.add(result);
  list.add(tauc);
  list.add(mask);
  
  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (m.ocean(i,j)) {
      result(i,j) = 0.0;
    } else {
      const double
        C = basal_sliding_law->drag(tauc(i,j), velocity(i,j).u, velocity(i,j).v),
        basal_stress_x = - C * velocity(i,j).u,
        basal_stress_y = - C * velocity(i,j).v;
      result(i,j) = - basal_stress_x * velocity(i,j).u - basal_stress_y * velocity(i,j).v;
    }
  }

  return 0;
}


//! \brief Compute eigenvalues of the horizontal, vertically-integrated strain rate tensor.
/*!
Calculates all components \f$D_{xx}, D_{yy}, D_{xy}=D_{yx}\f$ of the
vertically-averaged strain rate tensor \f$D\f$ [\ref SchoofStream].  Then computes
the eigenvalues `result(i,j,0)` = (maximum eigenvalue), `result(i,j,1)` = (minimum
eigenvalue).  Uses the provided thickness to make decisions (PIK) about computing
strain rates near calving front.

Note that `result(i,j,0)` >= `result(i,j,1)`, but there is no necessary relation between 
the magnitudes, and either principal strain rate could be negative or positive.

Result can be used in a calving law, for example in eigencalving (PIK).

Note: strain rates will be derived from SSA velocities, using ghosts when
necessary. Both implementations (SSAFD and SSAFEM) call
update_ghosts() to ensure that ghost values are up to date.
 */
PetscErrorCode ShallowStressBalance::compute_2D_principal_strain_rates(IceModelVec2V &velocity, IceModelVec2Int &mask,
                                                                       IceModelVec2 &result) {
  double    dx = grid.dx, dy = grid.dy;

  if (result.get_ndof() != 2) {
    throw RuntimeError("result.dof() == 2 is required");
  }

  IceModelVec::AccessList list;
  list.add(velocity);
  list.add(result);
  list.add(mask);

  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (ice_free(mask.as_int(i,j))) {
      result(i,j,0) = 0.0;
      result(i,j,1) = 0.0;
      continue;
    }

    planeStar<int> m = mask.int_star(i,j);
    planeStar<Vector2> U = velocity.star(i,j);

    // strain in units s-1
    double u_x = 0, u_y = 0, v_x = 0, v_y = 0,
      east = 1, west = 1, south = 1, north = 1;

    // Computes u_x using second-order centered finite differences written as
    // weighted sums of first-order one-sided finite differences.
    //
    // Given the cell layout
    // *----n----*
    // |         |
    // |         |
    // w         e
    // |         |
    // |         |
    // *----s----*
    // east == 0 if the east neighbor of the current cell is ice-free. In
    // this case we use the left- (west-) sided difference.
    //
    // If both neighbors in the east-west (x) direction are ice-free the
    // x-derivative is set to zero (see u_x, v_x initialization above).
    //
    // Similarly in other directions.
    if (ice_free(m.e))
      east = 0;
    if (ice_free(m.w))
      west = 0;
    if (ice_free(m.n))
      north = 0;
    if (ice_free(m.s))
      south = 0;

    if (west + east > 0) {
      u_x = 1.0 / (dx * (west + east)) * (west * (U.ij.u - U[West].u) + east * (U[East].u - U.ij.u));
      v_x = 1.0 / (dx * (west + east)) * (west * (U.ij.v - U[West].v) + east * (U[East].v - U.ij.v));
    }

    if (south + north > 0) {
      u_y = 1.0 / (dy * (south + north)) * (south * (U.ij.u - U[South].u) + north * (U[North].u - U.ij.u));
      v_y = 1.0 / (dy * (south + north)) * (south * (U.ij.v - U[South].v) + north * (U[North].v - U.ij.v));
    }

    const double A = 0.5 * (u_x + v_y),  // A = (1/2) trace(D)
      B   = 0.5 * (u_x - v_y),
      Dxy = 0.5 * (v_x + u_y),  // B^2 = A^2 - u_x v_y
      q   = sqrt(PetscSqr(B) + PetscSqr(Dxy));
    result(i,j,0) = A + q;
    result(i,j,1) = A - q; // q >= 0 so e1 >= e2

  }


  return 0;
}

//! \brief Compute 2D deviatoric stresses.
/*! Note: IceModelVec2 result has to have dof == 3. */
PetscErrorCode ShallowStressBalance::compute_2D_stresses(IceModelVec2V &velocity, IceModelVec2Int &mask,
                                                         IceModelVec2 &result) {
  double    dx = grid.dx, dy = grid.dy;

  if (result.get_ndof() != 3) {
    throw RuntimeError("result.get_dof() == 3 is required");
  }

  // NB: uses constant ice hardness; choice is to use SSA's exponent; see issue #285
  double hardness = pow(config.get("ice_softness"),-1.0/config.get("ssa_Glen_exponent"));

  IceModelVec::AccessList list;
  list.add(velocity);
  list.add(result);
  list.add(mask);

  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (ice_free(mask.as_int(i,j))) {
      result(i,j,0) = 0.0;
      result(i,j,1) = 0.0;
      result(i,j,2) = 0.0;
      continue;
    }

    planeStar<int> m = mask.int_star(i,j);
    planeStar<Vector2> U = velocity.star(i,j);

    // strain in units s-1
    double u_x = 0, u_y = 0, v_x = 0, v_y = 0,
      east = 1, west = 1, south = 1, north = 1;

    // Computes u_x using second-order centered finite differences written as
    // weighted sums of first-order one-sided finite differences.
    //
    // Given the cell layout
    // *----n----*
    // |         |
    // |         |
    // w         e
    // |         |
    // |         |
    // *----s----*
    // east == 0 if the east neighbor of the current cell is ice-free. In
    // this case we use the left- (west-) sided difference.
    //
    // If both neighbors in the east-west (x) direction are ice-free the
    // x-derivative is set to zero (see u_x, v_x initialization above).
    //
    // Similarly in y-direction.
    if (ice_free(m.e))
      east = 0;
    if (ice_free(m.w))
      west = 0;
    if (ice_free(m.n))
      north = 0;
    if (ice_free(m.s))
      south = 0;

    if (west + east > 0) {
      u_x = 1.0 / (dx * (west + east)) * (west * (U.ij.u - U[West].u) + east * (U[East].u - U.ij.u));
      v_x = 1.0 / (dx * (west + east)) * (west * (U.ij.v - U[West].v) + east * (U[East].v - U.ij.v));
    }

    if (south + north > 0) {
      u_y = 1.0 / (dy * (south + north)) * (south * (U.ij.u - U[South].u) + north * (U[North].u - U.ij.u));
      v_y = 1.0 / (dy * (south + north)) * (south * (U.ij.v - U[South].v) + north * (U[North].v - U.ij.v));
    }

    double nu;
    flow_law->effective_viscosity(hardness,
                                  secondInvariant_2D(u_x, u_y, v_x, v_y),
                                  &nu, NULL);

    //get deviatoric stresses
    result(i,j,0) = nu*u_x;
    result(i,j,1) = nu*v_y;
    result(i,j,2) = 0.5*nu*(u_y+v_x);

  }


  return 0;
}

SSB_taud::SSB_taud(ShallowStressBalance *m, IceGrid &g, Vars &my_vars)
  : Diag<ShallowStressBalance>(m, g, my_vars) {

  dof = 2;

  // set metadata:
  vars.push_back(NCSpatialVariable(grid.get_unit_system(), "taud_x", grid));
  vars.push_back(NCSpatialVariable(grid.get_unit_system(), "taud_y", grid));

  set_attrs("X-component of the driving shear stress at the base of ice", "",
            "Pa", "Pa", 0);
  set_attrs("Y-component of the driving shear stress at the base of ice", "",
            "Pa", "Pa", 1);

  for (int k = 0; k < dof; ++k)
    vars[k].set_string("comment",
                       "this field is purely diagnostic (not used by the model)");
}

/*!
 * The driving stress computed here is not used by the model, so this
 * implementation intentionally does not use the eta-transformation or special
 * cases at ice margins.
 */
PetscErrorCode SSB_taud::compute(IceModelVec* &output) {
  PetscErrorCode ierr;
  IceModelVec2S *thickness, *surface;

  IceModelVec2V *result = new IceModelVec2V;
  ierr = result->create(grid, "result", WITHOUT_GHOSTS); CHKERRQ(ierr);
  result->metadata() = vars[0];
  result->metadata(1) = vars[1];

  thickness = variables.get_2d_scalar("land_ice_thickness");
  surface = variables.get_2d_scalar("surface_altitude");

  double standard_gravity = grid.config.get("standard_gravity"),
    ice_density = grid.config.get("ice_density");

  IceModelVec::AccessList list;
  list.add(*result);
  list.add(*surface);
  list.add(*thickness);

  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    double pressure = ice_density * standard_gravity * (*thickness)(i,j);
    if (pressure <= 0.0) {
      (*result)(i,j).u = 0.0;
      (*result)(i,j).v = 0.0;
    } else {
      (*result)(i,j).u = - pressure * surface->diff_x_p(i,j);
      (*result)(i,j).v = - pressure * surface->diff_y_p(i,j);
    }
  }

  output = result;
  return 0;
}

SSB_taud_mag::SSB_taud_mag(ShallowStressBalance *m, IceGrid &g, Vars &my_vars)
  : Diag<ShallowStressBalance>(m, g, my_vars) {

  // set metadata:
  vars.push_back(NCSpatialVariable(grid.get_unit_system(), "taud_mag", grid));

  set_attrs("magnitude of the gravitational driving stress at the base of ice", "",
            "Pa", "Pa", 0);
  vars[0].set_string("comment",
                     "this field is purely diagnostic (not used by the model)");
}

PetscErrorCode SSB_taud_mag::compute(IceModelVec* &output) {
  PetscErrorCode ierr;

  // Allocate memory:
  IceModelVec2S *result = new IceModelVec2S;
  ierr = result->create(grid, "taud_mag", WITHOUT_GHOSTS); CHKERRQ(ierr);
  result->metadata() = vars[0];
  result->write_in_glaciological_units = true;

  IceModelVec* tmp;
  SSB_taud diag(model, grid, variables);

  ierr = diag.compute(tmp);

  IceModelVec2V *taud = dynamic_cast<IceModelVec2V*>(tmp);
  if (taud == NULL) {
    delete tmp;
    throw RuntimeError("expected an IceModelVec2V, but dynamic_cast failed");
  }

  ierr = taud->magnitude(*result); CHKERRQ(ierr);

  delete tmp;

  output = result;
  return 0;
}

SSB_taub::SSB_taub(ShallowStressBalance *m, IceGrid &g, Vars &my_vars)
  : Diag<ShallowStressBalance>(m, g, my_vars) {
  dof = 2;

  // set metadata:
  vars.push_back(NCSpatialVariable(grid.get_unit_system(), "taub_x", grid));
  vars.push_back(NCSpatialVariable(grid.get_unit_system(), "taub_y", grid));

  set_attrs("X-component of the shear stress at the base of ice", "",
            "Pa", "Pa", 0);
  set_attrs("Y-component of the shear stress at the base of ice", "",
            "Pa", "Pa", 1);

  for (int k = 0; k < dof; ++k)
    vars[k].set_string("comment",
                       "this field is purely diagnostic (not used by the model)");
}


PetscErrorCode SSB_taub::compute(IceModelVec* &output) {
  PetscErrorCode ierr;

  IceModelVec2V *result = new IceModelVec2V;
  ierr = result->create(grid, "result", WITHOUT_GHOSTS); CHKERRQ(ierr);
  result->metadata() = vars[0];
  result->metadata(1) = vars[1];

  IceModelVec2V *velocity;
  ierr = model->get_2D_advective_velocity(velocity); CHKERRQ(ierr);

  IceModelVec2V &vel = *velocity;
  IceModelVec2S *tauc = variables.get_2d_scalar("tauc");
  IceModelVec2Int *mask = variables.get_2d_mask("mask");

  const IceBasalResistancePlasticLaw *basal_sliding_law = model->get_sliding_law();

  MaskQuery m(*mask);

  IceModelVec::AccessList list;
  list.add(*result);
  list.add(*tauc);
  list.add(vel);
  list.add(*mask);
  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (m.grounded_ice(i,j)) {
      double beta = basal_sliding_law->drag((*tauc)(i,j), vel(i,j).u, vel(i,j).v);
      (*result)(i,j).u = - beta * vel(i,j).u;
      (*result)(i,j).v = - beta * vel(i,j).v;
    } else {
      (*result)(i,j).u = 0.0;
      (*result)(i,j).v = 0.0;
    }
  }

  output = result;

  return 0;
}

SSB_taub_mag::SSB_taub_mag(ShallowStressBalance *m, IceGrid &g, Vars &my_vars)
  : Diag<ShallowStressBalance>(m, g, my_vars) {

  // set metadata:
  vars.push_back(NCSpatialVariable(grid.get_unit_system(), "taub_mag", grid));

  set_attrs("magnitude of the basal shear stress at the base of ice", "",
            "Pa", "Pa", 0);
  vars[0].set_string("comment",
                     "this field is purely diagnostic (not used by the model)");
}

PetscErrorCode SSB_taub_mag::compute(IceModelVec* &output) {
  PetscErrorCode ierr;

  // Allocate memory:
  IceModelVec2S *result = new IceModelVec2S;
  ierr = result->create(grid, "taub_mag", WITHOUT_GHOSTS); CHKERRQ(ierr);
  result->metadata() = vars[0];
  result->write_in_glaciological_units = true;

  IceModelVec* tmp;
  SSB_taub diag(model, grid, variables);

  ierr = diag.compute(tmp);

  IceModelVec2V *taub = dynamic_cast<IceModelVec2V*>(tmp);
  if (taub == NULL) {
    delete tmp;
    throw RuntimeError("expected an IceModelVec2V, but dynamic_cast failed");
  }

  ierr = taub->magnitude(*result); CHKERRQ(ierr);

  delete tmp;

  output = result;
  return 0;
}

/**
 * Shallow stress balance class that reads `u` and `v` fields from a
 * file and holds them constant.
 *
 * The only use I can think of right now is testing.
 */
PrescribedSliding::PrescribedSliding(IceGrid &g, EnthalpyConverter &e, const Config &conf)
  : ZeroSliding(g, e, conf) {
  // empty
}

PrescribedSliding::~PrescribedSliding() {
  // empty
}

PetscErrorCode PrescribedSliding::update(bool fast, IceModelVec2S &melange_back_pressure) {
  PetscErrorCode ierr;
  if (fast == true)
    return 0;

  (void) melange_back_pressure;

  ierr = basal_frictional_heating.set(0.0); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode PrescribedSliding::init(Vars &vars) {
  PetscErrorCode ierr;
  ierr = ShallowStressBalance::init(vars); CHKERRQ(ierr);

  bool flag;
  std::string input_filename;
  ierr = OptionsString("-prescribed_sliding_file", "name of the file to read velocity fields from",
                           input_filename, flag); CHKERRQ(ierr);
  if (flag == false) {
    throw RuntimeError("option -prescribed_sliding_file is required.");
  }

  ierr = m_velocity.regrid(input_filename, CRITICAL); CHKERRQ(ierr);

  return 0;
}

SSB_beta::SSB_beta(ShallowStressBalance *m, IceGrid &g, Vars &my_vars)
  : Diag<ShallowStressBalance>(m, g, my_vars) {
  // set metadata:
  vars.push_back(NCSpatialVariable(grid.get_unit_system(), "beta", grid));

  set_attrs("basal drag coefficient", "", "Pa s / m", "Pa s / m", 0);
}

PetscErrorCode SSB_beta::compute(IceModelVec* &output) {
  PetscErrorCode ierr;

  // Allocate memory:
  IceModelVec2S *result = new IceModelVec2S;
  ierr = result->create(grid, "beta", WITHOUT_GHOSTS); CHKERRQ(ierr);
  result->metadata() = vars[0];
  result->write_in_glaciological_units = true;

  IceModelVec2S *tauc = variables.get_2d_scalar("tauc");

  const IceBasalResistancePlasticLaw *basal_sliding_law = model->get_sliding_law();

  IceModelVec2V *velocity;
  ierr = model->get_2D_advective_velocity(velocity); CHKERRQ(ierr);
  IceModelVec2V &vel = *velocity;

  IceModelVec::AccessList list;
  list.add(*result);
  list.add(*tauc);
  list.add(*velocity);
  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    (*result)(i,j) =  basal_sliding_law->drag((*tauc)(i,j), vel(i,j).u, vel(i,j).v);
  }

  output = result;
  return 0;
}

} // end of namespace pism
