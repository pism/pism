// Copyright (C) 2010, 2011, 2012, 2013 Constantine Khroulev and Ed Bueler
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

#include "ShallowStressBalance.hh"
#include "Mask.hh"
#include "PISMVars.hh"

//! \brief Allocate a shallow stress balance object.
PetscErrorCode ShallowStressBalance::allocate() {
  PetscErrorCode ierr;

  ierr = m_velocity.create(grid, "bar", true); CHKERRQ(ierr); // components ubar, vbar
  ierr = m_velocity.set_attrs("model_state",
                            "thickness-advective ice velocity (x-component)", 
                            "m s-1", "", 0); CHKERRQ(ierr);
  ierr = m_velocity.set_attrs("model_state",
                            "thickness-advective ice velocity (y-component)",
                            "m s-1", "", 1); CHKERRQ(ierr);
  ierr = m_velocity.set_glaciological_units("m year-1"); CHKERRQ(ierr);
  m_velocity.write_in_glaciological_units = true;

  ierr = basal_frictional_heating.create(grid, "bfrict", false); CHKERRQ(ierr);
  ierr = basal_frictional_heating.set_attrs("diagnostic",
                                            "basal frictional heating",
                                            "W m-2", ""); CHKERRQ(ierr);
  ierr = basal_frictional_heating.set_glaciological_units("mW m-2"); CHKERRQ(ierr);
  basal_frictional_heating.write_in_glaciological_units = true;

  ierr = D2.create(grid, "D2", true); CHKERRQ(ierr);
  ierr = D2.set_attrs("internal",
                      "(partial) square of the Frobenius norm of D_{ij}, the combined strain rates",
                      "", ""); CHKERRQ(ierr);

  return 0;
}

//! \brief Update the trivial shallow stress balance object.
PetscErrorCode SSB_Trivial::update(bool fast) {
  PetscErrorCode ierr;
  if (fast) return 0;

  ierr = m_velocity.set(0.0); CHKERRQ(ierr);

  max_u = max_v = 0.0;

  ierr = basal_frictional_heating.set(0.0); CHKERRQ(ierr);

  ierr = D2.set(0.0); CHKERRQ(ierr);
  
  return 0;
}

//! \brief Compute eigenvalues of the horizontal, vertically-integrated strain rate tensor.
/*!
Calculates all components \f$D_{xx}, D_{yy}, D_{xy}=D_{yx}\f$ of the
vertically-averaged strain rate tensor \f$D\f$ [\ref SchoofStream].  Then computes
the eigenvalues \c result(i,j,0) = (maximum eigenvalue), \c result(i,j,1) = (minimum
eigenvalue).  Uses the provided thickness to make decisions (PIK) about computing
strain rates near calving front.

Though there are two eigenvalues, such do not form a vector, so the output is not
an IceModelVec2V, though it could be a std::vector<IceModelVec2S> or such.

Note that \c result(i,j,0) >= \c result(i,j,1), but there is no necessary relation between 
the magnitudes, and either principal strain rate could be negative or positive.

Result can be used in a calving law, for example in eigencalving (PIK).

Note: strain rates will be derived from SSA velocities, using ghosts when
necessary. Both implementations (SSAFD and SSAFEM) call
update_ghosts() to ensure that ghost values are up to date.
 */
PetscErrorCode ShallowStressBalance::compute_2D_principal_strain_rates(IceModelVec2V &velocity, IceModelVec2Int &mask,
                                                                       IceModelVec2 &result) {
  PetscErrorCode ierr;
  PetscScalar    dx = grid.dx, dy = grid.dy;
  Mask M;

  if (result.get_dof() != 2)
    SETERRQ(grid.com, 1, "result.get_dof() == 2 is required");

  ierr = velocity.begin_access(); CHKERRQ(ierr);
  ierr = result.begin_access(); CHKERRQ(ierr);
  ierr = mask.begin_access(); CHKERRQ(ierr);

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {

      if (M.ice_free(mask.as_int(i,j))) {
        result(i,j,0) = 0.0;
        result(i,j,1) = 0.0;
        continue;
      }

      planeStar<int> m = mask.int_star(i,j);
      planeStar<PISMVector2> U = velocity.star(i,j);

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
      if (M.ice_free(m.e))
        east = 0;
      if (M.ice_free(m.w))
        west = 0;
      if (M.ice_free(m.n))
        north = 0;
      if (M.ice_free(m.s))
        south = 0;

      if (west + east > 0) {
        u_x = 1.0 / (dx * (west + east)) * (west * (U.ij.u - U[West].u) + east * (U[East].u - U.ij.u));
        v_x = 1.0 / (dx * (west + east)) * (west * (U.ij.v - U[West].v) + east * (U[East].v - U.ij.v));
      }

      if (south + north > 0) {
        u_y = 1.0 / (dy * (south + north)) * (south * (U.ij.u - U[South].u) + north * (U[North].u - U.ij.u));
        v_y = 1.0 / (dy * (south + north)) * (south * (U.ij.v - U[South].v) + north * (U[North].v - U.ij.v));
      }

      const PetscScalar A = 0.5 * (u_x + v_y),  // A = (1/2) trace(D)
        B   = 0.5 * (u_x - v_y),
        Dxy = 0.5 * (v_x + u_y),  // B^2 = A^2 - u_x v_y
        q   = sqrt(PetscSqr(B) + PetscSqr(Dxy));
      result(i,j,0) = A + q;
      result(i,j,1) = A - q; // q >= 0 so e1 >= e2

    } // j
  }   // i

  ierr = velocity.end_access(); CHKERRQ(ierr);
  ierr = result.end_access(); CHKERRQ(ierr);
  ierr = mask.end_access(); CHKERRQ(ierr);

  return 0;
}

//! \brief Compute 2D deviatoric stresses.
/*! Note: IceModelVec2 result has to have dof == 3. */
PetscErrorCode ShallowStressBalance::compute_2D_stresses(IceModelVec2V &velocity, IceModelVec2Int &mask,
                                                         IceModelVec2 &result) {
  PetscErrorCode ierr;
  PetscScalar    dx = grid.dx, dy = grid.dy;
  Mask M;

  if (result.get_dof() != 3)
    SETERRQ(grid.com, 1, "result.get_dof() == 3 is required");

  // NB: uses constant ice hardness.
  double hardness = pow(config.get("ice_softness"),-1.0/config.get("Glen_exponent"));

  ierr = velocity.begin_access(); CHKERRQ(ierr);
  ierr = result.begin_access(); CHKERRQ(ierr);
  ierr = mask.begin_access(); CHKERRQ(ierr);

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {

      if (M.ice_free(mask.as_int(i,j))) {
        result(i,j,0) = 0.0;
        result(i,j,1) = 0.0;
        result(i,j,2) = 0.0;
        continue;
      }

      planeStar<int> m = mask.int_star(i,j);
      planeStar<PISMVector2> U = velocity.star(i,j);

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
      if (M.ice_free(m.e))
        east = 0;
      if (M.ice_free(m.w))
        west = 0;
      if (M.ice_free(m.n))
        north = 0;
      if (M.ice_free(m.s))
        south = 0;

      if (west + east > 0) {
        u_x = 1.0 / (dx * (west + east)) * (west * (U.ij.u - U[West].u) + east * (U[East].u - U.ij.u));
        v_x = 1.0 / (dx * (west + east)) * (west * (U.ij.v - U[West].v) + east * (U[East].v - U.ij.v));
      }

      if (south + north > 0) {
        u_y = 1.0 / (dy * (south + north)) * (south * (U.ij.u - U[South].u) + north * (U[North].u - U.ij.u));
        v_y = 1.0 / (dy * (south + north)) * (south * (U.ij.v - U[South].v) + north * (U[North].v - U.ij.v));
      }

      double nu = flow_law->effective_viscosity(hardness, u_x, u_y, v_x, v_y);

      //get deviatoric stresses
      result(i,j,0) = nu*u_x;
      result(i,j,1) = nu*v_y;
      result(i,j,2) = 0.5*nu*(u_y+v_x);

    } // j
  }   // i

  ierr = velocity.end_access(); CHKERRQ(ierr);
  ierr = result.end_access(); CHKERRQ(ierr);
  ierr = mask.end_access(); CHKERRQ(ierr);

  return 0;
}

SSB_taud::SSB_taud(ShallowStressBalance *m, IceGrid &g, PISMVars &my_vars)
  : PISMDiag<ShallowStressBalance>(m, g, my_vars) {

  dof = 2;
  vars.resize(dof);
  // set metadata:
  vars[0].init_2d("taud_x", grid);
  vars[1].init_2d("taud_y", grid);

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
  ierr = result->create(grid, "result", false); CHKERRQ(ierr);
  ierr = result->set_metadata(vars[0], 0); CHKERRQ(ierr);
  ierr = result->set_metadata(vars[1], 1); CHKERRQ(ierr);

  thickness = dynamic_cast<IceModelVec2S*>(variables.get("land_ice_thickness"));
  if (thickness == NULL) SETERRQ(grid.com, 1, "land_ice_thickness is not available");

  surface = dynamic_cast<IceModelVec2S*>(variables.get("surface_altitude"));
  if (surface == NULL) SETERRQ(grid.com, 1, "surface_altitude is not available");

  PetscReal standard_gravity = grid.config.get("standard_gravity"),
    ice_density = grid.config.get("ice_density");

  ierr =    result->begin_access(); CHKERRQ(ierr);
  ierr =   surface->begin_access(); CHKERRQ(ierr);
  ierr = thickness->begin_access(); CHKERRQ(ierr);

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      PetscScalar pressure = ice_density * standard_gravity * (*thickness)(i,j);
      if (pressure <= 0.0) {
        (*result)(i,j).u = 0.0;
        (*result)(i,j).v = 0.0;
      } else {
        (*result)(i,j).u = - pressure * surface->diff_x_p(i,j);
        (*result)(i,j).v = - pressure * surface->diff_y_p(i,j);
      }
    }
  }

  ierr = thickness->end_access(); CHKERRQ(ierr);
  ierr =   surface->end_access(); CHKERRQ(ierr);
  ierr =    result->end_access(); CHKERRQ(ierr);

  output = result;
  return 0;
}

SSB_taud_mag::SSB_taud_mag(ShallowStressBalance *m, IceGrid &g, PISMVars &my_vars)
  : PISMDiag<ShallowStressBalance>(m, g, my_vars) {

  // set metadata:
  vars[0].init_2d("taud_mag", grid);

  set_attrs("magnitude of the driving shear stress at the base of ice", "",
            "Pa", "Pa", 0);
  vars[0].set_string("comment",
                     "this field is purely diagnostic (not used by the model)");
}

PetscErrorCode SSB_taud_mag::compute(IceModelVec* &output) {
  PetscErrorCode ierr;

  // Allocate memory:
  IceModelVec2S *result = new IceModelVec2S;
  ierr = result->create(grid, "taud_mag", false); CHKERRQ(ierr);
  ierr = result->set_metadata(vars[0], 0); CHKERRQ(ierr);
  result->write_in_glaciological_units = true;

  IceModelVec* tmp;
  SSB_taud diag(model, grid, variables);

  ierr = diag.compute(tmp);

  IceModelVec2V *taud = dynamic_cast<IceModelVec2V*>(tmp);
  if (taud == NULL)
    SETERRQ(grid.com, 1, "expected an IceModelVec2V, but dynamic_cast failed");

  ierr = taud->magnitude(*result); CHKERRQ(ierr);

  delete tmp;

  output = result;
  return 0;
}
