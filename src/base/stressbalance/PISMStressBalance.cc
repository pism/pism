// Copyright (C) 2010 Constantine Khroulev and Ed Bueler
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

#include "PISMStressBalance.hh"

PISMStressBalance::PISMStressBalance(IceGrid &g, IceFlowLaw &i,
                                     const NCConfigVariable &conf)
  : grid(g), ice(i), config(conf) {

  basal_melt_rate = NULL;
  stress_balance  = NULL;
  modifier        = NULL;
}

PISMStressBalance::~PISMStressBalance() {
  delete stress_balance;
  delete modifier;
}

//! \brief Initialize the PISMStressBalance object.
PetscErrorCode PISMStressBalance::init(PISMVars &vars) {
  PetscErrorCode ierr;

  // allocate the vertical velocity field:

  ierr = w.create(grid, "wvel_rel", false); CHKERRQ(ierr);
  ierr = w.set_attrs("diagnostic",
                     "vertical velocity of ice, relative to base of ice directly below",
                     "m s-1", ""); CHKERRQ(ierr);
  w.time_independent = false;
  ierr = w.set_glaciological_units("m year-1"); CHKERRQ(ierr);
  w.write_in_glaciological_units = true;

  if (!stress_balance) SETERRQ(1, "stress_balance == NULL");
  if (!modifier)       SETERRQ(1, "modifier == NULL");

  ierr = stress_balance->init(vars); CHKERRQ(ierr);   
  ierr = modifier->init(vars); CHKERRQ(ierr); 

  return 0;
}

PetscErrorCode PISMStressBalance::set_initial_guess(IceModelVec2V &guess) {
  PetscErrorCode ierr;
  ierr = stress_balance->set_initial_guess(guess); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode PISMStressBalance::read_initial_guess(string filename) {
  PetscErrorCode ierr;
  ierr = stress_balance->read_initial_guess(filename); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode PISMStressBalance::save_initial_guess(string filename) {
  PetscErrorCode ierr;
  ierr = stress_balance->save_initial_guess(filename); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode PISMStressBalance::set_boundary_conditions(IceModelVec2Mask &locations,
                                                          IceModelVec2V &velocities) {
  PetscErrorCode ierr;
  ierr = stress_balance->set_boundary_conditions(locations, velocities); CHKERRQ(ierr);
  return 0;
}

//! \brief Set the basal melt rate.
PetscErrorCode PISMStressBalance::set_basal_melt_rate(IceModelVec2S &bmr_input) {
  basal_melt_rate = &bmr_input;
  return 0;
}

//! \brief Performs the shallow stress balance computation.
PetscErrorCode PISMStressBalance::update(bool fast) {
  PetscErrorCode ierr;
  IceModelVec2V *velocity_2d;
  IceModelVec2S *D2;
  IceModelVec3  *u, *v;

  ierr = stress_balance->update(fast); CHKERRQ(ierr);

  ierr = stress_balance->get_advective_2d_velocity(velocity_2d); CHKERRQ(ierr); 
  ierr = stress_balance->get_D2(D2); CHKERRQ(ierr);

  ierr = modifier->update(velocity_2d, D2, fast); CHKERRQ(ierr);

  if (!fast) {
    ierr = modifier->get_horizontal_3d_velocity(u, v); CHKERRQ(ierr);

    ierr = compute_vertical_velocity(u, v, basal_melt_rate, w); CHKERRQ(ierr); 
  }

  return 0;
}

PetscErrorCode PISMStressBalance::get_advective_2d_velocity(IceModelVec2V* &result) {
  PetscErrorCode ierr;
  ierr = stress_balance->get_advective_2d_velocity(result); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode PISMStressBalance::get_diffusive_flux(IceModelVec2Stag* &result) {
  PetscErrorCode ierr;
  ierr = modifier->get_diffusive_flux(result); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode PISMStressBalance::get_max_diffusivity(PetscReal &D) {
  PetscErrorCode ierr;
  ierr = modifier->get_max_diffusivity(D); CHKERRQ(ierr);  
  return 0;
}

PetscErrorCode PISMStressBalance::get_max_2d_velocity(PetscReal &u_max, PetscReal &v_max) {
  PetscErrorCode ierr;
  ierr = stress_balance->get_max_2d_velocity(u_max, v_max); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode PISMStressBalance::get_3d_velocity(IceModelVec3* &u, IceModelVec3* &v, IceModelVec3* &w_out) {
  PetscErrorCode ierr;
  ierr = modifier->get_horizontal_3d_velocity(u, v); CHKERRQ(ierr);
  w_out = &w;
  return 0;
}

PetscErrorCode PISMStressBalance::get_max_3d_velocity(PetscReal &u, PetscReal &v, PetscReal &w_out) {
  PetscErrorCode ierr;
  ierr = modifier->get_max_horizontal_velocity(u, v); CHKERRQ(ierr);
  w_out = w_max;
  return 0;
}

PetscErrorCode PISMStressBalance::get_basal_frictional_heating(IceModelVec2S* &result) {
  PetscErrorCode ierr;
  ierr = stress_balance->get_basal_frictional_heating(result); CHKERRQ(ierr); 
  return 0;
}

//! \brief Extend the grid vertically.
PetscErrorCode PISMStressBalance::extend_the_grid(PetscInt old_Mz) {
  PetscErrorCode ierr;

  ierr = w.extend_vertically(old_Mz, 0.0); CHKERRQ(ierr);

  ierr = stress_balance->extend_the_grid(old_Mz); CHKERRQ(ierr);

  ierr = modifier->extend_the_grid(old_Mz); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PISMStressBalance::compute_vertical_velocity(IceModelVec3 *u, IceModelVec3 *v,
                                                            IceModelVec2S *bmr,
                                                            IceModelVec3 &result) {
  PetscErrorCode ierr;
  const PetscScalar dx = grid.dx, dy = grid.dy;

  ierr = u->begin_access(); CHKERRQ(ierr);
  ierr = v->begin_access(); CHKERRQ(ierr);
  ierr = result.begin_access(); CHKERRQ(ierr);

  if (bmr) {
    ierr = bmr->begin_access(); CHKERRQ(ierr);
  }

  PetscScalar *w, *u_im1, *u_ip1, *v_jm1, *v_jp1;

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      ierr = result.getInternalColumn(i,j,&w); CHKERRQ(ierr);

      ierr = u->getInternalColumn(i-1,j,&u_im1); CHKERRQ(ierr);
      ierr = u->getInternalColumn(i+1,j,&u_ip1); CHKERRQ(ierr);

      ierr = v->getInternalColumn(i,j-1,&v_jm1); CHKERRQ(ierr);
      ierr = v->getInternalColumn(i,j+1,&v_jp1); CHKERRQ(ierr);

      if (bmr) {
        w[0] = - (*bmr)(i,j);
      } else {
        w[0] = 0.0;
      }

      PetscScalar OLDintegrand
             = (u_ip1[0] - u_im1[0]) / (2.0*dx) + (v_jp1[0] - v_jm1[0]) / (2.0*dy);
      for (PetscInt k = 1; k < grid.Mz; ++k) {
        const PetscScalar NEWintegrand
             = (u_ip1[k] - u_im1[k]) / (2.0*dx) + (v_jp1[k] - v_jm1[k]) / (2.0*dy);
        const PetscScalar dz = grid.zlevels[k] - grid.zlevels[k-1];
        w[k] = w[k-1] - 0.5 * (NEWintegrand + OLDintegrand) * dz;
        OLDintegrand = NEWintegrand;
      }
    }
  }

  if (bmr) {
    ierr = bmr->end_access(); CHKERRQ(ierr);
  }

  ierr = u->end_access(); CHKERRQ(ierr);
  ierr = v->end_access(); CHKERRQ(ierr);
  ierr = result.end_access(); CHKERRQ(ierr);
  
  return 0;
}
