// Copyright (C) 2012  David Maxwell
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

#include "LogRelativeFunctional.hh"

PetscErrorCode LogRelativeFunctional::normalize(PetscReal scale) {
  PetscErrorCode   ierr;

  // The local value of the weights
  PetscReal value = 0;

  for( PetscInt i=m_grid.xs; i<m_grid.xs+m_grid.xm; i++) {
    for( PetscInt j=m_grid.ys; j<m_grid.ys+m_grid.ym; j++) {
      value += 1;
    }
  }
  
  ierr = PISMGlobalSum(&value, &m_normalization, m_grid.com); CHKERRQ(ierr);
  m_normalization *= scale;
  return 0;
}

PetscErrorCode LogRelativeFunctional::valueAt(IceModelVec2V &x, PetscReal *OUTPUT)  {
  PetscErrorCode   ierr;

  // The value of the objective
  PetscReal value = 0;

  PISMVector2 **x_a;
  ierr = x.get_array(x_a); CHKERRQ(ierr);

  PetscReal v_eps = m_grid.config.get("inv_ssa_velocity_scale")*1e-4;

  PISMVector2 **u_obs_a;
  ierr = m_u_observed.get_array(u_obs_a); CHKERRQ(ierr);
  for( PetscInt i=m_grid.xs; i<m_grid.xs+m_grid.xm; i++) {
    for( PetscInt j=m_grid.ys; j<m_grid.ys+m_grid.ym; j++) {
      PISMVector2 &x_ij = x_a[i][j];
      PISMVector2 &u_obs_ij = u_obs_a[i][j];
      PetscReal obsMagSq = u_obs_ij.u*u_obs_ij.u + u_obs_ij.v*u_obs_ij.v + v_eps*v_eps;
      value += log( 1 + (x_ij.u*x_ij.u + x_ij.v*x_ij.v)/obsMagSq);
    }
  }
  ierr = m_u_observed.end_access(); CHKERRQ(ierr);
  
  value /= m_normalization;
  
  ierr = PISMGlobalSum(&value, OUTPUT, m_grid.com); CHKERRQ(ierr);

  ierr = x.end_access(); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode LogRelativeFunctional::gradientAt(IceModelVec2V &x, IceModelVec2V &gradient)  {
  PetscErrorCode   ierr;

  gradient.set(0);

  PISMVector2 **x_a;
  ierr = x.get_array(x_a); CHKERRQ(ierr);

  PISMVector2 **gradient_a;
  ierr = gradient.get_array(gradient_a); CHKERRQ(ierr);

  PetscReal v_eps = m_grid.config.get("inv_ssa_velocity_scale")*1e-4;

  PISMVector2 **u_obs_a;
  ierr = m_u_observed.get_array(u_obs_a); CHKERRQ(ierr);
  for( PetscInt i=m_grid.xs; i<m_grid.xs+m_grid.xm; i++) {
    for( PetscInt j=m_grid.ys; j<m_grid.ys+m_grid.ym; j++) {
      PISMVector2 &x_ij = x_a[i][j];
      PISMVector2 &u_obs_ij = u_obs_a[i][j];
      PetscReal obsMagSq = u_obs_ij.u*u_obs_ij.u + u_obs_ij.v*u_obs_ij.v + v_eps*v_eps;
      PetscReal dJdu =  1./( obsMagSq + x_ij.u*x_ij.u + x_ij.v*x_ij.v);

      gradient_a[i][j].u = dJdu*2*x_a[i][j].u/m_normalization;
      gradient_a[i][j].v = dJdu*2*x_a[i][j].v/m_normalization;
    }
  }
  ierr = m_u_observed.end_access(); CHKERRQ(ierr);

  ierr = x.end_access(); CHKERRQ(ierr);
  ierr = gradient.end_access(); CHKERRQ(ierr);

  return 0;
}
