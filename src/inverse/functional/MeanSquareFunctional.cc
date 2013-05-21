// Copyright (C) 2012  David Maxwell
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

#include "MeanSquareFunctional.hh"

PetscErrorCode MeanSquareFunctional2V::normalize(PetscReal scale) {
  PetscErrorCode   ierr;

  // The local value of the weights
  PetscReal value = 0;

  if(m_weights) {
    PetscReal **w_a;
    ierr = m_weights->get_array(w_a); CHKERRQ(ierr);
    for( PetscInt i=m_grid.xs; i<m_grid.xs+m_grid.xm; i++) {
      for( PetscInt j=m_grid.ys; j<m_grid.ys+m_grid.ym; j++) {
        value += w_a[i][j];
      }
    }
    ierr = m_weights->end_access(); CHKERRQ(ierr);
  } else {
    for( PetscInt i=m_grid.xs; i<m_grid.xs+m_grid.xm; i++) {
      for( PetscInt j=m_grid.ys; j<m_grid.ys+m_grid.ym; j++) {
        value += 1;
      }
    }
  }
  
  ierr = PISMGlobalSum(&value, &m_normalization, m_grid.com); CHKERRQ(ierr);
  m_normalization *= (scale*scale);
  return 0;
}

PetscErrorCode MeanSquareFunctional2V::valueAt(IceModelVec2V &x, PetscReal *OUTPUT)  {
  PetscErrorCode   ierr;

  // The value of the objective
  PetscReal value = 0;

  PISMVector2 **x_a;
  ierr = x.get_array(x_a); CHKERRQ(ierr);

  if(m_weights) {
    PetscReal **w_a;
    ierr = m_weights->get_array(w_a); CHKERRQ(ierr);
    for( PetscInt i=m_grid.xs; i<m_grid.xs+m_grid.xm; i++) {
      for( PetscInt j=m_grid.ys; j<m_grid.ys+m_grid.ym; j++) {
        PISMVector2 &x_ij = x_a[i][j];
        value += (x_ij.u*x_ij.u + x_ij.v*x_ij.v)*w_a[i][j];
      }
    }
    ierr = m_weights->end_access(); CHKERRQ(ierr);
  } else {
    for( PetscInt i=m_grid.xs; i<m_grid.xs+m_grid.xm; i++) {
      for( PetscInt j=m_grid.ys; j<m_grid.ys+m_grid.ym; j++) {
        PISMVector2 &x_ij = x_a[i][j];
        value += (x_ij.u*x_ij.u + x_ij.v*x_ij.v);
      }
    }
  }
  value /= m_normalization;
  
  ierr = PISMGlobalSum(&value, OUTPUT, m_grid.com); CHKERRQ(ierr);

  ierr = x.end_access(); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode MeanSquareFunctional2V::dot(IceModelVec2V &a, IceModelVec2V &b, PetscReal *OUTPUT)  {
  PetscErrorCode   ierr;

  // The value of the objective
  PetscReal value = 0;

  PISMVector2 **a_a;
  ierr = a.get_array(a_a); CHKERRQ(ierr);

  PISMVector2 **b_a;
  ierr = b.get_array(b_a); CHKERRQ(ierr);

  if(m_weights) {
    PetscReal **w_a;
    ierr = m_weights->get_array(w_a); CHKERRQ(ierr);
    for( PetscInt i=m_grid.xs; i<m_grid.xs+m_grid.xm; i++) {
      for( PetscInt j=m_grid.ys; j<m_grid.ys+m_grid.ym; j++) {
        PISMVector2 &a_ij = a_a[i][j];
        PISMVector2 &b_ij = b_a[i][j];
        value += (a_ij.u*b_ij.u + a_ij.v*b_ij.v)*w_a[i][j];
      }
    }
    ierr = m_weights->end_access(); CHKERRQ(ierr);
  } else {
    for( PetscInt i=m_grid.xs; i<m_grid.xs+m_grid.xm; i++) {
      for( PetscInt j=m_grid.ys; j<m_grid.ys+m_grid.ym; j++) {
        PISMVector2 &a_ij = a_a[i][j];
        PISMVector2 &b_ij = b_a[i][j];
        value += (a_ij.u*b_ij.u + a_ij.v*b_ij.v);
      }
    }
  }
  value /= m_normalization;
  
  ierr = PISMGlobalSum(&value, OUTPUT, m_grid.com); CHKERRQ(ierr);

  ierr = a.end_access(); CHKERRQ(ierr);
  ierr = b.end_access(); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode MeanSquareFunctional2V::gradientAt(IceModelVec2V &x, IceModelVec2V &gradient)  {
  PetscErrorCode   ierr;

  gradient.set(0);

  PISMVector2 **x_a;
  ierr = x.get_array(x_a); CHKERRQ(ierr);

  PISMVector2 **gradient_a;
  ierr = gradient.get_array(gradient_a); CHKERRQ(ierr);

  if(m_weights) {
    PetscReal **w_a;
    ierr = m_weights->get_array(w_a); CHKERRQ(ierr);
    for( PetscInt i=m_grid.xs; i<m_grid.xs+m_grid.xm; i++) {
      for( PetscInt j=m_grid.ys; j<m_grid.ys+m_grid.ym; j++) {
        gradient_a[i][j].u = 2*x_a[i][j].u*w_a[i][j]/m_normalization;
        gradient_a[i][j].v = 2*x_a[i][j].v*w_a[i][j]/m_normalization;
      }
    }
    ierr = m_weights->end_access(); CHKERRQ(ierr);
  } else {
    for( PetscInt i=m_grid.xs; i<m_grid.xs+m_grid.xm; i++) {
      for( PetscInt j=m_grid.ys; j<m_grid.ys+m_grid.ym; j++) {
        gradient_a[i][j].u = 2*x_a[i][j].u/m_normalization;
        gradient_a[i][j].v = 2*x_a[i][j].v/m_normalization;
      }
    }
  }

  ierr = x.end_access(); CHKERRQ(ierr);
  ierr = gradient.end_access(); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode MeanSquareFunctional2S::normalize(PetscReal scale) {
  PetscErrorCode   ierr;

  // The local value of the weights
  PetscReal value = 0;

  if(m_weights) {
    PetscReal **w_a;
    ierr = m_weights->get_array(w_a); CHKERRQ(ierr);
    for( PetscInt i=m_grid.xs; i<m_grid.xs+m_grid.xm; i++) {
      for( PetscInt j=m_grid.ys; j<m_grid.ys+m_grid.ym; j++) {
        value += w_a[i][j];
      }
    }
    ierr = m_weights->end_access(); CHKERRQ(ierr);
  } else {
    for( PetscInt i=m_grid.xs; i<m_grid.xs+m_grid.xm; i++) {
      for( PetscInt j=m_grid.ys; j<m_grid.ys+m_grid.ym; j++) {
        value += 1;
      }
    }
  }
  
  ierr = PISMGlobalSum(&value, &m_normalization, m_grid.com); CHKERRQ(ierr);
  m_normalization *= (scale*scale);
  return 0;
}

PetscErrorCode MeanSquareFunctional2S::valueAt(IceModelVec2S &x, PetscReal *OUTPUT)  {
  PetscErrorCode   ierr;

  // The value of the objective
  PetscReal value = 0;

  PetscReal **x_a;
  ierr = x.get_array(x_a); CHKERRQ(ierr);

  if(m_weights) {
    PetscReal **w_a;
    ierr = m_weights->get_array(w_a); CHKERRQ(ierr);
    for( PetscInt i=m_grid.xs; i<m_grid.xs+m_grid.xm; i++) {
      for( PetscInt j=m_grid.ys; j<m_grid.ys+m_grid.ym; j++) {
        PetscReal &x_ij = x_a[i][j];
        value += x_ij*x_ij*w_a[i][j];
      }
    }
    ierr = m_weights->end_access(); CHKERRQ(ierr);
  } else {
    for( PetscInt i=m_grid.xs; i<m_grid.xs+m_grid.xm; i++) {
      for( PetscInt j=m_grid.ys; j<m_grid.ys+m_grid.ym; j++) {
        PetscReal &x_ij = x_a[i][j];
        value += x_ij*x_ij;
      }
    }
  }
  value /= m_normalization;
  
  ierr = PISMGlobalSum(&value, OUTPUT, m_grid.com); CHKERRQ(ierr);

  ierr = x.end_access(); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode MeanSquareFunctional2S::dot(IceModelVec2S &a, IceModelVec2S &b, PetscReal *OUTPUT)  {
  PetscErrorCode   ierr;

  // The value of the objective
  PetscReal value = 0;

  PetscReal **a_a;
  ierr = a.get_array(a_a); CHKERRQ(ierr);

  PetscReal **b_a;
  ierr = b.get_array(b_a); CHKERRQ(ierr);

  if(m_weights) {
    PetscReal **w_a;
    ierr = m_weights->get_array(w_a); CHKERRQ(ierr);
    for( PetscInt i=m_grid.xs; i<m_grid.xs+m_grid.xm; i++) {
      for( PetscInt j=m_grid.ys; j<m_grid.ys+m_grid.ym; j++) {
        value += (a_a[i][j]*b_a[i][j])*w_a[i][j];
      }
    }
    ierr = m_weights->end_access(); CHKERRQ(ierr);
  } else {
    for( PetscInt i=m_grid.xs; i<m_grid.xs+m_grid.xm; i++) {
      for( PetscInt j=m_grid.ys; j<m_grid.ys+m_grid.ym; j++) {
        value += (a_a[i][j]*b_a[i][j]);
      }
    }
  }
  value /= m_normalization;
  
  ierr = PISMGlobalSum(&value, OUTPUT, m_grid.com); CHKERRQ(ierr);

  ierr = a.end_access(); CHKERRQ(ierr);
  ierr = b.end_access(); CHKERRQ(ierr);

  return 0;
}


PetscErrorCode MeanSquareFunctional2S::gradientAt(IceModelVec2S &x, IceModelVec2S &gradient)  {
  PetscErrorCode   ierr;

  gradient.set(0);

  PetscReal **x_a;
  ierr = x.get_array(x_a); CHKERRQ(ierr);

  PetscReal **gradient_a;
  ierr = gradient.get_array(gradient_a); CHKERRQ(ierr);

  if(m_weights) {
    PetscReal **w_a;
    ierr = m_weights->get_array(w_a); CHKERRQ(ierr);
    for( PetscInt i=m_grid.xs; i<m_grid.xs+m_grid.xm; i++) {
      for( PetscInt j=m_grid.ys; j<m_grid.ys+m_grid.ym; j++) {
        gradient_a[i][j] = 2*x_a[i][j]*w_a[i][j]/m_normalization;
      }
    }
    ierr = m_weights->end_access(); CHKERRQ(ierr);
  } else {
    for( PetscInt i=m_grid.xs; i<m_grid.xs+m_grid.xm; i++) {
      for( PetscInt j=m_grid.ys; j<m_grid.ys+m_grid.ym; j++) {
        gradient_a[i][j] = 2*x_a[i][j]/m_normalization;
      }
    }
  }

  ierr = x.end_access(); CHKERRQ(ierr);
  ierr = gradient.end_access(); CHKERRQ(ierr);

  return 0;
}
