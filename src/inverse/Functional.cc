// Copyright (C) 2012, 2013  David Maxwell
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

#include "Functional.hh"

PetscErrorCode gradientFD(Functional<IceModelVec2S> &f, IceModelVec2S &x, IceModelVec2S &gradient) {
  PetscErrorCode ierr;
  IceGrid &grid = *x.get_grid();
  PetscReal h = PETSC_SQRT_MACHINE_EPSILON;

  PetscReal F0,Fh;
  
  ierr = f.valueAt(x,&F0); CHKERRQ(ierr);
  
  ierr = gradient.begin_access(); CHKERRQ(ierr);
  for(PetscInt i=grid.xs; i< grid.xs+grid.xm; i++) {
    for(PetscInt j=grid.ys; j< grid.ys+grid.ym; j++) {
      x.begin_access(); CHKERRQ(ierr);
      x(i,j) += h;
      x.end_access(); CHKERRQ(ierr);
      x.update_ghosts();
      ierr = f.valueAt(x,&Fh); CHKERRQ(ierr);
      x.begin_access(); CHKERRQ(ierr);
      x(i,j) -= h;
      x.end_access(); CHKERRQ(ierr);
      x.update_ghosts();
      gradient(i,j) = (Fh-F0)/h;
    }
  }
  ierr = gradient.end_access(); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode gradientFD(Functional<IceModelVec2V> &f, IceModelVec2V &x, IceModelVec2V &gradient) {
  PetscErrorCode ierr;
  IceGrid &grid = *x.get_grid();
  PetscReal h = PETSC_SQRT_MACHINE_EPSILON;

  PetscReal F0,Fh;
  
  ierr = f.valueAt(x,&F0); CHKERRQ(ierr);
  
  ierr = gradient.begin_access(); CHKERRQ(ierr);
  for(PetscInt i=grid.xs; i< grid.xs+grid.xm; i++) {
    for(PetscInt j=grid.ys; j< grid.ys+grid.ym; j++) {
      x.begin_access(); CHKERRQ(ierr);
      x(i,j).u += h;
      x.end_access(); CHKERRQ(ierr);
      x.update_ghosts();
      ierr = f.valueAt(x,&Fh); CHKERRQ(ierr);
      x.begin_access(); CHKERRQ(ierr);
      x(i,j).u -= h;
      x.end_access(); CHKERRQ(ierr);
      x.update_ghosts();
      gradient(i,j).u = (Fh-F0)/h;

      x.begin_access(); CHKERRQ(ierr);
      x(i,j).v += h;
      x.end_access(); CHKERRQ(ierr);
      x.update_ghosts();
      ierr = f.valueAt(x,&Fh); CHKERRQ(ierr);
      x.begin_access(); CHKERRQ(ierr);
      x(i,j).v -= h;
      x.end_access(); CHKERRQ(ierr);
      x.update_ghosts();
      gradient(i,j).v = (Fh-F0)/h;
    }
  }
  ierr = gradient.end_access(); CHKERRQ(ierr);
  return 0;
}

// PetscErrorCode gradientFD(Functional<IceModelVec2V> &f, IceModelVec2V &x, IceModelVec2V &gradient);
