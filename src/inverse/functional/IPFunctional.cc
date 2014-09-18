// Copyright (C) 2012, 2013, 2014  David Maxwell
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

#include "IPFunctional.hh"

namespace pism {

PetscErrorCode gradientFD(IPFunctional<IceModelVec2S> &f, IceModelVec2S &x, IceModelVec2S &gradient) {
  PetscErrorCode ierr;
  IceGrid &grid = *x.get_grid();
  double h = PETSC_SQRT_MACHINE_EPSILON;

  double F0,Fh;
  
  ierr = f.valueAt(x,&F0); CHKERRQ(ierr);
  
  IceModelVec::AccessList list(gradient);

  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    {
      IceModelVec::AccessList access_x(x);
      x(i,j) += h;
    }
    x.update_ghosts();

    ierr = f.valueAt(x,&Fh); CHKERRQ(ierr);

    {
      IceModelVec::AccessList access_x(x);
      x(i,j) -= h;
    }
    x.update_ghosts();

    gradient(i,j) = (Fh-F0)/h;
  }
  return 0;
}

PetscErrorCode gradientFD(IPFunctional<IceModelVec2V> &f, IceModelVec2V &x, IceModelVec2V &gradient) {
  PetscErrorCode ierr;
  IceGrid &grid = *x.get_grid();
  double h = PETSC_SQRT_MACHINE_EPSILON;

  double F0,Fh;
  
  ierr = f.valueAt(x,&F0); CHKERRQ(ierr);
  
  IceModelVec::AccessList access_gradient(gradient);

  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    {
      IceModelVec::AccessList access_x(x);
      x(i,j).u += h;
    }
    x.update_ghosts();

    ierr = f.valueAt(x,&Fh); CHKERRQ(ierr);

    {
      IceModelVec::AccessList access_x(x);
      x(i,j).u -= h;
    }
    x.update_ghosts();

    gradient(i,j).u = (Fh-F0)/h;

    {
      IceModelVec::AccessList access_x(x);
      x(i,j).v += h;
    }
    x.update_ghosts();

    ierr = f.valueAt(x,&Fh); CHKERRQ(ierr);

    {
      IceModelVec::AccessList access_x(x);
      x(i,j).v -= h;
    }
    x.update_ghosts();

    gradient(i,j).v = (Fh-F0)/h;
  }
  return 0;
}

// PetscErrorCode gradientFD(IPFunctional<IceModelVec2V> &f, IceModelVec2V &x, IceModelVec2V &gradient);

} // end of namespace pism
