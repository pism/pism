// Copyright (C) 2012, 2013, 2014, 2015, 2017  David Maxwell
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
#include "pism/util/IceGrid.hh"
#include "pism/util/error_handling.hh"

namespace pism {
namespace inverse {

void gradientFD(IPFunctional<IceModelVec2S> &f, IceModelVec2S &x, IceModelVec2S &gradient) {
  const IceGrid &grid = *x.grid();
  double h = PETSC_SQRT_MACHINE_EPSILON;

  double F0,Fh;
  
  f.valueAt(x,&F0);
  
  IceModelVec::AccessList list(gradient);

  ParallelSection loop(grid.com);
  try {
    for (Points p(grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      {
        IceModelVec::AccessList access_x(x);
        x(i,j) += h;
      }
      x.update_ghosts();

      f.valueAt(x,&Fh);

      {
        IceModelVec::AccessList access_x(x);
        x(i,j) -= h;
      }
      x.update_ghosts();

      gradient(i,j) = (Fh-F0)/h;
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();
}

void gradientFD(IPFunctional<IceModelVec2V> &f, IceModelVec2V &x, IceModelVec2V &gradient) {
  const IceGrid &grid = *x.grid();
  double h = PETSC_SQRT_MACHINE_EPSILON;

  double F0,Fh;
  
  f.valueAt(x,&F0);
  
  IceModelVec::AccessList access_gradient(gradient);

  ParallelSection loop(grid.com);
  try {
    for (Points p(grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      {
        IceModelVec::AccessList access_x(x);
        x(i,j).u += h;
      }
      x.update_ghosts();

      f.valueAt(x,&Fh);

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

      f.valueAt(x,&Fh);

      {
        IceModelVec::AccessList access_x(x);
        x(i,j).v -= h;
      }
      x.update_ghosts();

      gradient(i,j).v = (Fh-F0)/h;
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();
}

// PetscErrorCode gradientFD(IPFunctional<IceModelVec2V> &f, IceModelVec2V &x, IceModelVec2V &gradient);

} // end of namespace inverse
} // end of namespace pism
