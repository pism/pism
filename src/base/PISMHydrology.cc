// Copyright (C) 2012 PISM Authors
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

#include "PISMHydrology.hh"
#include "pism_options.hh"

PetscErrorCode PISMHydrology::init_steady(IceModelVec2S W0) {
  PetscErrorCode ierr;

  ierr = W0.copy_to(W); CHKERRQ(ierr);

  thickness = dynamic_cast<IceModelVec2S*>(vars.get("land_ice_thickness"));
  if (thickness == NULL) SETERRQ(grid.com, 1, "land_ice_thickness is not available");
  // FIXME: same for Ubase

  ierr = thickness->begin_access(); CHKERRQ(ierr);
  ierr = Po.begin_access(); CHKERRQ(ierr);
  ierr = cbase.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      Po(i,j) = rhoi * g * (*thickness)(i,j);
      cbase(i,j) = sqrt((*Ubase)(i,j).u * (*Ubase)(i,j).u + (*Ubase)(i,j).v * (*Ubase)(i,j).v);
    }
  }
  ierr = thickness->end_access(); CHKERRQ(ierr);
  ierr = Po.end_access(); CHKERRQ(ierr);
  ierr = cbase.end_access(); CHKERRQ(ierr);

  ierr = W.begin_access(); CHKERRQ(ierr);
  ierr = P.begin_access(); CHKERRQ(ierr);
  ierr = Po.begin_access(); CHKERRQ(ierr);
  ierr = cbase.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      Ptmp = c1 * cbase(i,j) * PetscMax(0.0,Wr - W(i,j)) / (c2 * A);
      P(i,j) = PetscMax(0.0,Po(i,j) - pow(Ptmp,1.0/3.0));
    }
  }
  ierr = W.end_access(); CHKERRQ(ierr);
  ierr = P.end_access(); CHKERRQ(ierr);
  ierr = Po.end_access(); CHKERRQ(ierr);
  ierr = cbase.end_access(); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode PISMHydrology::update_water_and_pressure(PetscScalar dt) {
  PetscErrorCode ierr;

  ierr = W.beginGhostComm(); CHKERRQ(ierr);
  ierr = P.beginGhostComm(); CHKERRQ(ierr);
  ierr = W.endGhostComm(); CHKERRQ(ierr);
  ierr = P.endGhostComm(); CHKERRQ(ierr);

  ierr = get_V_components(); CHKERRQ(ierr);  // fills alph and beta

  PetscScalar **bwatnew; 
  ierr = vbwat.begin_access(); CHKERRQ(ierr);
  ierr = vWork2d[0].get_array(bwatnew); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      FIXME
      Wij = W(i,j);
      upe = up(alphV(i,j),   Wij,      W(i+1,j));
      upw = up(alphV(i-1,j), W(i-1,j), Wij);
      upn = up(betaV(i,j),   Wij,      W(i,j+1));
      ups = up(betaV(i,j-1), W(i,j-1), Wij);
      inputdepth = dt * Phi(i,j);
      dtlapW = mux * (Wea(i,j) * (W(i+1,j)-Wij) - Wea(i-1,j) * (Wij-W(i-1,j))) + ...
               muy * (Wno(i,j) * (W(i,j+1)-Wij) - Wno(i,j-1) * (Wij-W(i,j-1)));
      Wnew(i,j) = Wij - nux * (upe - upw) - nuy * (upn - ups) + dtlapW + ...
                  inputdepth;
    }
  }
  ierr = vWork2d[0].end_access(); CHKERRQ(ierr);
  ierr = vbwat.end_access(); CHKERRQ(ierr);

  // FIXME:  time step of P equation

  return 0;
}


PetscErrorCode get_water_thickness(IceModelVec2S &result) {
  PetscErrorCode ierr = W.copy_to(result); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode get_water_pressure(IceModelVec2S &result) {
  PetscErrorCode ierr = P.copy_to(result); CHKERRQ(ierr);
  return 0;
}

