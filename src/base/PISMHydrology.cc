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
#include "PISMVars.hh"
#include "pism_options.hh"
#include "Mask.hh"


PISMHydrology::PISMHydrology(IceGrid &g, const NCConfigVariable &conf)
    : PISMComponent_TS(g, conf)
{
    bed   = NULL;
    thk   = NULL;
    surf  = NULL;
    Ubase = NULL;

    if (allocate() != 0) {
      PetscPrintf(grid.com, "PISM ERROR: memory allocation failed in PISMHydrology constructor.\n");
      PISMEnd();
    }

    ice_density = config.get("ice_density");
    standard_gravity = config.get("standard_gravity");
    fresh_water_density = config.get("fresh_water_density");

    // initialize using constants from van Pelt & Bueler preprint
    // FIXME: should be configurable
    c1    = 0.500;      // m-1
    c2    = 0.040;      // [pure]
    K     = 1.0e-2;     // m s-1;  want Kmax or Kmin according to W > Wr
    Aglen = 3.1689e-24; // Pa-3 s-1; ice softness
    nglen = 3.0;
    Wr    = 1.0;        // m
    E0    = 1.0;        // m; what is optimal?
    Y0    = 0.001;      // m; regularization

    c0    = K / (fresh_water_density * standard_gravity); // constant in velocity formula
}


PetscErrorCode PISMHydrology::allocate() {
  PetscErrorCode ierr;

  ierr = W.create(grid, "bwat", true, grid.max_stencil_width); CHKERRQ(ierr);
  ierr = W.set_attrs("model_state",
                     "thickness of subglacial water layer",
                     "m", ""); CHKERRQ(ierr);
  ierr = W.set_attr("valid_min", 0.0); CHKERRQ(ierr);

  ierr = P.create(grid, "bwp", true, grid.max_stencil_width); CHKERRQ(ierr);
  ierr = P.set_attrs("model_state",
                     "pressure of water in subglacial layer",
                     "Pa", ""); CHKERRQ(ierr);
  ierr = P.set_attr("valid_min", 0.0); CHKERRQ(ierr);

  ierr = Po.create(grid, "ice-overburden-pressure", false); CHKERRQ(ierr);
  ierr = Po.set_attrs("internal",
                      "ice overburden pressure seen by subglacial water layer",
                      "Pa", ""); CHKERRQ(ierr);
  ierr = Po.set_attr("valid_min", 0.0); CHKERRQ(ierr);

  ierr = psi.create(grid, "hydraulic-potential", true, grid.max_stencil_width); CHKERRQ(ierr);
  ierr = psi.set_attrs("internal",
                       "hydraulic potential of water in subglacial layer",
                       "Pa", ""); CHKERRQ(ierr);

  ierr = cbase.create(grid, "ice-sliding-speed", false); CHKERRQ(ierr);
  ierr = cbase.set_attrs("internal",
                         "ice sliding speed seen by subglacial water layer",
                         "m s-1", ""); CHKERRQ(ierr);
  ierr = cbase.set_attr("valid_min", 0.0); CHKERRQ(ierr);

  ierr = V.create(grid, "water-velocity", true, grid.max_stencil_width); CHKERRQ(ierr);
  ierr = V.set_attrs("internal",
                     "cell face-centered (staggered) components of water velocity in subglacial water layer",
                     "m s-1", ""); CHKERRQ(ierr);

  ierr = Wnew.create(grid, "Wnew-internal", false); CHKERRQ(ierr);
  ierr = Wnew.set_attrs("internal",
                     "new thickness of subglacial water layer before update",
                     "m", ""); CHKERRQ(ierr);
  ierr = Wnew.set_attr("valid_min", 0.0); CHKERRQ(ierr);

  ierr = Pnew.create(grid, "Pnew-internal", false); CHKERRQ(ierr);
  ierr = Pnew.set_attrs("internal",
                     "new subglacial water pressure before update",
                     "Pa", ""); CHKERRQ(ierr);
  ierr = Pnew.set_attr("valid_min", 0.0); CHKERRQ(ierr);

  return 0;
}


PetscErrorCode PISMHydrology::init(PISMVars &vars) {
  PetscErrorCode ierr;

  variables = &vars;

  ierr = verbPrintf(2, grid.com,
    "* Initializing the vanPelt-Bueler subglacial hydrology model...\n"); CHKERRQ(ierr);

  bed = dynamic_cast<IceModelVec2S*>(vars.get("topg"));
  if (bed == NULL) SETERRQ(grid.com, 1, "topg is not available");

  thk = dynamic_cast<IceModelVec2S*>(vars.get("thk"));
  if (thk == NULL) SETERRQ(grid.com, 1, "thk is not available");

  surf = dynamic_cast<IceModelVec2S*>(vars.get("usurf"));
  if (surf == NULL) SETERRQ(grid.com, 1, "usurf is not available");

  Ubase = dynamic_cast<IceModelVec2V*>(vars.get("Ubase-NOTIONAL"));
  if (Ubase == NULL) SETERRQ(grid.com, 1, "Ubase-NOTIONAL is not available ... IT DOES NOT EXIST");

  // get the water layer thickness from the context if present
  IceModelVec2S *W_input = dynamic_cast<IceModelVec2S*>(vars.get("bwat"));
  if (W_input != NULL) {
    ierr = W.copy_from(*W_input); CHKERRQ(ierr);
  } else {
    ierr = W.set(0.0); CHKERRQ(ierr);
  }

  // get the water pressure from the context if present
  IceModelVec2S *P_input = dynamic_cast<IceModelVec2S*>(vars.get("bwp"));
  if (P_input != NULL) {
    ierr = P.copy_from(*P_input); CHKERRQ(ierr);
  } else {
    ierr = P_from_W_steady(P); CHKERRQ(ierr);
  }

  return 0;
}


PetscErrorCode PISMHydrology::P_from_W_steady(IceModelVec2S &result) {
  PetscErrorCode ierr;
  PetscReal CC = c1 / (c2 * Aglen),
            powglen = 1.0/nglen,
            sb, Wratio;

  ierr = W.begin_access(); CHKERRQ(ierr);
  ierr = Po.begin_access(); CHKERRQ(ierr);
  ierr = cbase.begin_access(); CHKERRQ(ierr);
  ierr = result.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      sb     = pow(CC * cbase(i,j),powglen);
      Wratio = PetscMax(0.0,Wr - W(i,j)) / (W(i,j) + Y0);
      // in cases where steady state is actually possible this will
      //   come out positive, but otherwise we should get underpressure P=0,
      //   and that is what it yields
      result(i,j) = PetscMax( 0.0,Po(i,j) - sb * pow(Wratio,powglen) );
    }
  }
  ierr = W.end_access(); CHKERRQ(ierr);
  ierr = Po.end_access(); CHKERRQ(ierr);
  ierr = cbase.end_access(); CHKERRQ(ierr);
  ierr = result.end_access(); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode PISMHydrology::V_components(IceModelVec2Stag &result) {
  PetscErrorCode ierr;
  PetscReal dbdx, dbdy, dPdx, dPdy;

  ierr = P.begin_access(); CHKERRQ(ierr);
  ierr = bed->begin_access(); CHKERRQ(ierr);
  ierr = result.begin_access(); CHKERRQ(ierr);
  for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
      dbdx = ((*bed)(i+1,j) - (*bed)(i,j)) / grid.dx;
      dbdy = ((*bed)(i,j+1) - (*bed)(i,j)) / grid.dy;
      dPdx = (P(i+1,j) - P(i,j)) / grid.dx;
      dPdy = (P(i,j+1) - P(i,j)) / grid.dy;
      result(i,j,0) = - c0 * dPdx - K * dbdx;
      result(i,j,1) = - c0 * dPdy - K * dbdy;
    }
  }
  ierr = P.end_access(); CHKERRQ(ierr);
  ierr = bed->end_access(); CHKERRQ(ierr);
  ierr = result.end_access(); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode PISMHydrology::update_ice_functions(IceModelVec2S &result_Po,
                                                   IceModelVec2S &result_cbase) {
  PetscErrorCode ierr;

  ierr = thk->begin_access(); CHKERRQ(ierr);
  ierr = Ubase->begin_access(); CHKERRQ(ierr);
  ierr = result_Po.begin_access(); CHKERRQ(ierr);
  ierr = result_cbase.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      result_Po(i,j) = ice_density * standard_gravity * (*thk)(i,j);
      result_cbase(i,j) = sqrt((*Ubase)(i,j).u * (*Ubase)(i,j).u + (*Ubase)(i,j).v * (*Ubase)(i,j).v);
    }
  }
  ierr = thk->end_access(); CHKERRQ(ierr);
  ierr = Ubase->end_access(); CHKERRQ(ierr);
  ierr = result_Po.end_access(); CHKERRQ(ierr);
  ierr = result_cbase.end_access(); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode PISMHydrology::update(PetscReal my_t, PetscReal my_dt) {
  PetscErrorCode ierr;

//FIXME: this first version is totally incorrect; see doublediff.m

  ierr = W.beginGhostComm(); CHKERRQ(ierr);
  ierr = P.beginGhostComm(); CHKERRQ(ierr);
  ierr = W.endGhostComm(); CHKERRQ(ierr);
  ierr = P.endGhostComm(); CHKERRQ(ierr);

  // from current ice geometry/velocity variables, initialize Po and cbase
  ierr = update_ice_functions(Po,cbase); CHKERRQ(ierr);

  for (PetscInt m=1; m<10; m++) { // FIXME: fake hydrology time-stepping loop
    ierr = V_components(V); CHKERRQ(ierr);  // fills alph and beta

    PetscScalar Wij;
    ierr = W.begin_access(); CHKERRQ(ierr);
    ierr = Wnew.begin_access(); CHKERRQ(ierr);
    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
        Wij = W(i,j);
        //inputdepth = dt * Phi(i,j);
        //dtlapW = mux * (Wea(i,j) * (W(i+1,j)-Wij) - Wea(i-1,j) * (Wij-W(i-1,j))) + ...
        //         muy * (Wno(i,j) * (W(i,j+1)-Wij) - Wno(i,j-1) * (Wij-W(i,j-1)));
        //Wnew(i,j) = Wij - FIXME + dtlapW + inputdepth;
        Wnew(i,j) = Wij;
      }
    }
    ierr = Wnew.end_access(); CHKERRQ(ierr);
    ierr = W.end_access(); CHKERRQ(ierr);

    // FIXME:  time step of P equation

  } // end of hydrology model time-stepping loop

  return 0;
}


PetscErrorCode PISMHydrology::water_layer_thickness(IceModelVec2S &result) {
  PetscErrorCode ierr = W.copy_to(result); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode PISMHydrology::water_pressure(IceModelVec2S &result) {
  PetscErrorCode ierr = P.copy_to(result); CHKERRQ(ierr);
  return 0;
}

