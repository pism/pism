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
    usurf  = NULL;
    bmelt  = NULL;
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

  ierr = cbase.create(grid, "ice-sliding-speed", false); CHKERRQ(ierr);
  ierr = cbase.set_attrs("internal",
                         "ice sliding speed seen by subglacial water layer",
                         "m s-1", ""); CHKERRQ(ierr);
  ierr = cbase.set_attr("valid_min", 0.0); CHKERRQ(ierr);

  ierr = psi.create(grid, "hydraulic-potential", true, grid.max_stencil_width); CHKERRQ(ierr);
  ierr = psi.set_attrs("internal",
                       "hydraulic potential of water in subglacial layer",
                       "Pa", ""); CHKERRQ(ierr);

  ierr = known.create(grid, "known-hydro-mask", true, grid.max_stencil_width); CHKERRQ(ierr);
  ierr = known.set_attrs("internal",
                       "mask for where subglacial hydrology state is known",
                       "", ""); CHKERRQ(ierr);

  ierr = V.create(grid, "water-velocity", true, grid.max_stencil_width); CHKERRQ(ierr);
  ierr = V.set_attrs("internal",
                     "cell face-centered (staggered) components of water velocity in subglacial water layer",
                     "m s-1", ""); CHKERRQ(ierr);

  ierr = Wstag.create(grid, "W-staggered", true, grid.max_stencil_width); CHKERRQ(ierr);
  ierr = Wstag.set_attrs("internal",
                     "cell face-centered (staggered) values of water layer thickness",
                     "m", ""); CHKERRQ(ierr);

  ierr = Qstag.create(grid, "advection-flux", true, grid.max_stencil_width); CHKERRQ(ierr);
  ierr = Qstag.set_attrs("internal",
                     "cell face-centered (staggered) components of advective subglacial water flux",
                     "m2 s-1", ""); CHKERRQ(ierr);

  ierr = Wnew.create(grid, "Wnew-internal", false); CHKERRQ(ierr);
  ierr = Wnew.set_attrs("internal",
                     "new thickness of subglacial water layer during update",
                     "m", ""); CHKERRQ(ierr);
  ierr = Wnew.set_attr("valid_min", 0.0); CHKERRQ(ierr);

  ierr = Pnew.create(grid, "Pnew-internal", false); CHKERRQ(ierr);
  ierr = Pnew.set_attrs("internal",
                     "new subglacial water pressure during update",
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

  usurf = dynamic_cast<IceModelVec2S*>(vars.get("usurf"));
  if (usurf == NULL) SETERRQ(grid.com, 1, "usurf is not available");

  bmelt = dynamic_cast<IceModelVec2S*>(vars.get("bmelt"));
  if (bmelt == NULL) SETERRQ(grid.com, 1, "bmelt is not available");

  Ubase = dynamic_cast<IceModelVec2V*>(vars.get("Ubase-NOTIONAL"));
  if (Ubase == NULL) SETERRQ(grid.com, 1, "Ubase-NOTIONAL is not available ... IT DOES NOT EXIST");

  // initialize water layer thickness from the context if present, otherwise zero
  IceModelVec2S *W_input = dynamic_cast<IceModelVec2S*>(vars.get("bwat"));
  if (W_input != NULL) {
    ierr = W.copy_from(*W_input); CHKERRQ(ierr);
  } else {
    ierr = W.set(0.0); CHKERRQ(ierr);
  }

  // initialize the water pressure from the context if present, otherwise steady P(W)
  IceModelVec2S *P_input = dynamic_cast<IceModelVec2S*>(vars.get("bwp"));
  if (P_input != NULL) {
    ierr = P.copy_from(*P_input); CHKERRQ(ierr);
  } else {
    ierr = P_from_W_steady(P); CHKERRQ(ierr);
  }

  return 0;
}


//! Check bounds on W and P and fail with message if not satisfied.
/*!
Checks \f$0 \le W\f$ and \f$0 \le P \le P_o\f$.
 */
PetscErrorCode PISMHydrology::check_bounds() {
  PetscErrorCode ierr;

  ierr = W.begin_access(); CHKERRQ(ierr);
  ierr = P.begin_access(); CHKERRQ(ierr);
  ierr = Po.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (W(i,j) < 0.0) {
        PetscPrintf(grid.com,
           "PISM ERROR: negative subglacial water layer thickness W(i,j) = %.3f m < 0.0 at (i,j)=(%d,%d)\n"
           "ENDING ... \n\n", W(i,j),i,j);
        PISMEnd();
      }
      if (P(i,j) < 0.0) {
        PetscPrintf(grid.com,
           "PISM ERROR: negative subglacial water pressure P(i,j) = %.3f Pa < 0.0 at (i,j)=(%d,%d)\n"
           "ENDING ... \n\n", P(i,j),i,j);
        PISMEnd();
      }
      if (P(i,j) > Po(i,j)) {
        PetscPrintf(grid.com,
           "PISM ERROR: subglacial water pressure P(i,j) = %.3f Pa exceeds overburden pressure\n"
           "            Po(i,j) = %.3f Pa at (i,j)=(%d,%d)\n"
           "ENDING ... \n\n", P(i,j),Po(i,j),i,j);
        PISMEnd();
      }
    }
  }
  ierr = W.end_access(); CHKERRQ(ierr);
  ierr = P.end_access(); CHKERRQ(ierr);
  ierr = Po.end_access(); CHKERRQ(ierr);
  return 0;
}


//! Compute functional relationship P(W) which applies only in steady state.
/*!
This will be used in initialization when P is otherwise unknown, and
in verification and/or reporting.  It is not used during time-dependent
model runs.  To be more complete, \f$P=P(W,P_o,|v_b|)\f$.
 */
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


//! Get the advection velocity V at the center of cell edges.
/*!
Computes the advection velocity \f$V=V(\nabla P,\nabla b)\f$ on the
staggered (face-centered) grid.  If V = (alpha,beta) in components
then we have <code> result(i,j,0) = alpha(i+1/2,j) </code> and
<code> result(i,j,1) = beta(i,j+1/2) </code>
 */
PetscErrorCode PISMHydrology::velocity_staggered(IceModelVec2Stag &result) {
  PetscErrorCode ierr;
  PetscReal dbdx, dbdy, dPdx, dPdy;

  ierr = P.begin_access(); CHKERRQ(ierr);
  ierr = bed->begin_access(); CHKERRQ(ierr);
  ierr = result.begin_access(); CHKERRQ(ierr);

SETERRQ(grid.com,1,"not implemented");
  // FIXME because this will be evaluated at i-1,j and i,j-1, compute as GHOSTS=1 ??
  for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
      dbdx = ((*bed)(i+1,j  ) - (*bed)(i,j)) / grid.dx;
      dbdy = ((*bed)(i  ,j+1) - (*bed)(i,j)) / grid.dy;
      dPdx = (P(i+1,j  ) - P(i,j)) / grid.dx;
      dPdy = (P(i  ,j+1) - P(i,j)) / grid.dy;
      result(i,j,0) = - c0 * dPdx - K * dbdx;
      result(i,j,1) = - c0 * dPdy - K * dbdy;
    }
  }
  ierr = P.end_access(); CHKERRQ(ierr);
  ierr = bed->end_access(); CHKERRQ(ierr);
  ierr = result.end_access(); CHKERRQ(ierr);
  return 0;
}


//! Average the regular grid water thickness to values at the center of cell edges.
PetscErrorCode PISMHydrology::water_thickness_staggered(IceModelVec2Stag &result) {
  PetscErrorCode ierr;
  ierr = W.begin_access(); CHKERRQ(ierr);
  ierr = result.begin_access(); CHKERRQ(ierr);

SETERRQ(grid.com,1,"not implemented");
  // FIXME because this will be evaluated at i-1,j and i,j-1, compute as GHOSTS=1 ??
  for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
      result(i,j,0) = 0.5 * (W(i,j) + W(i+1,j  ));
      result(i,j,1) = 0.5 * (W(i,j) + W(i  ,j+1));
    }
  }
  ierr = W.end_access(); CHKERRQ(ierr);
  ierr = result.end_access(); CHKERRQ(ierr);
  return 0;
}


//! Compute Q = V W at edge-centers (staggered grid) by first-order upwinding.
PetscErrorCode PISMHydrology::advective_fluxes(IceModelVec2Stag &result) {
  PetscErrorCode ierr;
  ierr = W.begin_access(); CHKERRQ(ierr);
  ierr = V.begin_access(); CHKERRQ(ierr);
  ierr = result.begin_access(); CHKERRQ(ierr);

SETERRQ(grid.com,1,"not implemented");
  // FIXME because this will be evaluated at i-1,j and i,j-1, compute as GHOSTS=1 ??
  for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
      result(i,j,0) = (V(i,j,0) >= 0.0) ? V(i,j,0) * W(i,j) :  V(i,j,0) * W(i+1,j  );
      result(i,j,1) = (V(i,j,1) >= 0.0) ? V(i,j,1) * W(i,j) :  V(i,j,1) * W(i,  j+1);
    }
  }
  ierr = W.end_access(); CHKERRQ(ierr);
  ierr = V.end_access(); CHKERRQ(ierr);
  ierr = result.end_access(); CHKERRQ(ierr);
  return 0;
}


//! Get the hydraulic potential from bedrock topography and current state variables.
/*!
Computes \f$\psi = P + \rho_w g (b + W)\f$.
 */
PetscErrorCode PISMHydrology::hydraulic_potential(IceModelVec2S &result) {
  PetscErrorCode ierr;
  ierr = W.begin_access(); CHKERRQ(ierr);
  ierr = P.begin_access(); CHKERRQ(ierr);
  ierr = bed->begin_access(); CHKERRQ(ierr);
  ierr = result.begin_access(); CHKERRQ(ierr);

SETERRQ(grid.com,1,"not implemented");
  // FIXME because this will be evaluated at i+1,j and i-1,j and i,j-1 and i,j+1 compute as GHOSTS=1 ??
  for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
      result(i,j) = P(i,j) + fresh_water_density * standard_gravity * ((*bed)(i,j) + W(i,j));
    }
  }
  ierr = W.end_access(); CHKERRQ(ierr);
  ierr = P.end_access(); CHKERRQ(ierr);
  ierr = bed->end_access(); CHKERRQ(ierr);
  ierr = result.end_access(); CHKERRQ(ierr);
  return 0;
}


//! Update the overburden pressure Po and the sliding speed |v_b| from quantities in IceModel.
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
      result_cbase(i,j) = sqrt((*Ubase)(i,j).u * (*Ubase)(i,j).u
                                + (*Ubase)(i,j).v * (*Ubase)(i,j).v);
    }
  }
  ierr = thk->end_access(); CHKERRQ(ierr);
  ierr = Ubase->end_access(); CHKERRQ(ierr);
  ierr = result_Po.end_access(); CHKERRQ(ierr);
  ierr = result_cbase.end_access(); CHKERRQ(ierr);
  return 0;
}


//! Compute dt = min[t_end-t_current, dtCFL, dtDIFFW, dtDIFFP], the adaptive time step.
PetscErrorCode PISMHydrology::adaptive_time_step(PetscReal t_current, PetscReal t_end, 
                                                 PetscReal &dt_result) {
  PetscErrorCode ierr;
  PetscReal dtmax, dtCFL, dtDIFFW, dtDIFFP, maxW, maxH;

SETERRQ(grid.com,1,"not implemented");

  // FIXME: need series of global maxes
  // FIXME?: dtCFL can be infinity if velocity is zero because P and b are constant
  //dtCFL = 0.5 / (max(max(abs(alphV)))/dx + max(max(abs(betaV)))/dy);

  //maxW = max(max(max(Wea)),max(max(Wno))) + 0.001;
  //dtDIFFW = 0.25 / (p.K * maxW * (1/dx^2 + 1/dy^2));

  //maxH = max(max(h-b)) + 2 * p.E0;  % regularized: forces dtDIFFP < dtDIFFW
  dtDIFFP = (fresh_water_density * E0 / (ice_density * maxH)) * dtDIFFW;

  dtmax = (1.0/12.0) * secpera;  // FIXME: need better dtmax

  // dt = min([te-t dtmax dtCFL dtDIFFW dtDIFFP]);
  dt_result = PetscMin(t_end - t_current, dtmax);
  dt_result = PetscMin(dt_result, dtCFL);
  dt_result = PetscMin(dt_result, dtDIFFW);
  dt_result = PetscMin(dt_result, dtDIFFP);

  return 0;
}


//! Update the model state variables W,P by running the subglacial hydrology model.
/*!
Runs the hydrology model from time icet to time icet + icedt.  Here [icet,icedt]
is generally on the order of months to years.  This hydrology model will take its
own shorter time steps, perhaps hours to weeks.
 */
PetscErrorCode PISMHydrology::update(PetscReal icet, PetscReal icedt) {
  PetscErrorCode ierr;

  // if asked for the identical time interval versus last time, then
  //   do nothing; otherwise assume that [my_t,my_t+my_dt] is the time
  //   interval on which we are solving
  if ((fabs(icet - t) < 1e-6) && (fabs(icedt - dt) < 1e-6))
    return 0;
  // update PISMComponent times: t = current time, t+dt = target time
  t = icet;
  dt = icedt;

  // make sure W,P have valid ghosts before starting hydrology steps
  ierr = W.beginGhostComm(); CHKERRQ(ierr);
  ierr = P.beginGhostComm(); CHKERRQ(ierr);
  ierr = W.endGhostComm(); CHKERRQ(ierr);
  ierr = P.endGhostComm(); CHKERRQ(ierr);

  // from current ice geometry/velocity variables, initialize Po and cbase
  ierr = update_ice_functions(Po,cbase); CHKERRQ(ierr);

  PetscReal ht, hdt; // hydrology model time and time step
  while (ht < t + dt) {
    ierr = check_bounds(); CHKERRQ(ierr);

    ierr = velocity_staggered(V); CHKERRQ(ierr);
    ierr = water_thickness_staggered(Wstag); CHKERRQ(ierr);

    ierr = adaptive_time_step(ht, t+dt, hdt); CHKERRQ(ierr);

    ierr = hydraulic_potential(psi); CHKERRQ(ierr);

SETERRQ(grid.com,1,"not implemented");
    // FIXME ierr = known_state_mask(known); CHKERRQ(ierr);

    // update Pnew from time step
    PetscReal  pux = c0 / (grid.dx * grid.dx),
               puy = c0 / (grid.dy * grid.dy),
               Open, Close, divflux, Ptmp;
    ierr = P.begin_access(); CHKERRQ(ierr);
    ierr = W.begin_access(); CHKERRQ(ierr);
    ierr = cbase.begin_access(); CHKERRQ(ierr);
    ierr = psi.begin_access(); CHKERRQ(ierr);
    ierr = Wstag.begin_access(); CHKERRQ(ierr);
    ierr = bmelt->begin_access(); CHKERRQ(ierr);
    ierr = known.begin_access(); CHKERRQ(ierr);
    ierr = Pnew.begin_access(); CHKERRQ(ierr);
    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
        // opening and closure terms in pressure equation
        Open = PetscMax(0.0,c1 * cbase(i,j) * (Wr - W(i,j)));
        Close = c2 * Aglen * pow(Po(i,j) - P(i,j),nglen) * (W(i,j) + Y0);
        // divergence of flux
        divflux = 0;
        if (!known.as_int(i+1,j) && !known.as_int(i-1,j))
          divflux += pux * ( Wstag(i,j,0) * (psi(i+1,j) - psi(i,j))
                         - Wstag(i-1,j,0) * (psi(i,j) - psi(i-1,j)) );
        if (!known.as_int(i,j+1) && !known.as_int(i,j-1))
          divflux += puy * ( Wstag(i,j,1) * (psi(i,j+1) - psi(i,j))
                         - Wstag(i,j-1,1) * (psi(i,j) - psi(i,j-1)) );
        // candidate for update
        Ptmp = P(i,j) + (hdt * Po(i,j) / E0) * ( divflux + Close - Open + (*bmelt)(i,j) );
        // projection:
        Pnew(i,j) = PetscMin(PetscMax(0.0, Ptmp), Po(i,j));
      }
    }
    ierr = P.end_access(); CHKERRQ(ierr);
    ierr = W.end_access(); CHKERRQ(ierr);
    ierr = cbase.end_access(); CHKERRQ(ierr);
    ierr = Pnew.end_access(); CHKERRQ(ierr);
    ierr = psi.end_access(); CHKERRQ(ierr);
    ierr = bmelt->end_access(); CHKERRQ(ierr);
    ierr = Wstag.end_access(); CHKERRQ(ierr);
    ierr = known.end_access(); CHKERRQ(ierr);

SETERRQ(grid.com,1,"not implemented");
    // FIXME: start transfer Pnew into P

    ierr = advective_fluxes(Qstag); CHKERRQ(ierr);

    // update Wnew from time step
    PetscReal  wux = K / (grid.dx * grid.dx),
               wuy = K / (grid.dy * grid.dy),
               divadflux, diffW;
    ierr = W.begin_access(); CHKERRQ(ierr);
    ierr = Wstag.begin_access(); CHKERRQ(ierr);
    ierr = Qstag.begin_access(); CHKERRQ(ierr);
    ierr = bmelt->begin_access(); CHKERRQ(ierr);
    ierr = Wnew.begin_access(); CHKERRQ(ierr);
    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
        divadflux =   (Qstag(i,j,0) - Qstag(i-1,j  ,0)) / grid.dx
                    + (Qstag(i,j,1) - Qstag(i,  j-1,1)) / grid.dy;
        diffW =   wux * (Wstag(i,j,0) * (W(i+1,j  ) - W(i,j)) - Wstag(i-1,j  ,0) * (W(i,j) - W(i-1,  j)))
                + wuy * (Wstag(i,j,1) * (W(i  ,j+1) - W(i,j)) - Wstag(i  ,j-1,1) * (W(i,j) - W(i  ,j-1)));
        Wnew(i,j) = W(i,j) + hdt * (- divadflux + diffW + (*bmelt)(i,j));
      }
    }
    ierr = W.end_access(); CHKERRQ(ierr);
    ierr = Wstag.end_access(); CHKERRQ(ierr);
    ierr = Qstag.end_access(); CHKERRQ(ierr);
    ierr = bmelt->end_access(); CHKERRQ(ierr);
    ierr = Wnew.end_access(); CHKERRQ(ierr);

SETERRQ(grid.com,1,"not implemented");
    // FIXME: start transfer Wnew into W

    ht += hdt;
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

