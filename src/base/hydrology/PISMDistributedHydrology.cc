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


/************************************/
/******** PISMLakesHydrology ********/
/************************************/

PISMLakesHydrology::PISMLakesHydrology(IceGrid &g, const NCConfigVariable &conf)
    : PISMHydrology(g, conf)
{
    bed   = NULL;
    thk   = NULL;
    usurf = NULL;
    mask  = NULL;
    bmelt = NULL;

    if (allocate() != 0) {
      PetscPrintf(grid.com, "PISM ERROR: memory allocation failed in PISMLakesHydrology constructor.\n");
      PISMEnd();
    }

    ice_density = config.get("ice_density");
    standard_gravity = config.get("standard_gravity");
    fresh_water_density = config.get("fresh_water_density");
    sea_water_density = config.get("sea_water_density");
}


PetscErrorCode PISMLakesHydrology::allocate() {
  PetscErrorCode ierr;

  // model state variables; need ghosts
  ierr = W.create(grid, "bwat", true, 1); CHKERRQ(ierr);
  ierr = W.set_attrs("model_state",
                     "thickness of subglacial water layer",
                     "m", ""); CHKERRQ(ierr);
  ierr = W.set_attr("valid_min", 0.0); CHKERRQ(ierr);

  // auxiliary variables which NEED ghosts
  ierr = psi.create(grid, "hydraulic_potential", true, 1); CHKERRQ(ierr);
  ierr = psi.set_attrs("internal",
                       "hydraulic potential of water in subglacial layer",
                       "Pa", ""); CHKERRQ(ierr);
  ierr = Wstag.create(grid, "W_staggered", true, 1); CHKERRQ(ierr);
  ierr = Wstag.set_attrs("internal",
                     "cell face-centered (staggered) values of water layer thickness",
                     "m", ""); CHKERRQ(ierr);
  ierr = Qstag.create(grid, "advection_flux", true, 1); CHKERRQ(ierr);
  ierr = Qstag.set_attrs("internal",
                     "cell face-centered (staggered) components of advective subglacial water flux",
                     "m2 s-1", ""); CHKERRQ(ierr);
  ierr = Pwork.create(grid, "water_pressure_workspace", true, t); CHKERRQ(ierr);
  ierr = Pwork.set_attrs("internal",
                      "work space for modeled subglacial water pressure",
                      "Pa", ""); CHKERRQ(ierr);
  ierr = Pwork.set_attr("valid_min", 0.0); CHKERRQ(ierr);

  // auxiliary variables which do not need ghosts
  ierr = V.create(grid, "water_velocity", false); CHKERRQ(ierr);
  ierr = V.set_attrs("internal",
                     "cell face-centered (staggered) components of water velocity in subglacial water layer",
                     "m s-1", ""); CHKERRQ(ierr);

  // temporaries during update; do not need ghosts
  //FIXME: ought to be able to do this (issue #109):
  //  ierr = Wnew.create(grid, "Wnew_internal", false); CHKERRQ(ierr);
  ierr = Wnew.create(grid, "Wnew_internal", true, 1); CHKERRQ(ierr);
  ierr = Wnew.set_attrs("internal",
                     "new thickness of subglacial water layer during update",
                     "m", ""); CHKERRQ(ierr);
  ierr = Wnew.set_attr("valid_min", 0.0); CHKERRQ(ierr);

  return 0;
}


PetscErrorCode PISMLakesHydrology::init(PISMVars &vars) {
  PetscErrorCode ierr;

  ierr = verbPrintf(2, grid.com,
    "* Initializing the subglacial-lakes-suitable subglacial hydrology model...\n"); CHKERRQ(ierr);

  bed = dynamic_cast<IceModelVec2S*>(vars.get("topg"));
  if (bed == NULL) SETERRQ(grid.com, 1, "topg is not available");

  thk = dynamic_cast<IceModelVec2S*>(vars.get("thk"));
  if (thk == NULL) SETERRQ(grid.com, 1, "thk is not available");

  usurf = dynamic_cast<IceModelVec2S*>(vars.get("usurf"));
  if (usurf == NULL) SETERRQ(grid.com, 1, "usurf is not available");

  mask = dynamic_cast<IceModelVec2Int*>(vars.get("mask"));
  if (mask == NULL) SETERRQ(grid.com, 1, "mask is not available");

  bmelt = dynamic_cast<IceModelVec2S*>(vars.get("bmelt"));
  if (bmelt == NULL) SETERRQ(grid.com, 1, "bmelt is not available");

  // initialize water layer thickness from the context if present, otherwise zero
  IceModelVec2S *W_input = dynamic_cast<IceModelVec2S*>(vars.get("bwat"));
  if (W_input != NULL) {
    ierr = W.copy_from(*W_input); CHKERRQ(ierr);
  } else {
    ierr = W.set(0.0); CHKERRQ(ierr);
  }

  if (vars.get("bwat") == NULL) {
    ierr = vars.add(W); CHKERRQ(ierr);
  }
  return 0;
}


void PISMLakesHydrology::add_vars_to_output(string /*keyword*/, map<string,NCSpatialVariable> &result) {
  result["bwat"] = W.get_metadata();
  IceModelVec2S tmp;
  tmp.create(grid, "bwp", false);
  tmp.set_attrs("diagnostic","pressure of water in subglacial layer",
                "Pa", "");
  tmp.set_attr("valid_min", 0.0);
  result["bwp"] = tmp.get_metadata();
  // destructor called on tmp when we go out of scope here
}


PetscErrorCode PISMLakesHydrology::define_variables(set<string> vars, const PIO &nc,
                                                 PISM_IO_Type nctype) {
  PetscErrorCode ierr;
  if (set_contains(vars, "bwat")) {
    ierr = W.define(nc, nctype); CHKERRQ(ierr);
  }
  if (set_contains(vars, "bwp")) {
    IceModelVec2S tmp;
    ierr = tmp.create(grid, "bwp", false); CHKERRQ(ierr);
    ierr = tmp.set_attrs("diagnostic",
                     "pressure of water in subglacial layer",
                     "Pa", ""); CHKERRQ(ierr);
    ierr = tmp.set_attr("valid_min", 0.0); CHKERRQ(ierr);
    ierr = tmp.define(nc, nctype); CHKERRQ(ierr);
    // destructor called on tmp when we go out of scope here
  }
  return 0;
}


PetscErrorCode PISMLakesHydrology::write_variables(set<string> vars, const PIO &nc) {
  PetscErrorCode ierr;
  if (set_contains(vars, "bwat")) {
    ierr = W.write(nc); CHKERRQ(ierr);
  }
  if (set_contains(vars, "bwp")) {
    IceModelVec2S tmp;
    ierr = tmp.create(grid, "bwp", false); CHKERRQ(ierr);
    ierr = tmp.set_attrs("diagnostic",
                     "pressure of water in subglacial layer",
                     "Pa", ""); CHKERRQ(ierr);
    ierr = tmp.set_attr("valid_min", 0.0); CHKERRQ(ierr);
    ierr = water_pressure(tmp); CHKERRQ(ierr);
    ierr = tmp.write(nc); CHKERRQ(ierr);
    // destructor called on tmp when we go out of scope here
  }
  return 0;
}


//! Check W >= 0 and fails with message if not satisfied.
PetscErrorCode PISMLakesHydrology::check_Wpositive() {
  PetscErrorCode ierr;
  ierr = W.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (W(i,j) < 0.0) {
        PetscPrintf(grid.com,
           "PISM ERROR: disallowed negative subglacial water layer thickness\n"
           "    W(i,j) = %.6f m at (i,j)=(%d,%d)\n"
           "ENDING ... \n\n", W(i,j),i,j);
        PISMEnd();
      }
    }
  }
  ierr = W.end_access(); CHKERRQ(ierr);
  return 0;
}


//! Copies the W variable, the modeled water layer thickness.
PetscErrorCode PISMLakesHydrology::water_layer_thickness(IceModelVec2S &result) {
  PetscErrorCode ierr = W.copy_to(result); CHKERRQ(ierr);
  return 0;
}


//! Computes pressure diagnostically as fixed fraction of overburden.
/*!
Here
  \f[ P = \lambda P_o = \lambda (\rho_i g H) \f]
where \f$\lambda\f$=till_pw_fraction and \f$P_o\f$ is the overburden pressure.
 */
PetscErrorCode PISMLakesHydrology::water_pressure(IceModelVec2S &result) {
  PetscErrorCode ierr;
  ierr = update_overburden(result); CHKERRQ(ierr);
  ierr = result.scale(config.get("till_pw_fraction")); CHKERRQ(ierr);
  return 0;
}


//! Update the overburden pressure from ice thickness.
/*!
Accesses thk from PISMVars, which points into IceModel.
 */
PetscErrorCode PISMLakesHydrology::update_overburden(IceModelVec2S &result) {
  PetscErrorCode ierr;
  // P_o = rho_i g H
  ierr = result.copy_from(*thk); CHKERRQ(ierr);
  ierr = result.scale(ice_density * standard_gravity); CHKERRQ(ierr);
  return 0;
}


//! Get the hydraulic potential from bedrock topography and current state variables.
/*!
Computes \f$\psi = P + \rho_w g (b + W)\f$.

Calls water_pressure() method to get water pressure.
 */
PetscErrorCode PISMLakesHydrology::hydraulic_potential(IceModelVec2S &result) {
  PetscErrorCode ierr;
  ierr = water_pressure(Pwork); CHKERRQ(ierr);
  ierr = W.begin_access(); CHKERRQ(ierr);
  ierr = Pwork.begin_access(); CHKERRQ(ierr);
  ierr = bed->begin_access(); CHKERRQ(ierr);
  ierr = result.begin_access(); CHKERRQ(ierr);
  for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
      result(i,j) = Pwork(i,j)
                    + fresh_water_density * standard_gravity * ((*bed)(i,j) + W(i,j));
    }
  }
  ierr = W.end_access(); CHKERRQ(ierr);
  ierr = Pwork.end_access(); CHKERRQ(ierr);
  ierr = bed->end_access(); CHKERRQ(ierr);
  ierr = result.end_access(); CHKERRQ(ierr);
  return 0;
}


//! Get the advection velocity V at the center of cell edges.
/*!
Computes the advection velocity \f$\mathbf{V}=\mathbf{V}(\nabla P,\nabla b)\f$
on the staggered (face-centered) grid.  If V = (alpha,beta) in components
then we have <code> result(i,j,0) = alpha(i+1/2,j) </code> and
<code> result(i,j,1) = beta(i,j+1/2) </code>

The advection velocity is given by the formula
  \f[ \mathbf{V} = - \frac{K}{\rho_w g} \nabla P - K \nabla b \f]
where \f$\mathbf{V}\f$ is the lateral water velocity, \f$P\f$ is the water
pressure, and \f$b\f$ is the bedrock elevation.

Calls water_pressure() method to get water pressure.
 */
PetscErrorCode PISMLakesHydrology::velocity_staggered(IceModelVec2Stag &result) {
  PetscErrorCode ierr;
  PetscReal dbdx, dbdy, dPdx, dPdy, K, c0;
  // FIXME:  want Kmax or Kmin according to W > Wr ?
  K  = config.get("hydrology_hydraulic_conductivity");
  c0 = K / (fresh_water_density * standard_gravity);
  ierr = water_pressure(Pwork); CHKERRQ(ierr);  // does update ghosts
  ierr = Pwork.begin_access(); CHKERRQ(ierr);
  ierr = bed->begin_access(); CHKERRQ(ierr);
  ierr = result.begin_access(); CHKERRQ(ierr);
  for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
      dPdx = (Pwork(i+1,j  ) - Pwork(i,j)) / grid.dx;
      dPdy = (Pwork(i  ,j+1) - Pwork(i,j)) / grid.dy;
      dbdx = ((*bed)(i+1,j  ) - (*bed)(i,j)) / grid.dx;
      dbdy = ((*bed)(i  ,j+1) - (*bed)(i,j)) / grid.dy;
      result(i,j,0) = - c0 * dPdx - K * dbdx;
      result(i,j,1) = - c0 * dPdy - K * dbdy;
    }
  }
  ierr = Pwork.end_access(); CHKERRQ(ierr);
  ierr = bed->end_access(); CHKERRQ(ierr);
  ierr = result.end_access(); CHKERRQ(ierr);
  return 0;
}


//! Average the regular grid water thickness to values at the center of cell edges.
PetscErrorCode PISMLakesHydrology::water_thickness_staggered(IceModelVec2Stag &result) {
  PetscErrorCode ierr;
  ierr = W.begin_access(); CHKERRQ(ierr);
  ierr = result.begin_access(); CHKERRQ(ierr);
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
/*!
The field W must have valid ghost values, but V does not need them.
 */
PetscErrorCode PISMLakesHydrology::advective_fluxes(IceModelVec2Stag &result) {
  PetscErrorCode ierr;
  ierr = W.begin_access(); CHKERRQ(ierr);
  ierr = V.begin_access(); CHKERRQ(ierr);
  ierr = result.begin_access(); CHKERRQ(ierr);
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


//! Compute the adaptive time step for evolution of W.
PetscErrorCode PISMLakesHydrology::adaptive_for_W_evolution(
                  PetscReal t_current, PetscReal t_end, PetscReal &dt_result,
                  PetscReal &dt_DIFFW_result) {
  PetscErrorCode ierr;
  PetscReal dtmax, dtCFL, maxW, K;
  PetscReal tmp[2];
  // FIXME:  want Kmax or Kmin according to W > Wr ?
  K  = config.get("hydrology_hydraulic_conductivity");
  // fixme?: dtCFL can be infinity if velocity is zero because P and b are constant
  // Matlab: dtCFL = 0.5 / (max(max(abs(alphV)))/dx + max(max(abs(betaV)))/dy);
  ierr = V.absmaxcomponents(tmp); CHKERRQ(ierr);
  dtCFL = 0.5 / (tmp[0]/grid.dx + tmp[1]/grid.dy);
  // Matlab: maxW = max(max(max(Wea)),max(max(Wno))) + 0.001;
  ierr = Wstag.absmaxcomponents(tmp); CHKERRQ(ierr);
  maxW = PetscMax(tmp[0],tmp[1]) + 0.001;
  // Matlab: dtDIFFW = 0.25 / (p.K * maxW * (1/dx^2 + 1/dy^2));
  dt_DIFFW_result = 1.0/(grid.dx*grid.dx) + 1.0/(grid.dy*grid.dy);
  dt_DIFFW_result = 0.25 / (K * maxW * dt_DIFFW_result);
  dtmax = (1.0/12.0) * secpera;  // fixme?: need better dtmax than fixed at 1 month
  // dt = min([te-t dtmax dtCFL dtDIFFW]);
  dt_result = PetscMin(t_end - t_current, dtmax);
  dt_result = PetscMin(dt_result, dtCFL);
  dt_result = PetscMin(dt_result, dt_DIFFW_result);
  return 0;
}


/*!
Normally call this version if you don't need the dt associated to the diffusion term.
 */
PetscErrorCode PISMLakesHydrology::adaptive_for_W_evolution(
                  PetscReal t_current, PetscReal t_end, PetscReal &dt_result) {
  PetscReal discard;
  PetscErrorCode ierr = adaptive_for_W_evolution(
                             t_current, t_end, dt_result, discard); CHKERRQ(ierr);
  return 0;
}


//! Update the model state variable W by running the subglacial hydrology model.
/*!
Runs the hydrology model from time icet to time icet + icedt.  Here [icet,icedt]
is generally on the order of months to years.  This hydrology model will take its
own shorter time steps, perhaps hours to weeks.
 */
PetscErrorCode PISMLakesHydrology::update(PetscReal icet, PetscReal icedt) {
  PetscErrorCode ierr;

  // if asked for the identical time interval versus last time, then
  //   do nothing; otherwise assume that [my_t,my_t+my_dt] is the time
  //   interval on which we are solving
  if ((fabs(icet - t) < 1e-12) && (fabs(icedt - dt) < 1e-12))
    return 0;
  // update PISMComponent times: t = current time, t+dt = target time
  t = icet;
  dt = icedt;

  // make sure W has valid ghosts before starting hydrology steps
  ierr = W.beginGhostComm(); CHKERRQ(ierr);
  ierr = W.endGhostComm(); CHKERRQ(ierr);

  MaskQuery M(*mask);

  PetscReal ht = t, hdt, // hydrology model time and time step
            K = config.get("hydrology_hydraulic_conductivity"), // FIXME: want Kmax or Kmin according to W > Wr ?
            my_icefreelost = 0.0, my_oceanlost = 0.0, my_negativegain = 0.0,
            dA = grid.dx * grid.dy;
  bool      use_const   = config.get_flag("hydrology_use_const_bmelt");
  PetscReal const_bmelt = config.get("hydrology_const_bmelt");
  PetscInt hydrocount = 0; // count hydrology time steps

//ierr = PetscPrintf(grid.com, "starting PISMLakesHydrology::update() time-stepping loop\n"); CHKERRQ(ierr);
  while (ht < t + dt) {
    hydrocount++;
    ierr = check_Wpositive(); CHKERRQ(ierr);

    ierr = hydraulic_potential(psi); CHKERRQ(ierr);
    ierr = psi.beginGhostComm(); CHKERRQ(ierr);
    ierr = psi.endGhostComm(); CHKERRQ(ierr);

    ierr = velocity_staggered(V); CHKERRQ(ierr);

    ierr = water_thickness_staggered(Wstag); CHKERRQ(ierr);
    ierr = Wstag.beginGhostComm(); CHKERRQ(ierr);
    ierr = Wstag.endGhostComm(); CHKERRQ(ierr);

    // to get Qstag, W needs valid ghosts
    ierr = advective_fluxes(Qstag); CHKERRQ(ierr);
    ierr = Qstag.beginGhostComm(); CHKERRQ(ierr);
    ierr = Qstag.endGhostComm(); CHKERRQ(ierr);

    ierr = adaptive_for_W_evolution(ht, t+dt, hdt); CHKERRQ(ierr);
//ierr = PetscPrintf(grid.com, "adaptive...() reports hdt = %.6f seconds = %.6f years\n",hdt,hdt/secpera); CHKERRQ(ierr);

    // update Wnew from time step
    PetscReal  wux = K / (grid.dx * grid.dx),  // FIXME if K is variable
               wuy = K / (grid.dy * grid.dy),
               divadflux, diffW;
    ierr = W.begin_access(); CHKERRQ(ierr);
    ierr = mask->begin_access(); CHKERRQ(ierr);
    ierr = Wstag.begin_access(); CHKERRQ(ierr);
    ierr = Qstag.begin_access(); CHKERRQ(ierr);
    ierr = bmelt->begin_access(); CHKERRQ(ierr);
    ierr = Wnew.begin_access(); CHKERRQ(ierr);
    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
        divadflux =   (Qstag(i,j,0) - Qstag(i-1,j  ,0)) / grid.dx
                    + (Qstag(i,j,1) - Qstag(i,  j-1,1)) / grid.dy;
        diffW =   wux * (  Wstag(i,j,0) * (W(i+1,j  ) - W(i,j))
                         - Wstag(i-1,j  ,0) * (W(i,j) - W(i-1,  j)) )
                + wuy * ( Wstag(i,j,1) * (W(i  ,j+1) - W(i,j))
                         - Wstag(i  ,j-1,1) * (W(i,j) - W(i  ,j-1)) );
        PetscReal input = (use_const) ? const_bmelt : (*bmelt)(i,j);
        Wnew(i,j) = W(i,j) + hdt * (- divadflux + diffW + input);
        if (M.ice_free_land(i,j)) {
          my_icefreelost += Wnew(i,j) * dA * fresh_water_density; // FIXME: mult by cell area?
          Wnew(i,j) = 0.0;
        } else if (M.ocean(i,j)) {
          my_oceanlost += Wnew(i,j) * dA * fresh_water_density; // FIXME: mult by cell area?
          Wnew(i,j) = 0.0;
        } else if (Wnew(i,j) < 0.0) {
          my_negativegain += -Wnew(i,j) * dA * fresh_water_density; // FIXME: mult by cell area?
          Wnew(i,j) = 0.0;
        }
      }
    }
    ierr = W.end_access(); CHKERRQ(ierr);
    ierr = mask->end_access(); CHKERRQ(ierr);
    ierr = Wstag.end_access(); CHKERRQ(ierr);
    ierr = Qstag.end_access(); CHKERRQ(ierr);
    ierr = bmelt->end_access(); CHKERRQ(ierr);
    ierr = Wnew.end_access(); CHKERRQ(ierr);

    // transfer Wnew into W
    ierr = Wnew.beginGhostComm(W); CHKERRQ(ierr);
    ierr = Wnew.endGhostComm(W); CHKERRQ(ierr);

    ht += hdt;
  } // end of hydrology model time-stepping loop

  // lost mass over whole glacier/ice sheet
  PetscReal icefreelost, oceanlost, negativegain;
  ierr = PISMGlobalSum(&my_icefreelost, &icefreelost, grid.com); CHKERRQ(ierr);
  ierr = PISMGlobalSum(&my_oceanlost, &oceanlost, grid.com); CHKERRQ(ierr);
  ierr = PISMGlobalSum(&my_negativegain, &negativegain, grid.com); CHKERRQ(ierr);

  ierr = verbPrintf(2, grid.com,
    " 'lakes' distributed hydrology:  %d steps with average dt = %.6f years;  mass losses:\n"
    "     ice free land lost = %.3e kg, ocean lost = %.3e kg, negative bmelt gain = %.3e kg\n",
    hydrocount, (dt/hydrocount)/secpera, icefreelost, oceanlost, negativegain); CHKERRQ(ierr);
  return 0;
}



/************************************/
/***** PISMDistributedHydrology *****/
/************************************/

PISMDistributedHydrology::PISMDistributedHydrology(IceGrid &g, const NCConfigVariable &conf,
                                                   PISMStressBalance *sb)
    : PISMLakesHydrology(g, conf)
{
    stressbalance = sb;

    if (allocatePstuff() != 0) {
      PetscPrintf(grid.com, "PISM ERROR: memory allocation failed in PISMDistributedHydrology constructor.\n");
      PISMEnd();
    }

    ice_density = config.get("ice_density");
    standard_gravity = config.get("standard_gravity");
    fresh_water_density = config.get("fresh_water_density");
    sea_water_density = config.get("sea_water_density");

    // initialize using constants from van Pelt & Bueler preprint
    // FIXME: should be configurable
    c1    = 0.500;      // m-1
    c2    = 0.040;      // [pure]
    Aglen = 3.1689e-24; // Pa-3 s-1; ice softness
    nglen = 3.0;
    Wr    = 1.0;        // m
    E0    = 1.0;        // m; what is optimal?
    Y0    = 0.001;      // m; regularization
}


PetscErrorCode PISMDistributedHydrology::allocatePstuff() {
  PetscErrorCode ierr;

  // additional variables beyond PISMLakesHydrology::allocate()
  ierr = P.create(grid, "bwpdistributed", true, 1); CHKERRQ(ierr);
  ierr = P.set_attrs("model_state",
                     "pressure of water in subglacial layer",
                     "Pa", ""); CHKERRQ(ierr);
  ierr = P.set_attr("valid_min", 0.0); CHKERRQ(ierr);
  ierr = cbase.create(grid, "ice_sliding_speed", false); CHKERRQ(ierr);
  ierr = cbase.set_attrs("internal",
                         "ice sliding speed seen by subglacial water layer",
                         "m s-1", ""); CHKERRQ(ierr);
  ierr = cbase.set_attr("valid_min", 0.0); CHKERRQ(ierr);
  ierr = Pnew.create(grid, "Pnew_internal", false); CHKERRQ(ierr);
  ierr = Pnew.set_attrs("internal",
                     "new subglacial water pressure during update",
                     "Pa", ""); CHKERRQ(ierr);
  ierr = Pnew.set_attr("valid_min", 0.0); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode PISMDistributedHydrology::init(PISMVars &vars) {
  PetscErrorCode ierr;

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

  if (vars.get("bwat") == NULL) {
    ierr = vars.add(W); CHKERRQ(ierr);
  }
  if (vars.get("bwp") == NULL) {
    ierr = vars.add(P); CHKERRQ(ierr);
  }
  return 0;
}


void PISMDistributedHydrology::add_vars_to_output(string /*keyword*/, map<string,NCSpatialVariable> &result) {
  result["bwatdistributed"] = W.get_metadata();
  result["bwpdistributed"]  = P.get_metadata();
}


PetscErrorCode PISMDistributedHydrology::define_variables(set<string> vars, const PIO &nc,
                                                 PISM_IO_Type nctype) {
  PetscErrorCode ierr;
  if (set_contains(vars, "bwatdistributed")) {
    ierr = W.define(nc, nctype); CHKERRQ(ierr);
  }
  if (set_contains(vars, "bwpdistributed")) {
    ierr = P.define(nc, nctype); CHKERRQ(ierr);
  }
  return 0;
}


PetscErrorCode PISMDistributedHydrology::write_variables(set<string> vars, const PIO &nc) {
  PetscErrorCode ierr;
  if (set_contains(vars, "bwatdistributed")) {
    ierr = W.write(nc); CHKERRQ(ierr);
  }
  if (set_contains(vars, "bwpdistributed")) {
    ierr = P.write(nc); CHKERRQ(ierr);
  }
  return 0;
}


//! Copies the P state variable which is the modeled water pressure.
PetscErrorCode PISMDistributedHydrology::water_pressure(IceModelVec2S &result) {
  PetscErrorCode ierr = P.copy_to(result); CHKERRQ(ierr);
  return 0;
}


//! Check bounds on W and P and fail with message if not satisfied.
/*!
Checks \f$0 \le W\f$ and \f$0 \le P \le P_o\f$.
 */
PetscErrorCode PISMDistributedHydrology::check_bounds() {
  PetscErrorCode ierr;
  ierr = check_Wpositive(); CHKERRQ(ierr);
  ierr = P.begin_access(); CHKERRQ(ierr);
  ierr = Po.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (P(i,j) < 0.0) {
        PetscPrintf(grid.com,
           "PISM ERROR: disallowed negative subglacial water pressure\n"
           "    P(i,j) = %.6f Pa\n at (i,j)=(%d,%d)\n"
           "ENDING ... \n\n", P(i,j),i,j);
        PISMEnd();
      }
      if (P(i,j) > Po(i,j)) {
        PetscPrintf(grid.com,
           "PISM ERROR: subglacial water pressure P(i,j) = %.6f Pa exceeds\n"
           "    overburden pressure Po(i,j) = %.6f Pa at (i,j)=(%d,%d)\n"
           "ENDING ... \n\n", P(i,j),Po(i,j),i,j);
        PISMEnd();
      }
    }
  }
  ierr = P.end_access(); CHKERRQ(ierr);
  ierr = Po.end_access(); CHKERRQ(ierr);
  return 0;
}


//! Compute functional relationship P(W) which applies only in steady state.
/*!
In steady state in this model, water pressure is determined by a balance of
cavitation (opening) caused by sliding and creep closure.

This will be used in initialization when P is otherwise unknown, and
in verification and/or reporting.  It is not used during time-dependent
model runs.  To be more complete, \f$P=P(W,P_o,|v_b|)\f$.
 */
PetscErrorCode PISMDistributedHydrology::P_from_W_steady(IceModelVec2S &result) {
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


//! Update the the sliding speed |v_b| from ice quantities.
/*!
Calls a PISMStressBalance method to get the vector basal velocity of the ice,
and then computes the magnitude of that.

fixme:  Is taking Ubase from the PISMStressBalance the correct method?
 */
PetscErrorCode PISMDistributedHydrology::update_cbase(IceModelVec2S &result_cbase) {
  PetscErrorCode ierr;
  IceModelVec2V* Ubase; // ice sliding velocity
  // cbase = |v_b|
  ierr = stressbalance->get_2D_advective_velocity(Ubase); CHKERRQ(ierr);
  ierr = Ubase->magnitude(result_cbase); CHKERRQ(ierr);
  return 0;
}


//! Computes the adaptive time step for this (W,P) state space model.
PetscErrorCode PISMDistributedHydrology::adaptive_for_WandP_evolution(
                  PetscReal t_current, PetscReal t_end, PetscReal &dt_result) {
  PetscErrorCode ierr;
  PetscReal dtDIFFW, dtDIFFP, maxH;
  ierr = adaptive_for_W_evolution(t_current,t_end,dt_result,dtDIFFW); CHKERRQ(ierr);

  // Matlab: dtDIFFP = (p.rhow * p.E0 / (p.rhoi * maxH)) * dtDIFFW;
  ierr = thk->max(maxH); CHKERRQ(ierr);
  maxH += 2.0 * E0; // regularized: forces dtDIFFP < dtDIFFW
  dtDIFFP = (fresh_water_density * E0 / (ice_density * maxH)) * dtDIFFW;

  // dt = min([te-t dtmax dtCFL dtDIFFW dtDIFFP]);
  dt_result = PetscMin(dt_result, dtDIFFP);
  return 0;
}


//! Update the model state variables W,P by running the subglacial hydrology model.
/*!
Runs the hydrology model from time icet to time icet + icedt.  Here [icet,icedt]
is generally on the order of months to years.  This hydrology model will take its
own shorter time steps, perhaps hours to weeks.
 */
PetscErrorCode PISMDistributedHydrology::update(PetscReal icet, PetscReal icedt) {
  PetscErrorCode ierr;

  // if asked for the identical time interval versus last time, then
  //   do nothing; otherwise assume that [my_t,my_t+my_dt] is the time
  //   interval on which we are solving
  if ((fabs(icet - t) < 1e-12) && (fabs(icedt - dt) < 1e-12))
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
  ierr = update_overburden(Po); CHKERRQ(ierr);
  ierr = update_cbase(cbase); CHKERRQ(ierr);

  PetscReal ht, hdt, // hydrology model time and time step
            K, c0;
  // FIXME:  want Kmax or Kmin according to W > Wr ?
  K  = config.get("hydrology_hydraulic_conductivity");
  c0 = K / (fresh_water_density * standard_gravity);

  bool      use_const   = config.get_flag("hydrology_use_const_bmelt");
  PetscReal const_bmelt = config.get("hydrology_const_bmelt");

  while (ht < t + dt) {
    ierr = check_bounds(); CHKERRQ(ierr);

    ierr = hydraulic_potential(psi); CHKERRQ(ierr);
    ierr = psi.beginGhostComm(); CHKERRQ(ierr);

    ierr = velocity_staggered(V); CHKERRQ(ierr);

    ierr = water_thickness_staggered(Wstag); CHKERRQ(ierr);
    ierr = Wstag.beginGhostComm(); CHKERRQ(ierr);

    // to get Qstag, W needs valid ghosts
    ierr = advective_fluxes(Qstag); CHKERRQ(ierr);
    ierr = Qstag.beginGhostComm(); CHKERRQ(ierr);

    ierr = adaptive_for_WandP_evolution(ht, t+dt, hdt); CHKERRQ(ierr);

    ierr = psi.endGhostComm(); CHKERRQ(ierr);
    ierr = Wstag.endGhostComm(); CHKERRQ(ierr);

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
    ierr = mask->begin_access(); CHKERRQ(ierr);
    ierr = Pnew.begin_access(); CHKERRQ(ierr);
    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
        // opening and closure terms in pressure equation
        Open = PetscMax(0.0,c1 * cbase(i,j) * (Wr - W(i,j)));
        Close = c2 * Aglen * pow(Po(i,j) - P(i,j),nglen) * (W(i,j) + Y0);
        // divergence of flux
        divflux = 0;
        //FIXME if (!known.as_int(i+1,j) && !known.as_int(i-1,j))
          divflux += pux * ( Wstag(i,j,0) * (psi(i+1,j) - psi(i,j))
                         - Wstag(i-1,j,0) * (psi(i,j) - psi(i-1,j)) );
        //FIXME if (!known.as_int(i,j+1) && !known.as_int(i,j-1))
          divflux += puy * ( Wstag(i,j,1) * (psi(i,j+1) - psi(i,j))
                         - Wstag(i,j-1,1) * (psi(i,j) - psi(i,j-1)) );
        // candidate for update
        PetscReal input = (use_const) ? const_bmelt : (*bmelt)(i,j);
        Ptmp = P(i,j) + (hdt * Po(i,j) / E0) * ( divflux + Close - Open + input );
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
    ierr = mask->end_access(); CHKERRQ(ierr);

    // start transfer Pnew into P; note Wstag, Qstag unaffected in Wnew update below
    ierr = Pnew.beginGhostComm(P); CHKERRQ(ierr);

    ierr = Qstag.endGhostComm(); CHKERRQ(ierr);

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
        diffW =   wux * (  Wstag(i,j,0) * (W(i+1,j  ) - W(i,j))
                         - Wstag(i-1,j  ,0) * (W(i,j) - W(i-1,  j)) )
                + wuy * ( Wstag(i,j,1) * (W(i  ,j+1) - W(i,j))
                         - Wstag(i  ,j-1,1) * (W(i,j) - W(i  ,j-1)) );
        Wnew(i,j) = W(i,j) + hdt * (- divadflux + diffW + (*bmelt)(i,j));
      }
    }
    ierr = W.end_access(); CHKERRQ(ierr);
    ierr = Wstag.end_access(); CHKERRQ(ierr);
    ierr = Qstag.end_access(); CHKERRQ(ierr);
    ierr = bmelt->end_access(); CHKERRQ(ierr);
    ierr = Wnew.end_access(); CHKERRQ(ierr);

    // start transfer Wnew into W
    ierr = Wnew.beginGhostComm(W); CHKERRQ(ierr);

    // finalizes update of P, W
    ierr = Pnew.endGhostComm(P); CHKERRQ(ierr);
    ierr = Wnew.endGhostComm(W); CHKERRQ(ierr);

    ht += hdt;
  } // end of hydrology model time-stepping loop

  return 0;
}

