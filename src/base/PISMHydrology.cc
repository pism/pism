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


PISMTillCanHydrology::PISMTillCanHydrology(IceGrid &g, const NCConfigVariable &conf, bool Whasghosts)
    : PISMHydrology(g, conf)
{
    thk   = NULL;
    bmelt  = NULL;
    if (allocate(Whasghosts) != 0) {
      PetscPrintf(grid.com,
        "PISM ERROR: allocation failed in PISMTillCanHydrology constructor.\n");
      PISMEnd();
    }
}


PetscErrorCode PISMTillCanHydrology::allocate(bool Whasghosts) {
  PetscErrorCode ierr;
  // model state variables
  if (Whasghosts) {
    ierr = W.create(grid, "bwattillcan", true, 1); CHKERRQ(ierr);
  } else {
    ierr = W.create(grid, "bwattillcan", false); CHKERRQ(ierr);
  }
  ierr = W.set_attrs("model_state",
                     "thickness of subglacial water layer",
                     "m", ""); CHKERRQ(ierr);
  ierr = W.set_attr("valid_min", 0.0); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode PISMTillCanHydrology::init(PISMVars &vars) {
  PetscErrorCode ierr;
  ierr = verbPrintf(2, grid.com,
    "* Initializing the till-can subglacial hydrology model...\n"); CHKERRQ(ierr);
  variables = &vars;
  thk = dynamic_cast<IceModelVec2S*>(vars.get("thk"));
  if (thk == NULL) SETERRQ(grid.com, 1, "thk is not available");
  bmelt = dynamic_cast<IceModelVec2S*>(vars.get("bmelt"));
  if (bmelt == NULL) SETERRQ(grid.com, 1, "bmelt is not available");
  // initialize water layer thickness from the context if present, otherwise zero
  IceModelVec2S *W_input = dynamic_cast<IceModelVec2S*>(vars.get("bwat"));
  if (W_input != NULL) {
    ierr = W.copy_from(*W_input); CHKERRQ(ierr);
  } else {
    ierr = W.set(0.0); CHKERRQ(ierr);
  }
  return 0;
}


void PISMTillCanHydrology::add_vars_to_output(string /*keyword*/, map<string,NCSpatialVariable> &result) {
  result["bwattillcan"] = W.get_metadata();
  IceModelVec2S tmp;
  tmp.create(grid, "bwptillcan", false);
  tmp.set_attrs("diagnostic",
                       "pressure of water in subglacial layer",
                       "Pa", "");
  tmp.set_attr("valid_min", 0.0);
  result["bwptillcan"] = tmp.get_metadata();
  // destructor called on tmp when we go out of scope here
}


PetscErrorCode PISMTillCanHydrology::define_variables(set<string> vars, const PIO &nc,
                                                 PISM_IO_Type nctype) {
  PetscErrorCode ierr;
  if (set_contains(vars, "bwattillcan")) {
    ierr = W.define(nc, nctype); CHKERRQ(ierr);
  }
  if (set_contains(vars, "bwptillcan")) {
    IceModelVec2S tmp;
    ierr = tmp.create(grid, "bwptillcan", false); CHKERRQ(ierr);
    ierr = tmp.set_attrs("diagnostic",
                     "pressure of water in subglacial layer",
                     "Pa", ""); CHKERRQ(ierr);
    ierr = tmp.set_attr("valid_min", 0.0); CHKERRQ(ierr);
    ierr = tmp.define(nc, nctype); CHKERRQ(ierr);
    // destructor called on tmp when we go out of scope here
  }
  return 0;
}


PetscErrorCode PISMTillCanHydrology::write_variables(set<string> vars, const PIO &nc) {
  PetscErrorCode ierr;
  if (set_contains(vars, "bwattillcan")) {
    ierr = W.write(nc); CHKERRQ(ierr);
  }
  if (set_contains(vars, "bwptillcan")) {
    IceModelVec2S tmp;
    ierr = tmp.create(grid, "bwptillcan", false); CHKERRQ(ierr);
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


PetscErrorCode PISMTillCanHydrology::water_layer_thickness(IceModelVec2S &result) {
  PetscErrorCode ierr = W.copy_to(result); CHKERRQ(ierr);
  return 0;
}


//! Computes pressure diagnostically.
/*!
  \f[ P = \lambda P_o \max\{1,W / W_{crit}\} \f]
where \f$\lambda\f$=till_pw_fraction, \f$P_o = \rho_i g H\f$, \f$W_{crit}\f$=bwat_max.
 */
PetscErrorCode PISMTillCanHydrology::water_pressure(IceModelVec2S &result) {
  PetscErrorCode ierr;
  ierr = result.set_attrs("diagnostic",
                     "pressure of water in subglacial layer",
                     "Pa", ""); CHKERRQ(ierr);
  ierr = result.set_attr("valid_min", 0.0); CHKERRQ(ierr);
  ierr = check_W_bounds(); CHKERRQ(ierr); // check:  W \le bwat_max = Wcrit
  ierr = result.copy_from(W); CHKERRQ(ierr);
  ierr = result.scale(1.0/config.get("bwat_max")); CHKERRQ(ierr); // = max{1,W/Wcrit}
  ierr = result.multiply_by((*thk)); CHKERRQ(ierr); // = H max{1,W/Wcrit}
  ierr = result.scale(config.get("ice_density") * config.get("standard_gravity")); CHKERRQ(ierr);
  ierr = result.scale(config.get("till_pw_fraction")); CHKERRQ(ierr); // P = lambda rhoi g H max{1,W/Wcrit}
  return 0;
}


/*!
Checks \f$0 \le W \le W_{crit} =\f$bwat_max.
 */
PetscErrorCode PISMTillCanHydrology::check_W_bounds() {
  PetscErrorCode ierr;
  PetscReal bwat_max = config.get("bwat_max");
  ierr = W.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (W(i,j) < 0.0) {
        PetscPrintf(grid.com,
           "PISMTillCanHydrology ERROR: disallowed negative subglacial water layer thickness W(i,j) = %.6f m\n"
           "            at (i,j)=(%d,%d)\n"
           "ENDING ... \n\n", W(i,j),i,j);
        PISMEnd();
      }
      if (W(i,j) > bwat_max) {
        PetscPrintf(grid.com,
           "PISMTillCanHydrology ERROR: subglacial water layer thickness W(i,j) = %.6f m exceeds\n"
           "            bmelt_max = %.6f at (i,j)=(%d,%d)\n"
           "ENDING ... \n\n", W(i,j),bwat_max,i,j);
        PISMEnd();
      }
    }
  }
  ierr = W.end_access(); CHKERRQ(ierr);
  return 0;
}


//! Update the water thickness from bmelt input, bwat_max, and decay rate.
/*!
Solves on explicit (forward Euler) step of the integration
  \f[ \frac{dW}{dt} = \text{bmelt} - C \f]
but subject to the inequalities
  \f[ 0 \le W \le W_{crit} \f]
where \f$C=\f$bwat_decay_rate and \f$W_{crit}\f$=bwat_max. 
 */
PetscErrorCode PISMTillCanHydrology::update(PetscReal icet, PetscReal icedt) {
  // if asked for the identical time interval as last time, then do nothing
  if ((fabs(icet - t) < 1e-6) && (fabs(icedt - dt) < 1e-6))
    return 0;
  t = icet;
  dt = icedt;

  PetscErrorCode ierr;
  PetscReal bwat_max        = config.get("bwat_max"),
            bwat_decay_rate = config.get("bwat_decay_rate");
  ierr = W.begin_access(); CHKERRQ(ierr);
  ierr = bmelt->begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      W(i,j) = W(i,j) + ((*bmelt)(i,j) - bwat_decay_rate) * icedt;
      W(i,j) = PetscMax(0.0, PetscMin(bwat_max, W(i,j)) );
    }
  }
  ierr = W.end_access(); CHKERRQ(ierr);
  ierr = bmelt->end_access(); CHKERRQ(ierr);
  return 0;
}


PISMDiffusebwatHydrology::PISMDiffusebwatHydrology(IceGrid &g, const NCConfigVariable &conf)
    : PISMTillCanHydrology(g, conf, true)
{
  if (allocateWnew() != 0) {
    PetscPrintf(grid.com,
      "PISM ERROR: allocation of Wnew failed in PISMDiffusebwatHydrology constructor.\n");
    PISMEnd();
  }
}


PetscErrorCode PISMDiffusebwatHydrology::init(PISMVars &vars) {
  PetscErrorCode ierr;
  ierr = PISMTillCanHydrology::init(vars); CHKERRQ(ierr);
  ierr = verbPrintf(2, grid.com,
    "  using the diffusive water layer variant ...\n"); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode PISMDiffusebwatHydrology::allocateWnew() {
  PetscErrorCode ierr;
  // also need temporary space during update
  //FIXME: shouldn't I be able to do this?  gives error
  //   "makes no sense to communicate ghosts for GLOBAL IceModelVec! (has name='Wnew-internal')!"
  //ierr = Wnew.create(grid, "Wnew-internal", false); CHKERRQ(ierr);
  ierr = Wnew.create(grid, "Wnew_internal", true, 1); CHKERRQ(ierr);
  ierr = Wnew.set_attrs("internal",
                     "new thickness of subglacial water layer during update",
                     "m", ""); CHKERRQ(ierr);
  ierr = Wnew.set_attr("valid_min", 0.0); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode PISMDiffusebwatHydrology::update(PetscReal icet, PetscReal icedt) {
  // if asked for the identical time interval as last time, then do nothing
  if ((fabs(icet - t) < 1e-6) && (fabs(icedt - dt) < 1e-6))
    return 0;
  t = icet;
  dt = icedt;

  PetscErrorCode ierr;
  const PetscReal L = config.get("bwat_diffusion_distance");
  if (L <= 0.0)  {
    ierr = PISMTillCanHydrology::update(icet,icedt); CHKERRQ(ierr);
    return 0;
  }

  const PetscReal
    diffusion_time  = config.get("bwat_diffusion_time", "years", "seconds"), // convert to seconds
    bwat_max        = config.get("bwat_max"),
    bwat_decay_rate = config.get("bwat_decay_rate"),
    K               = L * L / (2.0 * diffusion_time);
  PetscReal hdt;
  PetscInt NN;
  hdt = (1.0 / (grid.dx*grid.dx)) + (1.0 / (grid.dy*grid.dy));
  hdt = 1.0 / (2.0 * K * hdt);
  NN = ceil(dt / hdt);

  if (NN > 1) {
    verbPrintf(2,grid.com,
      "PISMDiffusebwatHydrology WARNING: more than one time step per ice dynamics time step\n"
      "   ... NN = %d > 1 ... THIS IS BELIEVED TO BE RARE\n",NN);
  }

  // assume ghosts are valid at start of time step

  hdt = dt / NN;
  PetscReal  Rx = K * dt / (grid.dx * grid.dx),
             Ry = K * dt / (grid.dy * grid.dy),
             oneM4R = 1.0 - 2.0 * Rx - 2.0 * Ry;
  for (PetscInt n=0; n<NN; ++n) {
    ierr = W.begin_access(); CHKERRQ(ierr);
    ierr = Wnew.begin_access(); CHKERRQ(ierr);
    ierr = bmelt->begin_access(); CHKERRQ(ierr);
    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
        Wnew(i,j) = oneM4R * W(i,j) + Rx * (W(i+1,j  ) + W(i-1,j  ))
                                    + Ry * (W(i  ,j+1) + W(i  ,j-1));
        Wnew(i,j) = Wnew(i,j) + ((*bmelt)(i,j) - bwat_decay_rate) * icedt;
        Wnew(i,j) = PetscMax(0.0, PetscMin(bwat_max, Wnew(i,j)) );
      }
    }
    ierr = W.end_access(); CHKERRQ(ierr);
    ierr = Wnew.end_access(); CHKERRQ(ierr);
    ierr = bmelt->end_access(); CHKERRQ(ierr);

    ierr = Wnew.beginGhostComm(W); CHKERRQ(ierr);
    ierr = Wnew.endGhostComm(W); CHKERRQ(ierr);
  }
  return 0;
}


PISMDistributedHydrology::PISMDistributedHydrology(IceGrid &g, const NCConfigVariable &conf)
    : PISMHydrology(g, conf)
{
    bed   = NULL;
    thk   = NULL;
    usurf  = NULL;
    bmelt  = NULL;

    if (allocate() != 0) {
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
    K     = 1.0e-2;     // m s-1;  want Kmax or Kmin according to W > Wr
    Aglen = 3.1689e-24; // Pa-3 s-1; ice softness
    nglen = 3.0;
    Wr    = 1.0;        // m
    E0    = 1.0;        // m; what is optimal?
    Y0    = 0.001;      // m; regularization

    c0    = K / (fresh_water_density * standard_gravity); // constant in velocity formula
}


PetscErrorCode PISMDistributedHydrology::allocate() {
  PetscErrorCode ierr;

  // model state variables; need ghosts
  ierr = W.create(grid, "bwatdistributed", true, 1); CHKERRQ(ierr);
  ierr = W.set_attrs("model_state",
                     "thickness of subglacial water layer",
                     "m", ""); CHKERRQ(ierr);
  ierr = W.set_attr("valid_min", 0.0); CHKERRQ(ierr);
  ierr = P.create(grid, "bwpdistributed", true, 1); CHKERRQ(ierr);
  ierr = P.set_attrs("model_state",
                     "pressure of water in subglacial layer",
                     "Pa", ""); CHKERRQ(ierr);
  ierr = P.set_attr("valid_min", 0.0); CHKERRQ(ierr);

  // auxiliary variables which NEED ghosts
  ierr = psi.create(grid, "hydraulic_potential", true, 1); CHKERRQ(ierr);
  ierr = psi.set_attrs("internal",
                       "hydraulic potential of water in subglacial layer",
                       "Pa", ""); CHKERRQ(ierr);
  ierr = known.create(grid, "known_hydro_mask", true, 1); CHKERRQ(ierr);
  ierr = known.set_attrs("internal",
                       "mask for where subglacial hydrology state is known",
                       "", ""); CHKERRQ(ierr);
  ierr = Wstag.create(grid, "W_staggered", true, 1); CHKERRQ(ierr);
  ierr = Wstag.set_attrs("internal",
                     "cell face-centered (staggered) values of water layer thickness",
                     "m", ""); CHKERRQ(ierr);
  ierr = Qstag.create(grid, "advection_flux", true, 1); CHKERRQ(ierr);
  ierr = Qstag.set_attrs("internal",
                     "cell face-centered (staggered) components of advective subglacial water flux",
                     "m2 s-1", ""); CHKERRQ(ierr);

  // auxiliary variables which do not need ghosts
  ierr = Po.create(grid, "ice_overburden_pressure", false); CHKERRQ(ierr);
  ierr = Po.set_attrs("internal",
                      "ice overburden pressure seen by subglacial water layer",
                      "Pa", ""); CHKERRQ(ierr);
  ierr = Po.set_attr("valid_min", 0.0); CHKERRQ(ierr);
  ierr = cbase.create(grid, "ice_sliding_speed", false); CHKERRQ(ierr);
  ierr = cbase.set_attrs("internal",
                         "ice sliding speed seen by subglacial water layer",
                         "m s-1", ""); CHKERRQ(ierr);
  ierr = cbase.set_attr("valid_min", 0.0); CHKERRQ(ierr);
  ierr = V.create(grid, "water_velocity", false); CHKERRQ(ierr);
  ierr = V.set_attrs("internal",
                     "cell face-centered (staggered) components of water velocity in subglacial water layer",
                     "m s-1", ""); CHKERRQ(ierr);

  // temporaries during update; do not need ghosts
  ierr = Wnew.create(grid, "Wnew_internal", false); CHKERRQ(ierr);
  ierr = Wnew.set_attrs("internal",
                     "new thickness of subglacial water layer during update",
                     "m", ""); CHKERRQ(ierr);
  ierr = Wnew.set_attr("valid_min", 0.0); CHKERRQ(ierr);
  ierr = Pnew.create(grid, "Pnew_internal", false); CHKERRQ(ierr);
  ierr = Pnew.set_attrs("internal",
                     "new subglacial water pressure during update",
                     "Pa", ""); CHKERRQ(ierr);
  ierr = Pnew.set_attr("valid_min", 0.0); CHKERRQ(ierr);

  return 0;
}


PetscErrorCode PISMDistributedHydrology::init(PISMVars &vars, PISMStressBalance &sb) {
  PetscErrorCode ierr;

  variables = &vars;
  stressbalance = &sb;

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


PetscErrorCode PISMDistributedHydrology::water_layer_thickness(IceModelVec2S &result) {
  PetscErrorCode ierr = W.copy_to(result); CHKERRQ(ierr);
  return 0;
}


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

  ierr = W.begin_access(); CHKERRQ(ierr);
  ierr = P.begin_access(); CHKERRQ(ierr);
  ierr = Po.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (W(i,j) < 0.0) {
        PetscPrintf(grid.com,
           "PISM ERROR: disallowed negative subglacial water layer thickness W(i,j) = %.6f m\n"
           "            at (i,j)=(%d,%d)\n"
           "ENDING ... \n\n", W(i,j),i,j);
        PISMEnd();
      }
      if (P(i,j) < 0.0) {
        PetscPrintf(grid.com,
           "PISM ERROR: disallowed negative subglacial water pressure P(i,j) = %.6f Pa\n"
           "            at (i,j)=(%d,%d)\n"
           "ENDING ... \n\n", P(i,j),i,j);
        PISMEnd();
      }
      if (P(i,j) > Po(i,j)) {
        PetscPrintf(grid.com,
           "PISM ERROR: subglacial water pressure P(i,j) = %.6f Pa exceeds\n"
           "            overburden pressure Po(i,j) = %.6f Pa at (i,j)=(%d,%d);\n"
           "            (not allowed)\n"
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


//! Get the advection velocity V at the center of cell edges.
/*!
Computes the advection velocity \f$V=V(\nabla P,\nabla b)\f$ on the
staggered (face-centered) grid.  If V = (alpha,beta) in components
then we have <code> result(i,j,0) = alpha(i+1/2,j) </code> and
<code> result(i,j,1) = beta(i,j+1/2) </code>
 */
PetscErrorCode PISMDistributedHydrology::velocity_staggered(IceModelVec2Stag &result) {
  PetscErrorCode ierr;
  PetscReal dbdx, dbdy, dPdx, dPdy;

  ierr = P.begin_access(); CHKERRQ(ierr);
  ierr = bed->begin_access(); CHKERRQ(ierr);
  ierr = result.begin_access(); CHKERRQ(ierr);
  for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
      dPdx = (P(i+1,j  ) - P(i,j)) / grid.dx;
      dPdy = (P(i  ,j+1) - P(i,j)) / grid.dy;
      dbdx = ((*bed)(i+1,j  ) - (*bed)(i,j)) / grid.dx;
      dbdy = ((*bed)(i  ,j+1) - (*bed)(i,j)) / grid.dy;
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
PetscErrorCode PISMDistributedHydrology::water_thickness_staggered(IceModelVec2Stag &result) {
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
PetscErrorCode PISMDistributedHydrology::advective_fluxes(IceModelVec2Stag &result) {
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


//! Get the hydraulic potential from bedrock topography and current state variables.
/*!
Computes \f$\psi = P + \rho_w g (b + W)\f$.
 */
PetscErrorCode PISMDistributedHydrology::hydraulic_potential(IceModelVec2S &result) {
  PetscErrorCode ierr;
  ierr = W.begin_access(); CHKERRQ(ierr);
  ierr = P.begin_access(); CHKERRQ(ierr);
  ierr = bed->begin_access(); CHKERRQ(ierr);
  ierr = result.begin_access(); CHKERRQ(ierr);
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


//! Compute a mask with states 0=unknown (active) location, 1=icefree, 2=floating.
PetscErrorCode PISMDistributedHydrology::known_state_mask(IceModelVec2Int &result) {
  PetscErrorCode ierr;
  PetscReal hij, bij, Hfloat,
            rr = ice_density / fresh_water_density,
            cc = 1.0 - rr;  // surface elevation is cc times thickness for floating ice
  ierr = usurf->begin_access(); CHKERRQ(ierr);
  ierr = bed->begin_access(); CHKERRQ(ierr);
  ierr = result.begin_access(); CHKERRQ(ierr);
  for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
      hij = (*usurf)(i,j);   // surface elevation
      bij = (*bed)(i,j);     // bedrock elevation
      // fixme: use PISM's sea level
      if ((hij > 0.0) && (hij < bij + 1.0))
        result(i,j) = 1.0;   // known: icefree
      else {
        Hfloat = hij / cc;
        if (ice_density * Hfloat < - sea_water_density * bij)
          result(i,j) = 2.0; // known: float
        else
          result(i,j) = 0.0; // not known, i.e. subglacial aquifer location
      }
    }
  }
  ierr = usurf->end_access(); CHKERRQ(ierr);
  ierr = bed->end_access(); CHKERRQ(ierr);
  ierr = result.end_access(); CHKERRQ(ierr);
  return 0;
}


//! Update the overburden pressure Po and the sliding speed |v_b| from ice quantities.
/*!
Accesses thk from PISMVars, which points into IceModel.  Calls a PISMStressBalance
method to get the vector basal velocity of the ice, and then computes the magnitude
of that.

fixme:  Is taking Ubase from the PISMStressBalance the correct method?
 */
PetscErrorCode PISMDistributedHydrology::update_ice_functions(IceModelVec2S &result_Po,
                                                   IceModelVec2S &result_cbase) {
  PetscErrorCode ierr;
  IceModelVec2V* Ubase; // ice sliding velocity
  // Po = rho_i g H
  ierr = result_Po.copy_from(*thk); CHKERRQ(ierr);
  ierr = result_Po.scale(ice_density * standard_gravity); CHKERRQ(ierr);
  // cbase = |v_b|
  ierr = stressbalance->get_2D_advective_velocity(Ubase); CHKERRQ(ierr);
  ierr = Ubase->magnitude(result_cbase); CHKERRQ(ierr);
  return 0;
}


//! Compute dt = min[t_end-t_current, dtCFL, dtDIFFW, dtDIFFP], the adaptive time step.
PetscErrorCode PISMDistributedHydrology::adaptive_time_step(PetscReal t_current, PetscReal t_end, 
                                                 PetscReal &dt_result) {
  PetscErrorCode ierr;
  PetscReal dtmax, dtCFL, dtDIFFW, dtDIFFP, maxW, maxH;
  PetscReal tmp[2];

  // fixme?: dtCFL can be infinity if velocity is zero because P and b are constant
  // Matlab: dtCFL = 0.5 / (max(max(abs(alphV)))/dx + max(max(abs(betaV)))/dy);
  ierr = V.absmaxcomponents(tmp); CHKERRQ(ierr);
  dtCFL = 0.5 / (tmp[0]/grid.dx + tmp[1]/grid.dy);

  // Matlab: maxW = max(max(max(Wea)),max(max(Wno))) + 0.001;
  ierr = Wstag.absmaxcomponents(tmp); CHKERRQ(ierr);
  maxW = PetscMax(tmp[0],tmp[1]) + 0.001;
  // Matlab: dtDIFFW = 0.25 / (p.K * maxW * (1/dx^2 + 1/dy^2));
  dtDIFFW = 1.0/(grid.dx*grid.dx) + 1.0/(grid.dy*grid.dy);
  dtDIFFW = 0.25 / (K * maxW * dtDIFFW);

  // Matlab: dtDIFFP = (p.rhow * p.E0 / (p.rhoi * maxH)) * dtDIFFW;
  ierr = thk->max(maxH); CHKERRQ(ierr);
  maxH += 2.0 * E0; // regularized: forces dtDIFFP < dtDIFFW
  dtDIFFP = (fresh_water_density * E0 / (ice_density * maxH)) * dtDIFFW;

  dtmax = (1.0/12.0) * secpera;  // fixme?: need better dtmax than fixed at 1 month

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
PetscErrorCode PISMDistributedHydrology::update(PetscReal icet, PetscReal icedt) {
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

    ierr = hydraulic_potential(psi); CHKERRQ(ierr);
    ierr = psi.beginGhostComm(); CHKERRQ(ierr);

    ierr = velocity_staggered(V); CHKERRQ(ierr);

    ierr = water_thickness_staggered(Wstag); CHKERRQ(ierr);
    ierr = Wstag.beginGhostComm(); CHKERRQ(ierr);

    ierr = known_state_mask(known); CHKERRQ(ierr);
    ierr = known.beginGhostComm(); CHKERRQ(ierr);

    // to get Qstag, W needs valid ghosts
    ierr = advective_fluxes(Qstag); CHKERRQ(ierr);
    ierr = Qstag.beginGhostComm(); CHKERRQ(ierr);

    ierr = adaptive_time_step(ht, t+dt, hdt); CHKERRQ(ierr);

    ierr = psi.endGhostComm(); CHKERRQ(ierr);
    ierr = Wstag.endGhostComm(); CHKERRQ(ierr);
    ierr = known.endGhostComm(); CHKERRQ(ierr);

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

