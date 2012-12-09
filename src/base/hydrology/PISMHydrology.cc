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


PetscErrorCode PISMHydrology::init(PISMVars &vars) {
  PetscErrorCode ierr;
  ierr = verbPrintf(4, grid.com,
    "entering initializer for base class PISMHydrology ...\n"); CHKERRQ(ierr);

  variables = &vars;

  thk = dynamic_cast<IceModelVec2S*>(vars.get("thk"));
  if (thk == NULL) SETERRQ(grid.com, 1, "thk is not available to PISMHydrology");

  bed = dynamic_cast<IceModelVec2S*>(vars.get("topg"));
  if (bed == NULL) SETERRQ(grid.com, 1, "topg is not available to PISMHydrology");

  bmelt = dynamic_cast<IceModelVec2S*>(vars.get("bmelt"));
  if (bmelt == NULL) SETERRQ(grid.com, 1, "bmelt is not available to PISMHydrology");

  mask = dynamic_cast<IceModelVec2Int*>(vars.get("mask"));
  if (mask == NULL) SETERRQ(grid.com, 1, "mask is not available to PISMHydrology");

  return 0;
}


void PISMHydrology::get_diagnostics(map<string, PISMDiagnostic*> &dict) {
  dict["bwp"] = new PISMHydrology_bwp(this, grid, *variables);
}


//! Compute the water input rate into the basal hydrology layer according to configuration and mask.
PetscErrorCode PISMHydrology::get_input_rate(IceModelVec2S &result) {
  PetscErrorCode ierr;
  bool      use_const   = config.get_flag("hydrology_use_const_bmelt");
  PetscReal const_bmelt = config.get("hydrology_const_bmelt");
  ierr = bmelt->begin_access(); CHKERRQ(ierr);
  ierr = mask->begin_access(); CHKERRQ(ierr);
  ierr = result.begin_access(); CHKERRQ(ierr);
  MaskQuery m(*mask);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (m.icy(i, j))
        result(i,j) = (use_const) ? const_bmelt : (*bmelt)(i,j);
      else
        result(i,j) = 0.0;
    }
  }
  ierr = bmelt->end_access(); CHKERRQ(ierr);
  ierr = mask->end_access(); CHKERRQ(ierr);
  ierr = result.end_access(); CHKERRQ(ierr);
  return 0;
}


PISMTillCanHydrology::PISMTillCanHydrology(IceGrid &g, const NCConfigVariable &conf,
                                           bool Whasghosts)
    : PISMHydrology(g, conf)
{
    if (allocate(Whasghosts) != 0) {
      PetscPrintf(grid.com,
        "PISM ERROR: allocation failed in PISMTillCanHydrology constructor.\n");
      PISMEnd();
    }
}


PetscErrorCode PISMTillCanHydrology::allocate(bool Whasghosts) {
  PetscErrorCode ierr;
  // workspace
  ierr = input.create(grid, "input_hydro", false); CHKERRQ(ierr);
  ierr = input.set_attrs("internal",
                         "workspace for input into subglacial water layer",
                         "m s-1", ""); CHKERRQ(ierr);
  // model state variables
  if (Whasghosts) {
    ierr = W.create(grid, "bwat", true, 1); CHKERRQ(ierr);
  } else {
    ierr = W.create(grid, "bwat", false); CHKERRQ(ierr);
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
  ierr = PISMHydrology::init(vars); CHKERRQ(ierr);
  // initialize water layer thickness from the context if present, otherwise zero
  IceModelVec2S *W_input = dynamic_cast<IceModelVec2S*>(vars.get("bwat"));
  if (W_input != NULL) {
    ierr = W.copy_from(*W_input); CHKERRQ(ierr);
    // FIXME: what about regrid case under -boot_file?
  } else {
    ierr = W.set(0.0); CHKERRQ(ierr);
  }
  if (vars.get("bwat") == NULL) { // since init() will get called twice, *we*
                                  //   might have already added "bwat"
    ierr = vars.add(W); CHKERRQ(ierr);
  }
  return 0;
}


void PISMTillCanHydrology::add_vars_to_output(string /*keyword*/, map<string,NCSpatialVariable> &result) {
  result["bwat"] = W.get_metadata();
}


PetscErrorCode PISMTillCanHydrology::define_variables(set<string> vars, const PIO &nc,
                                                 PISM_IO_Type nctype) {
  PetscErrorCode ierr;
  if (set_contains(vars, "bwat")) {
    ierr = W.define(nc, nctype); CHKERRQ(ierr);
  }
  return 0;
}


PetscErrorCode PISMTillCanHydrology::write_variables(set<string> vars, const PIO &nc) {
  PetscErrorCode ierr;
  if (set_contains(vars, "bwat")) {
    ierr = W.write(nc); CHKERRQ(ierr);
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
where \f$\lambda\f$=till_pw_fraction, \f$P_o = \rho_i g H\f$, \f$W_{crit}\f$=hydrology_bwat_max.
 */
PetscErrorCode PISMTillCanHydrology::water_pressure(IceModelVec2S &result) {
  PetscErrorCode ierr;

#if (PISM_DEBUG==1)
  ierr = check_W_bounds(); CHKERRQ(ierr); // check:  W \le bwat_max = Wcrit
#endif

  double bwat_max = config.get("hydrology_bwat_max"),
    till_pw_fraction = config.get("till_pw_fraction"),
    standard_gravity = config.get("standard_gravity"),
    ice_density = config.get("ice_density");

  ierr = W.begin_access(); CHKERRQ(ierr);
  ierr = result.begin_access(); CHKERRQ(ierr);
  ierr = thk->begin_access(); CHKERRQ(ierr);

  for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
      result(i,j) = till_pw_fraction * ice_density * standard_gravity * (*thk)(i,j) * (W(i,j) / bwat_max);
    }
  }

  ierr = thk->end_access(); CHKERRQ(ierr);
  ierr = result.end_access(); CHKERRQ(ierr);
  ierr = W.end_access(); CHKERRQ(ierr);
  return 0;
}


/*!
Checks \f$0 \le W \le W_{crit} =\f$hydrology_bwat_max.
 */
PetscErrorCode PISMTillCanHydrology::check_W_bounds() {
  PetscErrorCode ierr;
  PetscReal bwat_max = config.get("hydrology_bwat_max");
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
           "            hydrology_bwat_max = %.6f at (i,j)=(%d,%d)\n"
           "ENDING ... \n\n", W(i,j),bwat_max,i,j);
        PISMEnd();
      }
    }
  }
  ierr = W.end_access(); CHKERRQ(ierr);
  return 0;
}


//! Update the water thickness from input (melt and drainage from ice above), the upper bound on water amount, and the decay rate.
/*!
Solves on explicit (forward Euler) step of the integration
  \f[ \frac{dW}{dt} = \text{bmelt} - C \f]
but subject to the inequalities
  \f[ 0 \le W \le W_{crit} \f]
where \f$C=\f$hydrology_bwat_decay_rate and \f$W_{crit}\f$=hydrology_bwat_max. 
 */
PetscErrorCode PISMTillCanHydrology::update(PetscReal icet, PetscReal icedt) {
  // if asked for the identical time interval as last time, then do nothing
  if ((fabs(icet - t) < 1e-6) && (fabs(icedt - dt) < 1e-6))
    return 0;
  t = icet;
  dt = icedt;

  PetscErrorCode ierr;

  ierr = get_input_rate(input); CHKERRQ(ierr);

  PetscReal bwat_max        = config.get("hydrology_bwat_max"),
            bwat_decay_rate = config.get("hydrology_bwat_decay_rate");
  ierr = W.begin_access(); CHKERRQ(ierr);
  ierr = input.begin_access(); CHKERRQ(ierr);
  ierr = mask->begin_access(); CHKERRQ(ierr);
  MaskQuery m(*mask);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (m.grounded_ice(i, j)) {
        W(i,j) = W(i,j) + (input(i,j) - bwat_decay_rate) * icedt;
        W(i,j) = PetscMax(0.0, PetscMin(bwat_max, W(i,j)) );
      } else if (m.ice_free_land(i,j)) {
        W(i,j) = 0.0;
      } else { // floating or ocean cases
        W(i,j) = bwat_max;
      }
    }
  }
  ierr = W.end_access(); CHKERRQ(ierr);
  ierr = input.end_access(); CHKERRQ(ierr);
  ierr = mask->end_access(); CHKERRQ(ierr);
  return 0;
}


PISMDiffuseOnlyHydrology::PISMDiffuseOnlyHydrology(IceGrid &g, const NCConfigVariable &conf)
    : PISMTillCanHydrology(g, conf, true)
{
  if (allocateWnew() != 0) {
    PetscPrintf(grid.com,
      "PISM ERROR: allocation of Wnew failed in PISMDiffuseOnlyHydrology constructor.\n");
    PISMEnd();
  }
}


PetscErrorCode PISMDiffuseOnlyHydrology::allocateWnew() {
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


PetscErrorCode PISMDiffuseOnlyHydrology::init(PISMVars &vars) {
  PetscErrorCode ierr;
  ierr = PISMTillCanHydrology::init(vars); CHKERRQ(ierr);
  ierr = verbPrintf(2, grid.com,
    "  using the diffusive water layer variant ...\n"); CHKERRQ(ierr);
  return 0;
}


//! Explicit time step for diffusion of subglacial water layer bwat.
/*!
See equation (11) in \ref BBssasliding , namely
  \f[W_t = K \nabla^2 W.\f]
The diffusion constant \f$K\f$ is chosen so that the fundamental solution (Green's
function) of this equation has standard deviation \f$\sigma=L\f$ at time t=\c diffusion_time.
Note that \f$2 \sigma^2 = 4 K t\f$.

The time step restriction for the explicit method for this equation is believed
to be so rare that if it is triggered there is a stdout warning.
 */
PetscErrorCode PISMDiffuseOnlyHydrology::update(PetscReal icet, PetscReal icedt) {
  // if asked for the identical time interval as last time, then do nothing
  if ((fabs(icet - t) < 1e-6) && (fabs(icedt - dt) < 1e-6))
    return 0;
  t = icet;
  dt = icedt;

  PetscErrorCode ierr;
  const PetscReal L = config.get("hydrology_bwat_diffusion_distance");
  if (L <= 0.0)  {
    ierr = PISMTillCanHydrology::update(icet,icedt); CHKERRQ(ierr);
    return 0;
  }

  const PetscReal
    diffusion_time  = config.get("hydrology_bwat_diffusion_time", "years", "seconds"), // convert to seconds
    bwat_max        = config.get("hydrology_bwat_max"),
    bwat_decay_rate = config.get("hydrology_bwat_decay_rate"),
    K               = L * L / (2.0 * diffusion_time);

  PetscReal hdt;
  PetscInt NN;
  hdt = (1.0 / (grid.dx*grid.dx)) + (1.0 / (grid.dy*grid.dy));
  hdt = 1.0 / (2.0 * K * hdt);
  NN = ceil(dt / hdt);

  if (NN > 1) {
    verbPrintf(2,grid.com,
      "PISMDiffuseOnlyHydrology WARNING: more than one time step per ice dynamics time step\n"
      "   ... NN = %d > 1 ... THIS IS BELIEVED TO BE RARE\n",NN);
  }

  ierr = get_input_rate(input); CHKERRQ(ierr);

  hdt = dt / NN;
  PetscReal  Rx = K * dt / (grid.dx * grid.dx),
             Ry = K * dt / (grid.dy * grid.dy),
             oneM4R = 1.0 - 2.0 * Rx - 2.0 * Ry;
  for (PetscInt n=0; n<NN; ++n) {
    // time-splitting: first, Euler step on source terms
    ierr = W.begin_access(); CHKERRQ(ierr);
    ierr = input.begin_access(); CHKERRQ(ierr);
    ierr = mask->begin_access(); CHKERRQ(ierr);
    MaskQuery m(*mask);
    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
        if (m.grounded_ice(i, j)) {
          W(i,j) = W(i,j) + (input(i,j) - bwat_decay_rate) * icedt;
          W(i,j) = PetscMax(0.0, PetscMin(bwat_max, W(i,j)) );
        } else if (m.ice_free_land(i,j)) {
          W(i,j) = 0.0;
        } else { // floating or ocean cases
          W(i,j) = bwat_max;
        }
      }
    }
    ierr = W.end_access(); CHKERRQ(ierr);
    ierr = input.end_access(); CHKERRQ(ierr);
    ierr = mask->end_access(); CHKERRQ(ierr);

    // valid ghosts for diffusion below
    ierr = W.beginGhostComm(); CHKERRQ(ierr);
    ierr = W.endGhostComm(); CHKERRQ(ierr);

    // time-splitting: second, diffusion by first-order explicit
    ierr = W.begin_access(); CHKERRQ(ierr);
    ierr = Wnew.begin_access(); CHKERRQ(ierr);
    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
        Wnew(i,j) = oneM4R * W(i,j) + Rx * (W(i+1,j  ) + W(i-1,j  ))
                                    + Ry * (W(i  ,j+1) + W(i  ,j-1));
        // no check of bounds here because maximum principle applies to step
      }
    }
    ierr = W.end_access(); CHKERRQ(ierr);
    ierr = Wnew.end_access(); CHKERRQ(ierr);

    // maybe unneeded: valid ghosts for future actions
    ierr = Wnew.beginGhostComm(W); CHKERRQ(ierr);
    ierr = Wnew.endGhostComm(W); CHKERRQ(ierr);
  }
  return 0;
}


PISMHydrology_bwp::PISMHydrology_bwp(PISMHydrology *m, IceGrid &g, PISMVars &my_vars)
  : PISMDiag<PISMHydrology>(m, g, my_vars) {

  // set metadata:
  vars[0].init_2d("bwp", grid);

  set_attrs("pressure of water in subglacial layer", "",
            "Pa", "Pa", 0);
}


PetscErrorCode PISMHydrology_bwp::compute(IceModelVec* &output) {
  PetscErrorCode ierr;

  IceModelVec2S *result = new IceModelVec2S;
  ierr = result->create(grid, "bwp", false); CHKERRQ(ierr);
  ierr = result->set_metadata(vars[0], 0); CHKERRQ(ierr);
  result->write_in_glaciological_units = true;

  ierr = model->water_pressure(*result); CHKERRQ(ierr);

  output = result;
  return 0;
}

