// Copyright (C) 2012-2013 PISM Authors
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

#include "PISMHydrology.hh"
#include "PISMVars.hh"
#include "pism_options.hh"
#include "Mask.hh"
#include "PISMStressBalance.hh"


PISMDistributedHydrology::PISMDistributedHydrology(IceGrid &g, const NCConfigVariable &conf,
                                                   PISMStressBalance *sb)
    : PISMRoutingHydrology(g, conf)
{
    stressbalance = sb;
    if (allocate_englacial() != 0) {
      PetscPrintf(grid.com,
        "PISM ERROR: memory allocation failed in PISMDistributedHydrology constructor (englacial).\n");
      PISMEnd();
    }
    if (allocate_pressure() != 0) {
      PetscPrintf(grid.com,
        "PISM ERROR: memory allocation failed in PISMDistributedHydrology constructor (pressure).\n");
      PISMEnd();
    }
}


PetscErrorCode PISMDistributedHydrology::allocate_englacial() {
  PetscErrorCode ierr;

  // additional conserved (mass) variable
  ierr = Wen.create(grid, "enwat", false); CHKERRQ(ierr);
  ierr = Wen.set_attrs("model_state",
                       "effective thickness of englacial water",
                       "m", ""); CHKERRQ(ierr);
  ierr = Wen.set_attr("valid_min", 0.0); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode PISMDistributedHydrology::allocate_pressure() {
  PetscErrorCode ierr;

  // additional variables beyond PISMRoutingHydrology::allocate()
  ierr = P.create(grid, "bwp", true, 1); CHKERRQ(ierr);
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
  ierr = psi.create(grid, "hydraulic_potential", true, 1); CHKERRQ(ierr);
  ierr = psi.set_attrs("internal",
                       "hydraulic potential of water in subglacial layer",
                       "Pa", ""); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode PISMDistributedHydrology::init(PISMVars &vars) {
  PetscErrorCode ierr;
  ierr = verbPrintf(2, grid.com,
    "* Initializing the vanPelt-Bueler distributed (linked-cavities) subglacial hydrology model...\n");
    CHKERRQ(ierr);

  // initialize water layer thickness and wate pressure from the context if present,
  //   otherwise from -i or -boot_file, otherwise with constant values
  bool i_set, bootstrap_set, init_P_from_steady, stripset;
  ierr = PetscOptionsBegin(grid.com, "",
            "Options controlling the 'distributed' subglacial hydrology model", ""); CHKERRQ(ierr);
  {
    ierr = PISMOptionsIsSet("-i", "PISM input file", i_set); CHKERRQ(ierr);
    ierr = PISMOptionsIsSet("-boot_file", "PISM bootstrapping file",
                            bootstrap_set); CHKERRQ(ierr);
    ierr = PISMOptionsIsSet("-init_P_from_steady",
                            "initialize P from formula P(W) which applies in steady state",
                            init_P_from_steady); CHKERRQ(ierr);
    ierr = PISMOptionsReal("-hydrology_null_strip",
                           "set the width, in km, of the strip around the edge of the computational domain in which hydrology is inactivated",
                           stripwidth,stripset); CHKERRQ(ierr);
    if (stripset) stripwidth *= 1.0e3;
  }
  ierr = PetscOptionsEnd(); CHKERRQ(ierr);

  ierr = PISMHydrology::init(vars); CHKERRQ(ierr);

  ierr = PISMRoutingHydrology::init_bwat(vars,i_set,bootstrap_set); CHKERRQ(ierr);

  // prepare for -i or -bootstrap
  string filename;
  int start;
  if (i_set || bootstrap_set) {
    ierr = find_pism_input(filename, bootstrap_set, start); CHKERRQ(ierr);
  }

  // initialize Wen: present or -i file or -bootstrap file or set to constant;
  //   then overwrite by regrid
  IceModelVec2S *Wen_input = dynamic_cast<IceModelVec2S*>(vars.get("enwat"));
  if (Wen_input != NULL) { // a variable called "enwat" is already in context
    ierr = Wen.copy_from(*Wen_input); CHKERRQ(ierr);
  } else if (i_set || bootstrap_set) {
    if (i_set) {
      ierr = Wen.read(filename, start); CHKERRQ(ierr);
    } else {
      ierr = Wen.regrid(filename,
                        config.get("bootstrapping_enwat_value_no_var")); CHKERRQ(ierr);
    }
  } else {
    ierr = Wen.set(config.get("bootstrapping_enwat_value_no_var")); CHKERRQ(ierr);
  }
  ierr = regrid(Wen); CHKERRQ(ierr); //  we could be asked to regrid from file

  // initialize P: present or -i file or -bootstrap file or set to constant;
  //   then overwrite by regrid; then overwrite by -init_P_from_steady
  IceModelVec2S *P_input = dynamic_cast<IceModelVec2S*>(vars.get("bwp"));
  if (P_input != NULL) { // a variable called "bwp" is already in context
    ierr = P.copy_from(*P_input); CHKERRQ(ierr);
  } else if (i_set || bootstrap_set) {
    if (i_set) {
      ierr = P.read(filename, start); CHKERRQ(ierr);
    } else {
      ierr = P.regrid(filename,
                      config.get("bootstrapping_bwp_value_no_var")); CHKERRQ(ierr);
    }
  } else {
    ierr = P.set(config.get("bootstrapping_bwp_value_no_var")); CHKERRQ(ierr);
  }
  ierr = regrid(P); CHKERRQ(ierr); //  we could be asked to regrid from file
  if (init_P_from_steady) { // if so, overwrite all the other stuff
    ierr = P_from_W_steady(P); CHKERRQ(ierr);
  }

  // add variables to the context if not already there
  if (vars.get("enwat") == NULL) {
    ierr = vars.add(Wen); CHKERRQ(ierr);
  }
  if (vars.get("bwp") == NULL) {
    ierr = vars.add(P); CHKERRQ(ierr);
  }
  return 0;
}


void PISMDistributedHydrology::add_vars_to_output(string /*keyword*/, set<string> &result) {
  result.insert("bwat");
  result.insert("enwat");
  result.insert("bwp");
}


PetscErrorCode PISMDistributedHydrology::define_variables(set<string> vars, const PIO &nc,
                                                 PISM_IO_Type nctype) {
  PetscErrorCode ierr;
  if (set_contains(vars, "bwat")) {
    ierr = W.define(nc, nctype); CHKERRQ(ierr);
  }
  if (set_contains(vars, "enwat")) {
    ierr = Wen.define(nc, nctype); CHKERRQ(ierr);
  }
  if (set_contains(vars, "bwp")) {
    ierr = P.define(nc, nctype); CHKERRQ(ierr);
  }
  return 0;
}


PetscErrorCode PISMDistributedHydrology::write_variables(set<string> vars, const PIO &nc) {
  PetscErrorCode ierr;
  if (set_contains(vars, "bwat")) {
    ierr = W.write(nc); CHKERRQ(ierr);
  }
  if (set_contains(vars, "enwat")) {
    ierr = Wen.write(nc); CHKERRQ(ierr);
  }
  if (set_contains(vars, "bwp")) {
    ierr = P.write(nc); CHKERRQ(ierr);
  }
  return 0;
}


void PISMDistributedHydrology::get_diagnostics(map<string, PISMDiagnostic*> &dict) {
  dict["bwprel"] = new PISMHydrology_bwprel(this, grid, *variables);
  dict["effbwp"] = new PISMHydrology_effbwp(this, grid, *variables);
  dict["tillwp"] = new PISMHydrology_tillwp(this, grid, *variables);
  dict["hydroinput"] = new PISMHydrology_hydroinput(this, grid, *variables);
  dict["wallmelt"] = new PISMHydrology_wallmelt(this, grid, *variables);
  dict["bwatvel"] = new PISMRoutingHydrology_bwatvel(this, grid, *variables);
}


//! Copies the Wen state variable which is the modeled effective englacial water thickness.
PetscErrorCode PISMDistributedHydrology::englacial_water_thickness(IceModelVec2S &result) {
  PetscErrorCode ierr;
  ierr = Wen.copy_to(result); CHKERRQ(ierr);
  return 0;
}


//! Copies the P state variable which is the modeled water pressure.
PetscErrorCode PISMDistributedHydrology::subglacial_water_pressure(IceModelVec2S &result) {
  PetscErrorCode ierr;
  ierr = P.copy_to(result); CHKERRQ(ierr);
  return 0;
}


//! Check Wen >= 0 and fails with message if not satisfied.
PetscErrorCode PISMDistributedHydrology::check_Wen_nonnegative() {
  PetscErrorCode ierr;
  ierr = Wen.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (Wen(i,j) < 0.0) {
        PetscPrintf(grid.com,
           "PISM ERROR: disallowed negative englacial effective water layer thickness (enwat)\n"
           "    Wen(i,j) = %.6f m at (i,j)=(%d,%d)\n"
           "ENDING ... \n\n", Wen(i,j),i,j);
        PISMEnd();
      }
    }
  }
  ierr = Wen.end_access(); CHKERRQ(ierr);
  return 0;
}


//! Check bounds on P and fail with message if not satisfied.  Optionally, enforces the upper bound instead of checking it.
/*!
The bounds are \f$0 \le P \le P_o\f$ where \f$P_o\f$ is the overburden pressure.
 */
PetscErrorCode PISMDistributedHydrology::check_P_bounds(bool enforce_upper) {
  PetscErrorCode ierr;

  ierr = overburden_pressure(Pover); CHKERRQ(ierr);

  ierr = P.begin_access(); CHKERRQ(ierr);
  ierr = Pover.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (P(i,j) < 0.0) {
        PetscPrintf(grid.com,
           "PISM ERROR: disallowed negative subglacial water pressure\n"
           "    P = %.6f Pa\n at (i,j)=(%d,%d)\n"
           "ENDING ... \n\n", P(i,j),i,j);
        PISMEnd();
      }
      if (enforce_upper) {
        P(i,j) = PetscMin(P(i,j), Pover(i,j));
      } else if (P(i,j) > Pover(i,j) + 0.001) {
        PetscPrintf(grid.com,
           "PISM ERROR: subglacial water pressure P = %.16f Pa exceeds\n"
           "    overburden pressure Po = %.16f Pa at (i,j)=(%d,%d)\n"
           "ENDING ... \n\n", P(i,j),Pover(i,j),i,j);
        PISMEnd();
      }
    }
  }
  ierr = P.end_access(); CHKERRQ(ierr);
  ierr = Pover.end_access(); CHKERRQ(ierr);
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
  PetscReal CC = config.get("hydrology_cavitation_opening_coefficient") /
                    (config.get("hydrology_creep_closure_coefficient") * config.get("ice_softness")),
            powglen = 1.0 / config.get("Glen_exponent"),
            Wr = config.get("hydrology_roughness_scale"),
            Y0 = config.get("hydrology_lower_bound_creep_regularization"),
            sb, Wratio;
  ierr = overburden_pressure(Pover); CHKERRQ(ierr);
  ierr = W.begin_access(); CHKERRQ(ierr);
  ierr = Pover.begin_access(); CHKERRQ(ierr);
  ierr = cbase.begin_access(); CHKERRQ(ierr);
  ierr = result.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      sb     = pow(CC * cbase(i,j),powglen);
      Wratio = PetscMax(0.0,Wr - W(i,j)) / (W(i,j) + Y0);
      // in cases where steady state is actually possible this will
      //   come out positive, but otherwise we should get underpressure P=0,
      //   and that is what it yields
      result(i,j) = PetscMax( 0.0,Pover(i,j) - sb * pow(Wratio,powglen) );
    }
  }
  ierr = W.end_access(); CHKERRQ(ierr);
  ierr = Pover.end_access(); CHKERRQ(ierr);
  ierr = cbase.end_access(); CHKERRQ(ierr);
  ierr = result.end_access(); CHKERRQ(ierr);
  return 0;
}


//! Update the the sliding speed |v_b| from ice quantities.
/*!
Calls a PISMStressBalance method to get the vector basal velocity of the ice,
and then computes the magnitude of that.
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
                  PetscReal t_current, PetscReal t_end, PetscReal maxKW,
                  PetscReal &dt_result,
                  PetscReal &maxV_result, PetscReal &maxD_result,
                  PetscReal &PtoCFLratio) {
  PetscErrorCode ierr;
  PetscReal dtCFL, dtDIFFW, dtDIFFP;

  ierr = adaptive_for_W_evolution(t_current,t_end, maxKW,
              dt_result,maxV_result,maxD_result,dtCFL,dtDIFFW); CHKERRQ(ierr);

  const PetscReal phisum   = config.get("hydrology_englacial_porosity")
                               + config.get("hydrology_regularizing_porosity");
  dtDIFFP = 2.0 * phisum * dtDIFFW;

  // dt = min([te-t dtmax dtCFL dtDIFFW dtDIFFP]);
  dt_result = PetscMin(dt_result, dtDIFFP);

  if (dtDIFFP > 0.0)
    PtoCFLratio = PetscMax(1.0, dtCFL / dtDIFFP);
  else
    PtoCFLratio = 1.0;

  ierr = verbPrintf(3,grid.com,
                    "   [%.5e  %.7f  %.6f  %.9f  -->  dt = %.9f (a)  at  t = %.6f (a)]\n",
                    grid.convert(maxV_result, "m/second", "m/year"),
                    grid.convert(dtCFL,       "seconds",  "years"),
                    grid.convert(dtDIFFW,     "seconds",  "years"),
                    grid.convert(dtDIFFP,     "seconds",  "years"),
                    grid.convert(dt_result,   "seconds",  "years"),
                    grid.convert(t_current,   "seconds",  "years")); CHKERRQ(ierr);
  return 0;
}


//! Update englacial amount Wen using the new subglacial pressure and the new total water amount Wnew_tot.  At exit, Wnew_tot is the subglacial water.
/*!
Given the pressure Pnew at a location, the key function here is of the form
        \f[W_{en}^{l+1} = F((W+W_{en})^{l+1})\f]
where Wnew_tot = \f$(W+W_{en})^{l+1}\f$ is the already-updated total water amount.
In this implementation the function \f$F\f$ is piecewise linear,
\f$F\le (\phi/(\rho_w g)) P_{new}\f$, and gives
\f$W_{en}^{l+1} \le 0.5 (W+W_{en})^{l+1}\f$.  Once this function is computed,
the new value \f$W^{l+1}\f$ is also established in a mass conserving manner.
 */
PetscErrorCode PISMDistributedHydrology::update_englacial_storage(
                                    IceModelVec2S &myPnew, IceModelVec2S &Wnew_tot) {
  PetscErrorCode ierr;
  const PetscReal rg = config.get("fresh_water_density") * config.get("standard_gravity"),
                  porosity = config.get("hydrology_englacial_porosity"), // the true porosity
                  CCpor = porosity / rg;
  PetscReal Wen_max;
  ierr = myPnew.begin_access(); CHKERRQ(ierr);
  ierr = Wen.begin_access(); CHKERRQ(ierr);
  ierr = Wnew_tot.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      // in next line: Wen_new satisfies scaled version of bounds 0 <= P <= P_o
      Wen_max = CCpor * myPnew(i,j);
      if (Wnew_tot(i,j) > Wen_max) {  // if there is enough water to charge englacial
                                      //   system and still have subglacial water present
        if (Wnew_tot(i,j) > 2.0 * Wen_max) {  // if there is plenty of water
          Wen(i,j) = Wen_max;
        } else {
          Wen(i,j) = Wnew_tot(i,j) - Wen_max;
        }
        Wnew_tot(i,j) -= Wen(i,j); // so the same amount moved out of subglacial
      } else {
        Wen(i,j) = 0.0;
      }
      // now the meaning can revert:  Wnew(i,j) := Wnew_tot(i,j)
    }
  }
  ierr = myPnew.end_access(); CHKERRQ(ierr);
  ierr = Wen.end_access(); CHKERRQ(ierr);
  ierr = Wnew_tot.end_access(); CHKERRQ(ierr);
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
  ierr = W.update_ghosts(); CHKERRQ(ierr);
  ierr = P.update_ghosts(); CHKERRQ(ierr);

  // from current ice geometry/velocity variables, initialize Po and cbase
  ierr = update_cbase(cbase); CHKERRQ(ierr);

  const PetscReal
            rg = config.get("fresh_water_density") * config.get("standard_gravity"),
            nglen = config.get("Glen_exponent"),
            Aglen = config.get("ice_softness"),
            c1 = config.get("hydrology_cavitation_opening_coefficient"),
            c2 = config.get("hydrology_creep_closure_coefficient"),
            Wr = config.get("hydrology_roughness_scale"),
            Y0 = config.get("hydrology_lower_bound_creep_regularization"),
            phisum = config.get("hydrology_englacial_porosity")
                       + config.get("hydrology_regularizing_porosity");

  if (phisum <= 0.0) {
    PetscPrintf(grid.com,
        "PISM ERROR:  phisum = englacial_porosity + regularizing_porosity <= 0 ... ENDING\n");
    PISMEnd();
  }

  const PetscReal  omegax = 1.0 / (grid.dx * grid.dx),
                   omegay = 1.0 / (grid.dy * grid.dy);

  PetscReal ht = t, hdt, // hydrology model time and time step
            maxKW, maxV, maxD;
  PetscReal icefreelost = 0, oceanlost = 0, negativegain = 0, nullstriplost = 0,
            delta_icefree, delta_ocean, delta_neggain, delta_nullstrip;

  PetscReal PtoCFLratio,  // for reporting ratio of dtCFL to dtDIFFP
            cumratio = 0.0;
  PetscInt hydrocount = 0; // count hydrology time steps

  while (ht < t + dt) {
    hydrocount++;

    ierr = check_W_nonnegative(); CHKERRQ(ierr);
    ierr = check_Wen_nonnegative(); CHKERRQ(ierr);

    // note that ice dynamics can change overburden pressure, so we can only check P
    //   bounds if thk has not changed; if thk could have just changed, such as in the
    //   first time through the current loop, we enforce them
    ierr = check_P_bounds((hydrocount == 1)); CHKERRQ(ierr);

    ierr = subglacial_hydraulic_potential(psi); CHKERRQ(ierr);
    ierr = psi.update_ghosts(); CHKERRQ(ierr);

    ierr = water_thickness_staggered(Wstag); CHKERRQ(ierr);
    ierr = Wstag.update_ghosts(); CHKERRQ(ierr);

    ierr = conductivity_staggered(Kstag,maxKW); CHKERRQ(ierr);
    ierr = Kstag.update_ghosts(); CHKERRQ(ierr);

    ierr = velocity_staggered(V); CHKERRQ(ierr);

    // to get Qstag, W needs valid ghosts
    ierr = advective_fluxes(Qstag); CHKERRQ(ierr);
    ierr = Qstag.update_ghosts(); CHKERRQ(ierr);

    ierr = adaptive_for_WandP_evolution(ht, t+dt, maxKW, hdt, maxV, maxD, PtoCFLratio); CHKERRQ(ierr);
    cumratio += PtoCFLratio;

    if ((inputtobed != NULL) || (hydrocount==1)) {
      ierr = get_input_rate(ht,hdt,total_input); CHKERRQ(ierr);
    }

    // update Pnew from time step
    const PetscReal  CC = (rg * hdt) / phisum;
    PetscReal  Open, Close, divflux,
               dpsie, dpsiw, dpsin, dpsis;
    ierr = overburden_pressure(Pover); CHKERRQ(ierr);

    MaskQuery M(*mask);

    ierr = P.begin_access(); CHKERRQ(ierr);
    ierr = W.begin_access(); CHKERRQ(ierr);
    ierr = cbase.begin_access(); CHKERRQ(ierr);
    ierr = psi.begin_access(); CHKERRQ(ierr);
    ierr = Wstag.begin_access(); CHKERRQ(ierr);
    ierr = Kstag.begin_access(); CHKERRQ(ierr);
    ierr = total_input.begin_access(); CHKERRQ(ierr);
    ierr = mask->begin_access(); CHKERRQ(ierr);
    ierr = Pover.begin_access(); CHKERRQ(ierr);
    ierr = Pnew.begin_access(); CHKERRQ(ierr);
    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
        if (M.ice_free_land(i,j))
          Pnew(i,j) = 0.0;
        else if (M.ocean(i,j))
          Pnew(i,j) = Pover(i,j);
        else if (W(i,j) <= 0.0)
          Pnew(i,j) = Pover(i,j);
        else {
          // opening and closure terms in pressure equation
          Open = PetscMax(0.0,c1 * cbase(i,j) * (Wr - W(i,j)));
          Close = c2 * Aglen * pow(Pover(i,j) - P(i,j),nglen) * (W(i,j) + Y0);
          // divergence of flux
          const bool knowne = (M.ice_free_land(i+1,j) || M.ocean(i+1,j)),
                     knownw = (M.ice_free_land(i-1,j) || M.ocean(i-1,j)),
                     knownn = (M.ice_free_land(i,j+1) || M.ocean(i,j+1)),
                     knowns = (M.ice_free_land(i,j-1) || M.ocean(i,j-1));
          dpsie = psi(i+1,j) - psi(i,j);
          dpsiw = psi(i,j)   - psi(i-1,j);
          dpsin = psi(i,j+1) - psi(i,j);
          dpsis = psi(i,j)   - psi(i,j-1);
          if (stripwidth > 0.0) {
            const bool nullij = (in_null_strip(i,j));
            if (nullij || in_null_strip(i+1,j))
              dpsie = 0.0;
            if (nullij || in_null_strip(i-1,j))
              dpsiw = 0.0;
            if (nullij || in_null_strip(i,j+1))
              dpsin = 0.0;
            if (nullij || in_null_strip(i,j-1))
              dpsis = 0.0;
          }
          divflux = 0.0;
          if (!knowne && !knownw) {
            const PetscReal We = Wstag(i,  j,0),
                            Ww = Wstag(i-1,j,0),
                            Ke = Kstag(i,  j,0),
                            Kw = Kstag(i-1,j,0);
            divflux += omegax * ( Ke * We * dpsie - Kw * Ww * dpsiw );
          }
          if (!knownn && !knowns) {
            const PetscReal Wn = Wstag(i,j  ,1),
                            Ws = Wstag(i,j-1,1);
            const PetscReal Kn = Kstag(i,j  ,1),
                            Ks = Kstag(i,j-1,1);
            divflux += omegay * ( Kn * Wn * dpsin - Ks * Ws * dpsis );
          }

          // candidate for pressure update
          Pnew(i,j) = P(i,j) + CC * ( divflux + Close - Open + total_input(i,j) );
          // projection to enforce  0 <= P <= P_o
          Pnew(i,j) = PetscMin(PetscMax(0.0, Pnew(i,j)), Pover(i,j));
        }
      }
    }
    ierr = P.end_access(); CHKERRQ(ierr);
    ierr = W.end_access(); CHKERRQ(ierr);
    ierr = cbase.end_access(); CHKERRQ(ierr);
    ierr = Pnew.end_access(); CHKERRQ(ierr);
    ierr = Pover.end_access(); CHKERRQ(ierr);
    ierr = psi.end_access(); CHKERRQ(ierr);
    ierr = total_input.end_access(); CHKERRQ(ierr);
    ierr = Wstag.end_access(); CHKERRQ(ierr);
    ierr = Kstag.end_access(); CHKERRQ(ierr);
    ierr = mask->end_access(); CHKERRQ(ierr);

    // update Wtotnew from W, Wstag, Qstag, total_input; the physics is
    // subglacial water movement:
    //    Wnew^{l+1} = W + (subglacial fluxes) + dt * total_input
    ierr = PISMRoutingHydrology::raw_update_W(hdt); CHKERRQ(ierr);

    // now include Wen = englacial into total water supply at new state
    //    Wnew_tot = (Wnew + Wen)^{l+1}
    ierr = Wnew.add(1.0,Wen); CHKERRQ(ierr);

    // now update Wen from knowledge of Pnew and Wnew_tot:
    //    Wen = C Pnew    if there is sufficient water, less otherwise
    // and then Wnew_tot -= Wen
    // and then revert meaning: Wnew_tot --> Wnew^{l+1}
    ierr = update_englacial_storage(Pnew, Wnew); CHKERRQ(ierr);

    ierr = Pnew.update_ghosts(P); CHKERRQ(ierr);

    ierr = boundary_mass_changes_with_null(Wnew,delta_icefree, delta_ocean,
                                 delta_neggain, delta_nullstrip); CHKERRQ(ierr);
    icefreelost  += delta_icefree;
    oceanlost    += delta_ocean;
    negativegain += delta_neggain;
    nullstriplost+= delta_nullstrip;

    ierr = Wnew.update_ghosts(W); CHKERRQ(ierr);

    ht += hdt;
  } // end of hydrology model time-stepping loop

  if (report_mass_accounting) {
    const PetscReal dtavyears = grid.convert(dt/hydrocount, "seconds", "years");
    ierr = verbPrintf(2, grid.com,
                      " 'distributed' hydrology summary:\n"
                      "     %d hydrology sub-steps with average dt = %.7f years = %.2f s\n"
                      "        (average of %.2f steps per CFL time; last max |V| = %.2e m s-1; last max D = %.2e m^2 s-1)\n"
                      "     ice free land lost = %.3e kg, ocean lost = %.3e kg\n"
                      "     negative bmelt gain = %.3e kg, null strip lost = %.3e kg\n",
                      hydrocount, dtavyears,
                      grid.convert(dtavyears, "seconds", "years"),
                      cumratio/hydrocount, maxV, maxD,
                      icefreelost, oceanlost, negativegain, nullstriplost); CHKERRQ(ierr);
  }
  return 0;
}

