// Copyright (C) 2012-2013 PISM Authors
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
#include "flowlaws.hh"
#include "ShallowStressBalance.hh"
#include "PISMStressBalance.hh"


/************************************/
/******** PISMLakesHydrology ********/
/************************************/

PISMLakesHydrology::PISMLakesHydrology(IceGrid &g, const NCConfigVariable &conf)
    : PISMHydrology(g, conf)
{
  stripwidth = config.get("hydrology_null_strip_width");
  if (allocate() != 0) {
    PetscPrintf(grid.com, "PISM ERROR: memory allocation failed in PISMLakesHydrology constructor.\n");
    PISMEnd();
  }
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
  ierr = Wstag.create(grid, "W_staggered", true, 1); CHKERRQ(ierr);
  ierr = Wstag.set_attrs("internal",
                     "cell face-centered (staggered) values of water layer thickness",
                     "m", ""); CHKERRQ(ierr);
  ierr = Wstag.set_attr("valid_min", 0.0); CHKERRQ(ierr);
  ierr = Qstag.create(grid, "advection_flux", true, 1); CHKERRQ(ierr);
  ierr = Qstag.set_attrs("internal",
                     "cell face-centered (staggered) components of advective subglacial water flux",
                     "m2 s-1", ""); CHKERRQ(ierr);
  ierr = Pwork.create(grid, "water_pressure_workspace", true, 1); CHKERRQ(ierr);
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
  ierr = Wnew.create(grid, "Wnew_internal", false); CHKERRQ(ierr);
  ierr = Wnew.set_attrs("internal",
                     "new thickness of subglacial water layer during update",
                     "m", ""); CHKERRQ(ierr);
  ierr = Wnew.set_attr("valid_min", 0.0); CHKERRQ(ierr);

  return 0;
}


PetscErrorCode PISMLakesHydrology::init(PISMVars &vars) {
  PetscErrorCode ierr;
  ierr = verbPrintf(2, grid.com,
    "* Initializing the subglacial-lakes-type subglacial hydrology model...\n"); CHKERRQ(ierr);
  // initialize water layer thickness from the context if present,
  //   otherwise from -i or -boot_file, otherwise with constant value
  bool i, bootstrap, stripset;
  ierr = PetscOptionsBegin(grid.com, "",
            "Options controlling the 'lakes' subglacial hydrology model", ""); CHKERRQ(ierr);
  {
    ierr = PISMOptionsIsSet("-i", "PISM input file", i); CHKERRQ(ierr);
    ierr = PISMOptionsIsSet("-boot_file", "PISM bootstrapping file",
                            bootstrap); CHKERRQ(ierr);
    ierr = PISMOptionsReal("-hydrology_null_strip",
                           "set the width, in km, of the strip around the edge of the computational domain in which hydrology is inactivated",
                           stripwidth,stripset); CHKERRQ(ierr);
    if (stripset) stripwidth *= 1.0e3;
  }
  ierr = PetscOptionsEnd(); CHKERRQ(ierr);
  // does not produce confusing messages in derived classes:
  ierr = init_actions(vars,i,bootstrap); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode PISMLakesHydrology::init_actions(PISMVars &vars, bool i_set, bool bootstrap_set) {
  PetscErrorCode ierr;

  ierr = PISMHydrology::init(vars); CHKERRQ(ierr);

  IceModelVec2S *W_input = dynamic_cast<IceModelVec2S*>(vars.get("bwat"));
  if (W_input != NULL) { // a variable called "bwat" is already in context
    ierr = W.copy_from(*W_input); CHKERRQ(ierr);
  } else if (i_set || bootstrap_set) {
    string filename;
    int start;
    ierr = find_pism_input(filename, bootstrap_set, start); CHKERRQ(ierr);
    if (i_set) {
      ierr = W.read(filename, start); CHKERRQ(ierr);
    } else {
      ierr = W.regrid(filename,
                      config.get("bootstrapping_bwat_value_no_var")); CHKERRQ(ierr);
    }
  } else {
    ierr = W.set(config.get("bootstrapping_bwat_value_no_var")); CHKERRQ(ierr);
  }

  // whether or not we could initialize from file, we could be asked to regrid from file
  ierr = regrid(W); CHKERRQ(ierr);

  // add bwat to the variables in the context if it is not already there
  if (vars.get("bwat") == NULL) {
    ierr = vars.add(W); CHKERRQ(ierr);
  }
  return 0;
}


void PISMLakesHydrology::add_vars_to_output(string /*keyword*/, map<string,NCSpatialVariable> &result) {
  result["bwat"] = W.get_metadata();
}


PetscErrorCode PISMLakesHydrology::define_variables(set<string> vars, const PIO &nc,
                                                 PISM_IO_Type nctype) {
  PetscErrorCode ierr;
  if (set_contains(vars, "bwat")) {
    ierr = W.define(nc, nctype); CHKERRQ(ierr);
  }
  return 0;
}


PetscErrorCode PISMLakesHydrology::write_variables(set<string> vars, const PIO &nc) {
  PetscErrorCode ierr;
  if (set_contains(vars, "bwat")) {
    ierr = W.write(nc); CHKERRQ(ierr);
  }
  return 0;
}


void PISMLakesHydrology::get_diagnostics(map<string, PISMDiagnostic*> &dict) {
  PISMHydrology::get_diagnostics(dict);
  dict["bwatvel"] = new PISMLakesHydrology_bwatvel(this, grid, *variables);
}


//! Check W >= 0 and fails with message if not satisfied.
PetscErrorCode PISMLakesHydrology::check_W_nonnegative() {
  PetscErrorCode ierr;
  ierr = W.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (W(i,j) < 0.0) {
        PetscPrintf(grid.com,
           "PISM ERROR: disallowed negative subglacial water layer thickness (bwat)\n"
           "    W(i,j) = %.6f m at (i,j)=(%d,%d)\n"
           "ENDING ... \n\n", W(i,j),i,j);
        PISMEnd();
      }
    }
  }
  ierr = W.end_access(); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode PISMLakesHydrology::boundary_mass_changes_with_null(
            IceModelVec2S &myWnew,
            PetscReal &icefreelost, PetscReal &oceanlost,
            PetscReal &negativegain, PetscReal &nullstriplost) {

  PetscErrorCode ierr;

  ierr = PISMHydrology::boundary_mass_changes(myWnew,
            icefreelost,oceanlost,negativegain); CHKERRQ(ierr);
  if (stripwidth <= 0.0)  return 0;

  PetscReal fresh_water_density = config.get("fresh_water_density");
  PetscReal my_nullstriplost = 0.0;

  ierr = myWnew.begin_access(); CHKERRQ(ierr);
  ierr = cellarea->begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      const PetscReal dmassdz = (*cellarea)(i,j) * fresh_water_density; // kg m-1
      if (in_null_strip(i,j)) {
        my_nullstriplost += Wnew(i,j) * dmassdz;
        Wnew(i,j) = 0.0;
      }
    }
  }
  ierr = myWnew.end_access(); CHKERRQ(ierr);
  ierr = cellarea->end_access(); CHKERRQ(ierr);

  ierr = PISMGlobalSum(&my_nullstriplost, &nullstriplost, grid.com); CHKERRQ(ierr);
  ierr = verbPrintf(4, grid.com,
    "     null strip loss = %.3e kg\n",nullstriplost); CHKERRQ(ierr);
  return 0;
}


//! Copies the W variable, the modeled water layer thickness.
PetscErrorCode PISMLakesHydrology::subglacial_water_thickness(IceModelVec2S &result) {
  PetscErrorCode ierr = W.copy_to(result); CHKERRQ(ierr);
  return 0;
}


//! Computes pressure diagnostically as fixed fraction of overburden.
/*!
Here
  \f[ P = \lambda P_o = \lambda (\rho_i g H) \f]
where \f$\lambda\f$=till_pw_fraction and \f$P_o\f$ is the overburden pressure.
 */
PetscErrorCode PISMLakesHydrology::subglacial_water_pressure(IceModelVec2S &result) {
  PetscErrorCode ierr;
  // ierr = verbPrintf(1,grid.com,"in PISMLakesHydrology::water_pressure()\n"); CHKERRQ(ierr);
  ierr = overburden_pressure(result); CHKERRQ(ierr);
  ierr = result.scale(config.get("till_pw_fraction")); CHKERRQ(ierr);
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


//! Get the advection velocity V at the center of cell edges.
/*!
Computes the advection velocity \f$\mathbf{V}=\mathbf{V}(\nabla P,\nabla b)\f$
on the staggered (face-centered) grid.  If V = (alpha,beta) in components
then we have <code> result(i,j,0) = alpha(i+1/2,j) </code> and
<code> result(i,j,1) = beta(i,j+1/2) </code>

The advection velocity is given by the formula
  \f[ \mathbf{V} = - k W^{\alpha-1} \nabla \left(\frac{P}{\rho_w g} + b\right) \f]
where \f$\mathbf{V}\f$ is the lateral water velocity, \f$P\f$ is the water
pressure, and \f$b\f$ is the bedrock elevation.

If the corresponding staggered grid value of the water thickness is zero then
that component of V is set to zero.  This does not change the flux value (which
would be zero anyway) but it does provide the correct max velocity in the
CFL calculation.  We assume Wstag is up-to-date with up-to-date ghosts.

Calls water_pressure() method to get water pressure.
 */
PetscErrorCode PISMLakesHydrology::velocity_staggered(IceModelVec2Stag &result) {
  PetscErrorCode ierr;
  const PetscReal
    g    = config.get("standard_gravity"),
    rhow = config.get("fresh_water_density"),
    k    = config.get("hydrology_hydraulic_conductivity"),
    alpha = config.get("hydrology_thickness_power_in_flux");
  PetscReal dbdx, dbdy, dPdx, dPdy;

  ierr = subglacial_water_pressure(Pwork); CHKERRQ(ierr);  // yes, it updates ghosts

  ierr = Pwork.begin_access(); CHKERRQ(ierr);
  ierr = Wstag.begin_access(); CHKERRQ(ierr);
  ierr = bed->begin_access(); CHKERRQ(ierr);
  ierr = result.begin_access(); CHKERRQ(ierr);
  for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
      if (Wstag(i,j,0) > 0.0) {
        dPdx = (Pwork(i+1,j  ) - Pwork(i,j)) / grid.dx;
        dbdx = ((*bed)(i+1,j  ) - (*bed)(i,j)) / grid.dx;
        const PetscReal We = Wstag(i,j,0);
        result(i,j,0) = - k * pow(We,alpha-1) * (dPdx / (rhow * g) + dbdx);
      } else
        result(i,j,0) = 0.0;
      if (Wstag(i,j,1) > 0.0) {
        dPdy = (Pwork(i  ,j+1) - Pwork(i,j)) / grid.dy;
        dbdy = ((*bed)(i  ,j+1) - (*bed)(i,j)) / grid.dy;
        const PetscReal Wn = Wstag(i,j,1);
        result(i,j,1) = - k * pow(Wn,alpha-1) * (dPdy / (rhow * g) + dbdy);
      } else
        result(i,j,1) = 0.0;
      if (in_null_strip(i,j) || in_null_strip(i+1,j))
        result(i,j,0) = 0.0;
      if (in_null_strip(i,j) || in_null_strip(i,j+1))
        result(i,j,1) = 0.0;
    }
  }
  ierr = Pwork.end_access(); CHKERRQ(ierr);
  ierr = Wstag.end_access(); CHKERRQ(ierr);
  ierr = bed->end_access(); CHKERRQ(ierr);
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
                  PetscReal &maxV_result, PetscReal &dtCFL_result, PetscReal &dtDIFFW_result) {
  PetscErrorCode ierr;
  const PetscReal k     = config.get("hydrology_hydraulic_conductivity"),
                  alpha = config.get("hydrology_thickness_power_in_flux"),
                  dtmax = config.get("hydrology_maximum_time_step_years") * secpera;
  PetscReal maxW, tmp[2];
  // Matlab: dtCFL = 0.5 / (max(max(abs(alphV)))/dx + max(max(abs(betaV)))/dy);
  ierr = V.absmaxcomponents(tmp); CHKERRQ(ierr); // V could be zero if P is constant and bed is flat
  maxV_result = sqrt(tmp[0]*tmp[0] + tmp[1]*tmp[1]);
  // regularize with eps = (1 cm/a) / (100km) :
  dtCFL_result = 0.5 / (tmp[0]/grid.dx + tmp[1]/grid.dy); // is regularization needed?
  // Matlab: maxW = max(max(max(Wea)),max(max(Wno))) + 0.001;
  ierr = Wstag.absmaxcomponents(tmp); CHKERRQ(ierr);
  maxW = PetscMax(tmp[0],tmp[1]) + 0.001;
  // Matlab: dtDIFFW = 0.25 / (k * maxW^alpha * (1/dx^2 + 1/dy^2));
  dtDIFFW_result = 1.0/(grid.dx*grid.dx) + 1.0/(grid.dy*grid.dy);
  dtDIFFW_result = 0.25 / (k * pow(maxW,alpha) * dtDIFFW_result);
  // dt = min([te-t dtmax dtCFL dtDIFFW]);
  dt_result = PetscMin(t_end - t_current, dtmax);
  dt_result = PetscMin(dt_result, dtCFL_result);
  dt_result = PetscMin(dt_result, dtDIFFW_result);
  return 0;
}


/*!
Call this version if you don't need the dt associated to the diffusion term.
 */
PetscErrorCode PISMLakesHydrology::adaptive_for_W_evolution(
                  PetscReal t_current, PetscReal t_end, PetscReal &dt_result) {
  PetscReal discard1, discard2, discard3;
  PetscErrorCode ierr = adaptive_for_W_evolution(
                             t_current, t_end, dt_result,
                             discard1, discard2, discard3); CHKERRQ(ierr);
  return 0;
}


//! The computation of Wnew, called by update().
PetscErrorCode PISMLakesHydrology::raw_update_W(PetscReal hdt) {
    PetscErrorCode ierr;
    const PetscReal
      k    = config.get("hydrology_hydraulic_conductivity"),
      alpha = config.get("hydrology_thickness_power_in_flux"),
      wux  = 1.0 / (grid.dx * grid.dx),
      wuy  = 1.0 / (grid.dy * grid.dy);
    PetscReal divadflux, diffW;

    ierr = W.begin_access(); CHKERRQ(ierr);
    ierr = Wstag.begin_access(); CHKERRQ(ierr);
    ierr = Qstag.begin_access(); CHKERRQ(ierr);
    ierr = total_input.begin_access(); CHKERRQ(ierr);
    ierr = Wnew.begin_access(); CHKERRQ(ierr);
    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
        divadflux =   (Qstag(i,j,0) - Qstag(i-1,j  ,0)) / grid.dx
                    + (Qstag(i,j,1) - Qstag(i,  j-1,1)) / grid.dy;
        const PetscReal  We = Wstag(i,  j,0),
                         Ww = Wstag(i-1,j,0),
                         Wn = Wstag(i,j  ,1),
                         Ws = Wstag(i,j-1,1);
        diffW =   wux * (  k * pow(We,alpha) * (W(i+1,j) - W(i,j))
                         - k * pow(Ww,alpha) * (W(i,j) - W(i-1,j)) )
                + wuy * (  k * pow(Wn,alpha) * (W(i,j+1) - W(i,j))
                         - k * pow(Ws,alpha) * (W(i,j) - W(i,j-1)) );
        Wnew(i,j) = W(i,j) + hdt * (- divadflux + diffW + total_input(i,j));
      }
    }
    ierr = W.end_access(); CHKERRQ(ierr);
    ierr = Wstag.end_access(); CHKERRQ(ierr);
    ierr = Qstag.end_access(); CHKERRQ(ierr);
    ierr = total_input.end_access(); CHKERRQ(ierr);
    ierr = Wnew.end_access(); CHKERRQ(ierr);

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
  ierr = W.update_ghosts(); CHKERRQ(ierr);

  MaskQuery M(*mask);
  PetscReal ht = t, hdt; // hydrology model time and time step
  PetscReal icefreelost = 0, oceanlost = 0, negativegain = 0, nullstriplost = 0,
            delta_icefree, delta_ocean, delta_neggain, delta_nullstrip;
  PetscInt hydrocount = 0; // count hydrology time steps

  while (ht < t + dt) {
    hydrocount++;

    ierr = check_W_nonnegative(); CHKERRQ(ierr);

    ierr = water_thickness_staggered(Wstag); CHKERRQ(ierr);
    ierr = Wstag.update_ghosts(); CHKERRQ(ierr);

    ierr = velocity_staggered(V); CHKERRQ(ierr);

    // to get Qstag, W needs valid ghosts
    ierr = advective_fluxes(Qstag); CHKERRQ(ierr);
    ierr = Qstag.update_ghosts(); CHKERRQ(ierr);

    ierr = adaptive_for_W_evolution(ht, t+dt, hdt); CHKERRQ(ierr);

    if ((inputtobed != NULL) || (hydrocount==1)) {
      ierr = get_input_rate(ht,hdt,total_input); CHKERRQ(ierr);
    }

    // update Wnew (the actual step) from W, Wstag, Qstag, total_input
    ierr = raw_update_W(hdt); CHKERRQ(ierr);

    ierr = boundary_mass_changes_with_null(Wnew,delta_icefree, delta_ocean,
                                 delta_neggain, delta_nullstrip); CHKERRQ(ierr);
    icefreelost  += delta_icefree;
    oceanlost    += delta_ocean;
    negativegain += delta_neggain;
    nullstriplost+= delta_nullstrip;

    // transfer Wnew into W
    ierr = Wnew.update_ghosts(W); CHKERRQ(ierr);

    ht += hdt;
  } // end of hydrology model time-stepping loop

  if (report_mass_accounting) {
    ierr = verbPrintf(2, grid.com,
      " 'lakes' hydrology summary:\n"
      "     %d hydrology sub-steps with average dt = %.6f years\n"
      "     ice free land lost = %.3e kg, ocean lost = %.3e kg\n"
      "     negative bmelt gain = %.3e kg, null strip lost = %.3e kg\n",
      hydrocount, (dt/hydrocount)/secpera,
      icefreelost, oceanlost, negativegain, nullstriplost); CHKERRQ(ierr);
  }
  return 0;
}


PISMLakesHydrology_bwatvel::PISMLakesHydrology_bwatvel(PISMLakesHydrology *m, IceGrid &g, PISMVars &my_vars)
    : PISMDiag<PISMLakesHydrology>(m, g, my_vars) {

  // set metadata:
  dof = 2;
  vars.resize(2);
  vars[0].init_2d("bwatvel[0]", grid);
  vars[1].init_2d("bwatvel[1]", grid);

  set_attrs("velocity of water in subglacial layer, i-offset", "",
            "m s-1", "m year-1", 0);
  set_attrs("velocity of water in subglacial layer, j-offset", "",
            "m s-1", "m year-1", 1);
}


PetscErrorCode PISMLakesHydrology_bwatvel::compute(IceModelVec* &output) {
  PetscErrorCode ierr;

  IceModelVec2Stag *result = new IceModelVec2Stag;
  ierr = result->create(grid, "bwatvel", true); CHKERRQ(ierr);
  ierr = result->set_metadata(vars[0], 0); CHKERRQ(ierr);
  ierr = result->set_metadata(vars[1], 1); CHKERRQ(ierr);
  result->write_in_glaciological_units = true;

  ierr = model->velocity_staggered(*result); CHKERRQ(ierr);

  output = result;
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

  // additional variables beyond PISMLakesHydrology::allocate()
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

  ierr = PISMLakesHydrology::init_actions(vars,i_set,bootstrap_set); CHKERRQ(ierr);

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


void PISMDistributedHydrology::add_vars_to_output(string /*keyword*/, map<string,NCSpatialVariable> &result) {
  result["bwat"] = W.get_metadata();
  result["enwat"] = Wen.get_metadata();
  result["bwp"]  = P.get_metadata();
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
  dict["bwatvel"] = new PISMLakesHydrology_bwatvel(this, grid, *variables);
  dict["bwprel"] = new PISMHydrology_bwprel(this, grid, *variables);
  dict["effbwp"] = new PISMHydrology_effbwp(this, grid, *variables);
  dict["hydroinput"] = new PISMHydrology_hydroinput(this, grid, *variables);
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

  ierr = overburden_pressure(Pwork); CHKERRQ(ierr);

  ierr = P.begin_access(); CHKERRQ(ierr);
  ierr = Pwork.begin_access(); CHKERRQ(ierr);
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
        P(i,j) = PetscMin(P(i,j), Pwork(i,j));
      } else if (P(i,j) > Pwork(i,j) + 0.001) {
        PetscPrintf(grid.com,
           "PISM ERROR: subglacial water pressure P = %.16f Pa exceeds\n"
           "    overburden pressure Po = %.16f Pa at (i,j)=(%d,%d)\n"
           "ENDING ... \n\n", P(i,j),Pwork(i,j),i,j);
        PISMEnd();
      }
    }
  }
  ierr = P.end_access(); CHKERRQ(ierr);
  ierr = Pwork.end_access(); CHKERRQ(ierr);
  return 0;
}


//! Get the hydraulic potential from bedrock topography and current state variables.
/*!
Computes \f$\psi = P + \rho_w g (b + W)\f$ except where floating, where \f$\psi = P_o\f$.

Calls water_pressure() method to get water pressure.
 */
PetscErrorCode PISMDistributedHydrology::hydraulic_potential(IceModelVec2S &result) {
  PetscErrorCode ierr;
  PetscReal fresh_water_density = config.get("fresh_water_density"),
            standard_gravity    = config.get("standard_gravity");
  MaskQuery M(*mask);
  ierr = overburden_pressure(Pwork); CHKERRQ(ierr);
  ierr = W.begin_access(); CHKERRQ(ierr);
  ierr = P.begin_access(); CHKERRQ(ierr);
  ierr = Pwork.begin_access(); CHKERRQ(ierr);
  ierr = bed->begin_access(); CHKERRQ(ierr);
  ierr = mask->begin_access(); CHKERRQ(ierr);
  ierr = result.begin_access(); CHKERRQ(ierr);
  for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
      if (M.ocean(i,j))
        result(i,j) = Pwork(i,j);
      else
        result(i,j) = P(i,j) + fresh_water_density * standard_gravity * ((*bed)(i,j) + W(i,j));
    }
  }
  ierr = W.end_access(); CHKERRQ(ierr);
  ierr = P.end_access(); CHKERRQ(ierr);
  ierr = Pwork.end_access(); CHKERRQ(ierr);
  ierr = bed->end_access(); CHKERRQ(ierr);
  ierr = mask->end_access(); CHKERRQ(ierr);
  ierr = result.end_access(); CHKERRQ(ierr);
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
  ierr = overburden_pressure(Pwork); CHKERRQ(ierr);
  ierr = W.begin_access(); CHKERRQ(ierr);
  ierr = Pwork.begin_access(); CHKERRQ(ierr);
  ierr = cbase.begin_access(); CHKERRQ(ierr);
  ierr = result.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      sb     = pow(CC * cbase(i,j),powglen);
      Wratio = PetscMax(0.0,Wr - W(i,j)) / (W(i,j) + Y0);
      // in cases where steady state is actually possible this will
      //   come out positive, but otherwise we should get underpressure P=0,
      //   and that is what it yields
      result(i,j) = PetscMax( 0.0,Pwork(i,j) - sb * pow(Wratio,powglen) );
    }
  }
  ierr = W.end_access(); CHKERRQ(ierr);
  ierr = Pwork.end_access(); CHKERRQ(ierr);
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
                  PetscReal t_current, PetscReal t_end,
                  PetscReal &dt_result, PetscReal &PtoCFLratio) {
  PetscErrorCode ierr;
  PetscReal maxV, dtCFL, dtDIFFW, dtDIFFP;

  ierr = adaptive_for_W_evolution(t_current,t_end,
              dt_result,maxV,dtCFL,dtDIFFW); CHKERRQ(ierr);

  const PetscReal rhoratio = config.get("fresh_water_density") / config.get("ice_density"),
                  phisum   = config.get("hydrology_englacial_porosity") + config.get("hydrology_regularizing_porosity");
  dtDIFFP = rhoratio * phisum * dtDIFFW;

  // dt = min([te-t dtmax dtCFL dtDIFFW dtDIFFP]);
  dt_result = PetscMin(dt_result, dtDIFFP);

  if (dtDIFFP > 0.0)
    PtoCFLratio = PetscMax(1.0, dtCFL / dtDIFFP);
  else
    PtoCFLratio = 1.0;

  ierr = verbPrintf(3,grid.com,
            "   [%.5e  %.7f  %.6f  %.9f  -->  dt = %.9f (a)  at  t = %.6f (a)]\n",
            maxV*secpera, dtCFL/secpera, dtDIFFW/secpera, dtDIFFP/secpera,
            dt_result/secpera, t_current/secpera); CHKERRQ(ierr);
  return 0;
}


//! Update conserved model state variables W,Wen by computing the englacial-to-subglacial connection using the change in subglacial pressure.
PetscErrorCode PISMDistributedHydrology::update_amounts_englacial_connection(IceModelVec2S &myPnew) {
  PetscErrorCode ierr;
  const PetscReal rhow = config.get("fresh_water_density"),
                  g = config.get("standard_gravity"),
                  porosity = config.get("hydrology_englacial_porosity"), // the true porosity
                  CCpor = porosity / (rhow * g);
  PetscReal Wen_new;
  ierr = myPnew.begin_access(); CHKERRQ(ierr);
  ierr = P.begin_access(); CHKERRQ(ierr);
  ierr = Wen.begin_access(); CHKERRQ(ierr);
  ierr = Wnew.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      // in next line: (i)  Wen satisfies scaled version of bounds 0 <= P <= P_o
      //               (ii) Wen comes from connection to water at base so Wen=0 if W=0
      Wen_new = (Wnew(i,j) > 0) ? CCpor * Pnew(i,j) : 0.0;
      const PetscReal deltaWen = Wen_new - Wen(i,j); // this much water moved into englacial
      Wen(i,j) = Wen_new;
      Wnew(i,j) -= deltaWen; // so the same amount moved out of subglacial
    }
  }
  ierr = myPnew.end_access(); CHKERRQ(ierr);
  ierr = P.end_access(); CHKERRQ(ierr);
  ierr = Wen.end_access(); CHKERRQ(ierr);
  ierr = Wnew.end_access(); CHKERRQ(ierr);
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
            k    = config.get("hydrology_hydraulic_conductivity"),
            alpha = config.get("hydrology_thickness_power_in_flux"),
            rhow = config.get("fresh_water_density"),
            g = config.get("standard_gravity"),
            nglen = config.get("Glen_exponent"),
            Aglen = config.get("ice_softness"),
            c1 = config.get("hydrology_cavitation_opening_coefficient"),
            c2 = config.get("hydrology_creep_closure_coefficient"),
            Wr = config.get("hydrology_roughness_scale"),
            Y0 = config.get("hydrology_lower_bound_creep_regularization"),
            phisum = config.get("hydrology_englacial_porosity") + config.get("hydrology_regularizing_porosity");

  if (phisum <= 0.0) {
    PetscPrintf(grid.com,
        "PISM ERROR:  phisum = englacial_porosity + regularizing_porosity <= 0 ... ENDING\n");
    PISMEnd();
  }

  const PetscReal  omegax = 1.0 / (rhow * g * grid.dx * grid.dx),
                   omegay = 1.0 / (rhow * g * grid.dy * grid.dy);

  PetscReal ht = t, hdt; // hydrology model time and time step

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

    ierr = hydraulic_potential(psi); CHKERRQ(ierr);
    ierr = psi.update_ghosts(); CHKERRQ(ierr);

    ierr = water_thickness_staggered(Wstag); CHKERRQ(ierr);
    ierr = Wstag.update_ghosts(); CHKERRQ(ierr);

    ierr = velocity_staggered(V); CHKERRQ(ierr);

    // to get Qstag, W needs valid ghosts
    ierr = advective_fluxes(Qstag); CHKERRQ(ierr);
    ierr = Qstag.update_ghosts(); CHKERRQ(ierr);

    ierr = adaptive_for_WandP_evolution(ht, t+dt, hdt, PtoCFLratio); CHKERRQ(ierr);
    cumratio += PtoCFLratio;

    if ((inputtobed != NULL) || (hydrocount==1)) {
      ierr = get_input_rate(ht,hdt,total_input); CHKERRQ(ierr);
    }

    // update Pnew from time step
    const PetscReal  CC = (rhow * g * hdt) / phisum;
    PetscReal  Open, Close, divflux,
               dpsie, dpsiw, dpsin, dpsis;
    ierr = overburden_pressure(Pwork); CHKERRQ(ierr);

    MaskQuery M(*mask);

    ierr = P.begin_access(); CHKERRQ(ierr);
    ierr = W.begin_access(); CHKERRQ(ierr);
    ierr = cbase.begin_access(); CHKERRQ(ierr);
    ierr = psi.begin_access(); CHKERRQ(ierr);
    ierr = Wstag.begin_access(); CHKERRQ(ierr);
    ierr = total_input.begin_access(); CHKERRQ(ierr);
    ierr = mask->begin_access(); CHKERRQ(ierr);
    ierr = Pwork.begin_access(); CHKERRQ(ierr);
    ierr = Pnew.begin_access(); CHKERRQ(ierr);
    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
        if (M.ice_free_land(i,j))
          Pnew(i,j) = 0.0;
        else if (M.ocean(i,j))
          Pnew(i,j) = Pwork(i,j);
        else if (W(i,j) <= 0.0)
          Pnew(i,j) = Pwork(i,j);
        else {
          // opening and closure terms in pressure equation
          Open = PetscMax(0.0,c1 * cbase(i,j) * (Wr - W(i,j)));
          Close = c2 * Aglen * pow(Pwork(i,j) - P(i,j),nglen) * (W(i,j) + Y0);
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
                            Ww = Wstag(i-1,j,0);
            divflux += omegax * k * ( pow(We,alpha) * dpsie - pow(Ww,alpha) * dpsiw );
          }
          if (!knownn && !knowns) {
            const PetscReal Wn = Wstag(i,j  ,1),
                            Ws = Wstag(i,j-1,1);
            divflux += omegay * k * ( pow(Wn,alpha) * dpsin - pow(Ws,alpha) * dpsis );
          }

          // candidate for pressure update
          Pnew(i,j) = P(i,j) + CC * ( divflux + Close - Open + total_input(i,j) );
          // projection to enforce  0 <= P <= P_o
          Pnew(i,j) = PetscMin(PetscMax(0.0, Pnew(i,j)), Pwork(i,j));
        }
      }
    }
    ierr = P.end_access(); CHKERRQ(ierr);
    ierr = W.end_access(); CHKERRQ(ierr);
    ierr = cbase.end_access(); CHKERRQ(ierr);
    ierr = Pnew.end_access(); CHKERRQ(ierr);
    ierr = Pwork.end_access(); CHKERRQ(ierr);
    ierr = psi.end_access(); CHKERRQ(ierr);
    ierr = total_input.end_access(); CHKERRQ(ierr);
    ierr = Wstag.end_access(); CHKERRQ(ierr);
    ierr = mask->end_access(); CHKERRQ(ierr);

    // update Wnew (the actual step) from W, Wstag, Qstag, total_input
    ierr = PISMLakesHydrology::raw_update_W(hdt); CHKERRQ(ierr);

    // update Wnew and Wen from  Delta P = Pnew - P
    ierr = update_amounts_englacial_connection(Pnew); CHKERRQ(ierr);

    // transfer Pnew into P; note Wstag, Qstag unaffected in Wnew update below
    ierr = Pnew.update_ghosts(P); CHKERRQ(ierr);

//FIXME: Wen?
    ierr = boundary_mass_changes_with_null(Wnew,delta_icefree, delta_ocean,
                                 delta_neggain, delta_nullstrip); CHKERRQ(ierr);
    icefreelost  += delta_icefree;
    oceanlost    += delta_ocean;
    negativegain += delta_neggain;
    nullstriplost+= delta_nullstrip;

    // transfer Wnew into W
    ierr = Wnew.update_ghosts(W); CHKERRQ(ierr);

    ht += hdt;
  } // end of hydrology model time-stepping loop

  if (report_mass_accounting) {
//FIXME: Wen?
    ierr = verbPrintf(2, grid.com,
      " 'distributed' hydrology summary:\n"
      "     %d hydrology sub-steps with average dt = %.6f years\n"
      "        (with average of %.2f sub-steps per CFL-forced step)\n"
      "     ice free land lost = %.3e kg, ocean lost = %.3e kg\n"
      "     negative bmelt gain = %.3e kg, null strip lost = %.3e kg\n",
      hydrocount, (dt/hydrocount)/secpera, cumratio/hydrocount,
      icefreelost, oceanlost, negativegain, nullstriplost); CHKERRQ(ierr);
  }
  return 0;
}

