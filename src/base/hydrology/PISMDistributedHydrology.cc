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
/******** PISMRoutingHydrology ********/
/************************************/

PISMRoutingHydrology::PISMRoutingHydrology(IceGrid &g, const NCConfigVariable &conf)
    : PISMHydrology(g, conf)
{
  stripwidth = config.get("hydrology_null_strip_width");
  if (allocate() != 0) {
    PetscPrintf(grid.com, "PISM ERROR: memory allocation failed in PISMRoutingHydrology constructor.\n");
    PISMEnd();
  }
}


PetscErrorCode PISMRoutingHydrology::allocate() {
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
  ierr = Kstag.create(grid, "K_staggered", true, 1); CHKERRQ(ierr);
  ierr = Kstag.set_attrs("internal",
                     "cell face-centered (staggered) values of nonlinear conductivity",
                     "", ""); CHKERRQ(ierr);
  ierr = Kstag.set_attr("valid_min", 0.0); CHKERRQ(ierr);
  ierr = Qstag.create(grid, "advection_flux", true, 1); CHKERRQ(ierr);
  ierr = Qstag.set_attrs("internal",
                     "cell face-centered (staggered) components of advective subglacial water flux",
                     "m2 s-1", ""); CHKERRQ(ierr);
  ierr = R.create(grid, "potential_workspace", true, 1); CHKERRQ(ierr); // box stencil used
  ierr = R.set_attrs("internal",
                      "work space for modeled subglacial water hydraulic potential",
                      "Pa", ""); CHKERRQ(ierr);

  // auxiliary variables which do not need ghosts
  ierr = Pover.create(grid, "overburden_pressure_internal", false); CHKERRQ(ierr);
  ierr = Pover.set_attrs("internal",
                     "overburden pressure",
                     "Pa", ""); CHKERRQ(ierr);
  ierr = Pover.set_attr("valid_min", 0.0); CHKERRQ(ierr);
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


PetscErrorCode PISMRoutingHydrology::init(PISMVars &vars) {
  PetscErrorCode ierr;
  ierr = verbPrintf(2, grid.com,
    "* Initializing the routing subglacial hydrology model ...\n"); CHKERRQ(ierr);
  // initialize water layer thickness from the context if present,
  //   otherwise from -i or -boot_file, otherwise with constant value
  bool i, bootstrap, stripset;
  ierr = PetscOptionsBegin(grid.com, "",
            "Options controlling the 'routing' subglacial hydrology model", ""); CHKERRQ(ierr);
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


PetscErrorCode PISMRoutingHydrology::init_actions(PISMVars &vars, bool i_set, bool bootstrap_set) {
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


void PISMRoutingHydrology::add_vars_to_output(string /*keyword*/, map<string,NCSpatialVariable> &result) {
  result["bwat"] = W.get_metadata();
}


PetscErrorCode PISMRoutingHydrology::define_variables(set<string> vars, const PIO &nc,
                                                 PISM_IO_Type nctype) {
  PetscErrorCode ierr;
  if (set_contains(vars, "bwat")) {
    ierr = W.define(nc, nctype); CHKERRQ(ierr);
  }
  return 0;
}


PetscErrorCode PISMRoutingHydrology::write_variables(set<string> vars, const PIO &nc) {
  PetscErrorCode ierr;
  if (set_contains(vars, "bwat")) {
    ierr = W.write(nc); CHKERRQ(ierr);
  }
  return 0;
}


void PISMRoutingHydrology::get_diagnostics(map<string, PISMDiagnostic*> &dict) {
  PISMHydrology::get_diagnostics(dict);
  dict["bwatvel"] = new PISMRoutingHydrology_bwatvel(this, grid, *variables);
}


//! Check W >= 0 and fails with message if not satisfied.
PetscErrorCode PISMRoutingHydrology::check_W_nonnegative() {
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


PetscErrorCode PISMRoutingHydrology::boundary_mass_changes_with_null(
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
PetscErrorCode PISMRoutingHydrology::subglacial_water_thickness(IceModelVec2S &result) {
  PetscErrorCode ierr = W.copy_to(result); CHKERRQ(ierr);
  return 0;
}


//! Computes pressure diagnostically as fixed fraction of overburden.
/*!
Here
  \f[ P = \lambda P_o = \lambda (\rho_i g H) \f]
where \f$\lambda\f$=till_pw_fraction and \f$P_o\f$ is the overburden pressure.
 */
PetscErrorCode PISMRoutingHydrology::subglacial_water_pressure(IceModelVec2S &result) {
  PetscErrorCode ierr;
  ierr = overburden_pressure(result); CHKERRQ(ierr);
  ierr = result.scale(config.get("till_pw_fraction")); CHKERRQ(ierr);  //FIXME issue #127
  return 0;
}


//! Average the regular grid water thickness to values at the center of cell edges.
PetscErrorCode PISMRoutingHydrology::water_thickness_staggered(IceModelVec2Stag &result) {
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


//! Compute the nonlinear conductivity at the center of cell edges.
/*!
Computes
    \f[ K = K(W,\nabla P, \nabla b) = k W^{\alpha-1} |\nabla(P+\rho_w g b)|^{\beta-2} \f]
on the staggered grid.  We denote \f$R = P+\rho_w g b\f$ internally.  The quantity
    \f[ \Pi = |\nabla(P+\rho_w g b)|^2 = |\nabla R|^2 \f]
is computed on a staggered grid by a [\ref Mahaffy] -like scheme.  This requires
\f$R\f$ to be defined on a box stencil of width 1.

Also computes the maximum over all staggered points of \f$ K W \f$.
 */
PetscErrorCode PISMRoutingHydrology::conductivity_staggered(
                       IceModelVec2Stag &result, PetscReal &maxKW) {
  PetscErrorCode ierr;
  const PetscReal
    k     = config.get("hydrology_hydraulic_conductivity"),
    alpha = config.get("hydrology_thickness_power_in_flux"),
    beta  = config.get("hydrology_potential_gradient_power_in_flux"),
    rg    = config.get("standard_gravity") * config.get("fresh_water_density");
  if (alpha < 1.0) {
    PetscPrintf(grid.com,
           "PISM ERROR: alpha = %f < 1 which is not allowed\n"
           "ENDING ... \n\n", alpha);
    PISMEnd();
  }

  if (beta == 2.0) {
    ierr = verbPrintf(4, grid.com,
      "    in PISMRoutingHydrology::conductivity_staggered(): "
      "beta == 2.0 exactly; simplifying calculation\n"); CHKERRQ(ierr);
  } else {
    // general case where beta is used; put norm of square gradient temporarily
    //   in result
    ierr = subglacial_water_pressure(R); CHKERRQ(ierr);  // yes, it updates ghosts
    ierr = R.add(rg, (*bed)); CHKERRQ(ierr); // R  <-- P + rhow g b
    ierr = R.update_ghosts(); CHKERRQ(ierr);

    PetscReal dRdx, dRdy;
    ierr = R.begin_access(); CHKERRQ(ierr);
    ierr = result.begin_access(); CHKERRQ(ierr);
    for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
      for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
        dRdx = ( R(i+1,j) - R(i,j) ) / grid.dx;
        dRdy = ( R(i+1,j+1) + R(i,j+1) - R(i+1,j-1) - R(i,j-1) ) / (4.0 * grid.dy);
        result(i,j,0) = dRdx * dRdx + dRdy * dRdy;
        dRdx = ( R(i+1,j+1) + R(i+1,j) - R(i-1,j+1) - R(i-1,j) ) / (4.0 * grid.dx);
        dRdy = ( R(i,j+1) - R(i,j) ) / grid.dy;
        result(i,j,1) = dRdx * dRdx + dRdy * dRdy;
      }
    }
    ierr = R.end_access(); CHKERRQ(ierr);
    ierr = result.end_access(); CHKERRQ(ierr);
  }

  PetscReal mymaxKW = 0.0;
  ierr = R.begin_access(); CHKERRQ(ierr);
  ierr = Wstag.begin_access(); CHKERRQ(ierr);
  ierr = result.begin_access(); CHKERRQ(ierr);
  for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
      for (PetscInt o = 0; o < 2; ++o) {
        if (beta == 2.0)
          result(i,j,o) = k * pow(Wstag(i,j,o),alpha-1.0);
        else {
          if ((result(i,j,o) <= 0.0) && (beta < 2.0)) {
            result(i,j,o) = 1000.0 * k;  // FIXME: ad hoc
          } else {
            result(i,j,o) = k * pow(Wstag(i,j,o),alpha-1.0)
                                * pow(result(i,j,o),(beta-2.0)/2.0);
          }
        }
        mymaxKW = PetscMax( mymaxKW, result(i,j,o) * Wstag(i,j,o) );
      }
    }
  }
  ierr = R.end_access(); CHKERRQ(ierr);
  ierr = Wstag.end_access(); CHKERRQ(ierr);
  ierr = result.end_access(); CHKERRQ(ierr);

  ierr = PISMGlobalMax(&mymaxKW, &maxKW, grid.com); CHKERRQ(ierr);

  return 0;
}


//! Get the advection velocity V at the center of cell edges.
/*!
Computes the advection velocity \f$\mathbf{V}\f$ on the staggered
(edge-centered) grid.  If V = (u,v) in components then we have
<code> result(i,j,0) = u(i+1/2,j) </code> and
<code> result(i,j,1) = v(i,j+1/2) </code>

The advection velocity is given by the formula
  \f[ \mathbf{V} = - K \left(\nabla P + \rho_w g \nabla b\right) \f]
where \f$\mathbf{V}\f$ is the water velocity, \f$P\f$ is the water
pressure, and \f$b\f$ is the bedrock elevation.

If the corresponding staggered grid value of the water thickness is zero then
that component of V is set to zero.  This does not change the flux value (which
would be zero anyway) but it does provide the correct max velocity in the
CFL calculation.  We assume Wstag and Kstag are up-to-date.  We assume P and b
have valid ghosts.

Calls subglacial_water_pressure() method to get water pressure.
 */
PetscErrorCode PISMRoutingHydrology::velocity_staggered(IceModelVec2Stag &result) {
  PetscErrorCode ierr;
  const PetscReal  rg = config.get("standard_gravity") * config.get("fresh_water_density");
  PetscReal dbdx, dbdy, dPdx, dPdy;

  ierr = subglacial_water_pressure(R); CHKERRQ(ierr);  // R=P; yes, it updates ghosts

  ierr = R.begin_access(); CHKERRQ(ierr);
  ierr = Wstag.begin_access(); CHKERRQ(ierr);
  ierr = Kstag.begin_access(); CHKERRQ(ierr);
  ierr = bed->begin_access(); CHKERRQ(ierr);
  ierr = result.begin_access(); CHKERRQ(ierr);
  for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
      if (Wstag(i,j,0) > 0.0) {
        dPdx = (R(i+1,j) - R(i,j)) / grid.dx;
        dbdx = ((*bed)(i+1,j) - (*bed)(i,j)) / grid.dx;
        result(i,j,0) = - Kstag(i,j,0) * (dPdx + rg * dbdx);
      } else
        result(i,j,0) = 0.0;
      if (Wstag(i,j,1) > 0.0) {
        dPdy = (R(i,j+1) - R(i,j)) / grid.dy;
        dbdy = ((*bed)(i,j+1) - (*bed)(i,j)) / grid.dy;
        result(i,j,1) = - Kstag(i,j,1) * (dPdy + rg * dbdy);
      } else
        result(i,j,1) = 0.0;
      if (in_null_strip(i,j) || in_null_strip(i+1,j))
        result(i,j,0) = 0.0;
      if (in_null_strip(i,j) || in_null_strip(i,j+1))
        result(i,j,1) = 0.0;
    }
  }
  ierr = R.end_access(); CHKERRQ(ierr);
  ierr = Wstag.end_access(); CHKERRQ(ierr);
  ierr = Kstag.end_access(); CHKERRQ(ierr);
  ierr = bed->end_access(); CHKERRQ(ierr);
  ierr = result.end_access(); CHKERRQ(ierr);
  return 0;
}


//! Compute Q = V W at edge-centers (staggered grid) by first-order upwinding.
/*!
The field W must have valid ghost values, but V does not need them.

FIXME:  This could be re-implemented using the Koren (1993) flux-limiter.
 */
PetscErrorCode PISMRoutingHydrology::advective_fluxes(IceModelVec2Stag &result) {
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
PetscErrorCode PISMRoutingHydrology::adaptive_for_W_evolution(
                  PetscReal t_current, PetscReal t_end, PetscReal maxKW,
                  PetscReal &dt_result,
                  PetscReal &maxV_result, PetscReal &maxD_result,
                  PetscReal &dtCFL_result, PetscReal &dtDIFFW_result) {
  PetscErrorCode ierr;
  const PetscReal
    dtmax = config.get("hydrology_maximum_time_step_years") * secpera,
    rg    = config.get("standard_gravity") * config.get("fresh_water_density");
  PetscReal tmp[2];
  ierr = V.absmaxcomponents(tmp); CHKERRQ(ierr); // V could be zero if P is constant and bed is flat
  maxV_result = sqrt(tmp[0]*tmp[0] + tmp[1]*tmp[1]);
  maxD_result = rg * maxKW;
  dtCFL_result = 0.5 / (tmp[0]/grid.dx + tmp[1]/grid.dy); // is regularization needed?
  dtDIFFW_result = 1.0/(grid.dx*grid.dx) + 1.0/(grid.dy*grid.dy);
  dtDIFFW_result = 0.25 / (maxD_result * dtDIFFW_result);
  // dt = min { te-t, dtmax, dtCFL, dtDIFFW }
  dt_result = PetscMin(t_end - t_current, dtmax);
  dt_result = PetscMin(dt_result, dtCFL_result);
  dt_result = PetscMin(dt_result, dtDIFFW_result);
  return 0;
}


//! The computation of Wnew, called by update().
PetscErrorCode PISMRoutingHydrology::raw_update_W(PetscReal hdt) {
    PetscErrorCode ierr;
    const PetscReal
      wux  = 1.0 / (grid.dx * grid.dx),
      wuy  = 1.0 / (grid.dy * grid.dy),
      rg   = config.get("standard_gravity") * config.get("fresh_water_density");
    PetscReal divadflux, diffW;

    ierr = W.begin_access(); CHKERRQ(ierr);
    ierr = Wstag.begin_access(); CHKERRQ(ierr);
    ierr = Kstag.begin_access(); CHKERRQ(ierr);
    ierr = Qstag.begin_access(); CHKERRQ(ierr);
    ierr = total_input.begin_access(); CHKERRQ(ierr);
    ierr = Wnew.begin_access(); CHKERRQ(ierr);
    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
        divadflux =   (Qstag(i,j,0) - Qstag(i-1,j  ,0)) / grid.dx
                    + (Qstag(i,j,1) - Qstag(i,  j-1,1)) / grid.dy;
        const PetscReal  De = rg * Kstag(i,  j,0) * Wstag(i,  j,0),
                         Dw = rg * Kstag(i-1,j,0) * Wstag(i-1,j,0),
                         Dn = rg * Kstag(i,j  ,1) * Wstag(i,j  ,1),
                         Ds = rg * Kstag(i,j-1,1) * Wstag(i,j-1,1);
        diffW =   wux * (  De * (W(i+1,j) - W(i,j)) - Dw * (W(i,j) - W(i-1,j)) )
                + wuy * (  Dn * (W(i,j+1) - W(i,j)) - Ds * (W(i,j) - W(i,j-1)) );
        Wnew(i,j) = W(i,j) + hdt * (- divadflux + diffW + total_input(i,j));
      }
    }
    ierr = W.end_access(); CHKERRQ(ierr);
    ierr = Wstag.end_access(); CHKERRQ(ierr);
    ierr = Kstag.end_access(); CHKERRQ(ierr);
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
PetscErrorCode PISMRoutingHydrology::update(PetscReal icet, PetscReal icedt) {
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
  PetscReal ht = t, hdt, // hydrology model time and time step
            maxKW, maxV, maxD, dtCFL, dtDIFFW;
  PetscReal icefreelost = 0, oceanlost = 0, negativegain = 0, nullstriplost = 0,
            delta_icefree, delta_ocean, delta_neggain, delta_nullstrip;
  PetscInt hydrocount = 0; // count hydrology time steps

  while (ht < t + dt) {
    hydrocount++;

    ierr = check_W_nonnegative(); CHKERRQ(ierr);

    ierr = water_thickness_staggered(Wstag); CHKERRQ(ierr);
    ierr = Wstag.update_ghosts(); CHKERRQ(ierr);

    ierr = conductivity_staggered(Kstag,maxKW); CHKERRQ(ierr);
    ierr = Kstag.update_ghosts(); CHKERRQ(ierr);

    ierr = velocity_staggered(V); CHKERRQ(ierr);

    // to get Qstag, W needs valid ghosts
    ierr = advective_fluxes(Qstag); CHKERRQ(ierr);
    ierr = Qstag.update_ghosts(); CHKERRQ(ierr);

    ierr = adaptive_for_W_evolution(ht, t+dt, maxKW,
                                    hdt, maxV, maxD, dtCFL, dtDIFFW); CHKERRQ(ierr);

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
      " 'routing' hydrology summary:\n"
      "     %d hydrology sub-steps with average dt = %.6f years = %.2f s\n"
      "        (last max |V| = %.2e m s-1; last max D = %.2e m^2 s-1)\n"
      "     ice free land lost = %.3e kg, ocean lost = %.3e kg\n"
      "     negative bmelt gain = %.3e kg, null strip lost = %.3e kg\n",
      hydrocount, (dt/hydrocount)/secpera, dt/hydrocount, maxV, maxD,
      icefreelost, oceanlost, negativegain, nullstriplost); CHKERRQ(ierr);
  }
  return 0;
}


PISMRoutingHydrology_bwatvel::PISMRoutingHydrology_bwatvel(PISMRoutingHydrology *m, IceGrid &g, PISMVars &my_vars)
    : PISMDiag<PISMRoutingHydrology>(m, g, my_vars) {

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


PetscErrorCode PISMRoutingHydrology_bwatvel::compute(IceModelVec* &output) {
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

  ierr = PISMRoutingHydrology::init_actions(vars,i_set,bootstrap_set); CHKERRQ(ierr);

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
  dict["bwatvel"] = new PISMRoutingHydrology_bwatvel(this, grid, *variables);
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


//! Get the hydraulic potential from bedrock topography and current state variables.
/*!
Computes \f$\psi = P + \rho_w g (b + W)\f$ except where floating, where \f$\psi = P_o\f$.

Calls water_pressure() method to get water pressure.
 */
PetscErrorCode PISMDistributedHydrology::hydraulic_potential(IceModelVec2S &result) {
  PetscErrorCode ierr;
  const PetscReal
    rg = config.get("fresh_water_density") * config.get("standard_gravity");
  MaskQuery M(*mask);
  ierr = overburden_pressure(Pover); CHKERRQ(ierr);
  ierr = W.begin_access(); CHKERRQ(ierr);
  ierr = P.begin_access(); CHKERRQ(ierr);
  ierr = Pover.begin_access(); CHKERRQ(ierr);
  ierr = bed->begin_access(); CHKERRQ(ierr);
  ierr = mask->begin_access(); CHKERRQ(ierr);
  ierr = result.begin_access(); CHKERRQ(ierr);
  for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
      if (M.ocean(i,j))
        result(i,j) = Pover(i,j);
      else
        result(i,j) = P(i,j) + rg * ((*bed)(i,j) + W(i,j));
    }
  }
  ierr = W.end_access(); CHKERRQ(ierr);
  ierr = P.end_access(); CHKERRQ(ierr);
  ierr = Pover.end_access(); CHKERRQ(ierr);
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
            maxV_result*secpera, dtCFL/secpera, dtDIFFW/secpera, dtDIFFP/secpera,
            dt_result/secpera, t_current/secpera); CHKERRQ(ierr);
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

    ierr = hydraulic_potential(psi); CHKERRQ(ierr);
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
    const PetscReal dtavyears = (dt/hydrocount)/secpera;
    ierr = verbPrintf(2, grid.com,
      " 'distributed' hydrology summary:\n"
      "     %d hydrology sub-steps with average dt = %.7f years = %.2f s\n"
      "        (average of %.2f steps per CFL time; last max |V| = %.2e m s-1; last max D = %.2e m^2 s-1)\n"
      "     ice free land lost = %.3e kg, ocean lost = %.3e kg\n"
      "     negative bmelt gain = %.3e kg, null strip lost = %.3e kg\n",
      hydrocount, dtavyears, dtavyears * secpera,
      cumratio/hydrocount, maxV, maxD,
      icefreelost, oceanlost, negativegain, nullstriplost); CHKERRQ(ierr);
  }
  return 0;
}

