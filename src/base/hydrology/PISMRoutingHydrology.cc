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

#include "PISMVars.hh"
#include "pism_options.hh"
#include "Mask.hh"
#include "PISMHydrology.hh"
#include "hydrology_diagnostics.hh"


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
                     "thickness of transportable subglacial water layer",
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
                     "new thickness of transportable subglacial water layer during update",
                     "m", ""); CHKERRQ(ierr);
  ierr = Wnew.set_attr("valid_min", 0.0); CHKERRQ(ierr);
  ierr = Wtilnew.create(grid, "Wtilnew_internal", false); CHKERRQ(ierr);
  ierr = Wtilnew.set_attrs("internal",
                     "new thickness of till (subglacial) water layer during update",
                     "m", ""); CHKERRQ(ierr);
  ierr = Wtilnew.set_attr("valid_min", 0.0); CHKERRQ(ierr);

  return 0;
}


PetscErrorCode PISMRoutingHydrology::init(PISMVars &vars) {
  PetscErrorCode ierr;
  ierr = verbPrintf(2, grid.com,
    "* Initializing the routing subglacial hydrology model ...\n"); CHKERRQ(ierr);
  // initialize water layer thickness from the context if present,
  //   otherwise from -i or -boot_file, otherwise with constant value
  bool stripset;
  ierr = PetscOptionsBegin(grid.com, "",
            "Options controlling the 'routing' subglacial hydrology model", ""); CHKERRQ(ierr);
  {
    ierr = PISMOptionsIsSet("-report_mass_accounting",
      "Report to stdout on mass accounting in hydrology models", report_mass_accounting); CHKERRQ(ierr);
    ierr = PISMOptionsReal("-hydrology_null_strip",
                           "set the width, in km, of the strip around the edge of the computational domain in which hydrology is inactivated",
                           stripwidth,stripset); CHKERRQ(ierr);
    if (stripset) stripwidth *= 1.0e3;
  }
  ierr = PetscOptionsEnd(); CHKERRQ(ierr);

  ierr = PISMHydrology::init(vars); CHKERRQ(ierr);

  ierr = init_bwat(vars); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode PISMRoutingHydrology::init_bwat(PISMVars &vars) {
  PetscErrorCode ierr;

  // initialize water layer thickness from the context if present,
  //   otherwise from -i or -boot_file, otherwise with constant value
  bool i, bootstrap;
  ierr = PetscOptionsBegin(grid.com, "",
            "Options for initializing bwat in the 'routing' subglacial hydrology model", ""); CHKERRQ(ierr);
  {
    ierr = PISMOptionsIsSet("-i", "PISM input file", i); CHKERRQ(ierr);
    ierr = PISMOptionsIsSet("-boot_file", "PISM bootstrapping file",
                            bootstrap); CHKERRQ(ierr);
  }
  ierr = PetscOptionsEnd(); CHKERRQ(ierr);

  IceModelVec2S *W_input = dynamic_cast<IceModelVec2S*>(vars.get("bwat"));
  if (W_input != NULL) { // a variable called "bwat" is already in context
    ierr = W.copy_from(*W_input); CHKERRQ(ierr);
  } else if (i || bootstrap) {
    string filename;
    int start;
    ierr = find_pism_input(filename, bootstrap, start); CHKERRQ(ierr);
    if (i) {
      ierr = W.read(filename, start); CHKERRQ(ierr);
    } else {
      ierr = W.regrid(filename,
                      config.get("bootstrapping_bwat_value_no_var")); CHKERRQ(ierr);
    }
  } else {
    ierr = W.set(config.get("bootstrapping_bwat_value_no_var")); CHKERRQ(ierr);
  }

  // however we initialized it, we could be asked to regrid from file
  ierr = regrid(W); CHKERRQ(ierr);
  return 0;
}


void PISMRoutingHydrology::add_vars_to_output(string keyword, set<string> &result) {
  PISMHydrology::add_vars_to_output(keyword, result);
  result.insert("bwat");
}


PetscErrorCode PISMRoutingHydrology::define_variables(set<string> vars, const PIO &nc,
                                                 PISM_IO_Type nctype) {
  PetscErrorCode ierr;
  ierr = PISMHydrology::define_variables(vars, nc, nctype); CHKERRQ(ierr);
  if (set_contains(vars, "bwat")) {
    ierr = W.define(nc, nctype); CHKERRQ(ierr);
  }
  return 0;
}


PetscErrorCode PISMRoutingHydrology::write_variables(set<string> vars, const PIO &nc) {
  PetscErrorCode ierr;
  ierr = PISMHydrology::write_variables(vars, nc); CHKERRQ(ierr);
  if (set_contains(vars, "bwat")) {
    ierr = W.write(nc); CHKERRQ(ierr);
  }
  return 0;
}


void PISMRoutingHydrology::get_diagnostics(map<string, PISMDiagnostic*> &dict,
                                           map<string, PISMTSDiagnostic*> &/*ts_dict*/) {
  // bwat is state
  dict["bwp"] = new PISMHydrology_bwp(this, grid, *variables);
  dict["bwprel"] = new PISMHydrology_bwprel(this, grid, *variables);
  dict["effbwp"] = new PISMHydrology_effbwp(this, grid, *variables);
  dict["hydroinput"] = new PISMHydrology_hydroinput(this, grid, *variables);
  dict["wallmelt"] = new PISMHydrology_wallmelt(this, grid, *variables);
  // add diagnostic that only makes sense if transport is modeled
  dict["bwatvel"] = new PISMRoutingHydrology_bwatvel(this, grid, *variables);
}


//! Check thk >= 0 and fails with message if not satisfied.
PetscErrorCode PISMRoutingHydrology::check_water_thickness_nonnegative(IceModelVec2S &waterthk) {
  PetscErrorCode ierr;
  ierr = waterthk.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (waterthk(i,j) < 0.0) {
        PetscPrintf(grid.com,
           "PISMRoutingHydrology ERROR: disallowed negative water layer thickness\n"
           "    waterthk(i,j) = %.6f m at (i,j)=(%d,%d)\n"
           "ENDING ... \n\n", waterthk(i,j),i,j);
        PISMEnd();
      }
    }
  }
  ierr = waterthk.end_access(); CHKERRQ(ierr);
  return 0;
}


//! Correct the new water thickness based on boundary requirements.  Do mass accounting.
/*!
At ice free locations and ocean locations we require that water thicknesses
(i.e. both the transportable water thickness \f$W\f$ and the till water
thickness \f$W_{til}\f$) be zero at the end of each time step.  Also we require
that any negative water thicknesses be set to zero (i.e. we dp projection to
enforce lower bound).  This method does not enforce any upper bounds.

This method should be called once for each thickness field which needs to be
processed.  This method takes alters the "new" field newthk in-place and sums
the boundary removals.

This method does no reporting at stdout; the calling routine can do that.
 */
PetscErrorCode PISMRoutingHydrology::boundary_mass_changes(
            IceModelVec2S &newthk,
            PetscReal &icefreelost, PetscReal &oceanlost,
            PetscReal &negativegain, PetscReal &nullstriplost) {
  PetscErrorCode ierr;
  PetscReal fresh_water_density = config.get("fresh_water_density");
  PetscReal my_icefreelost = 0.0, my_oceanlost = 0.0, my_negativegain = 0.0;
  MaskQuery M(*mask);

  ierr = newthk.begin_access(); CHKERRQ(ierr);
  ierr = cellarea->begin_access(); CHKERRQ(ierr);
  ierr = mask->begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      const PetscReal dmassdz = (*cellarea)(i,j) * fresh_water_density; // kg m-1
      if (newthk(i,j) < 0.0) {
        my_negativegain += -newthk(i,j) * dmassdz;
        newthk(i,j) = 0.0;
      }
      if (M.ice_free_land(i,j) && (newthk(i,j) > 0.0)) {
        my_icefreelost += newthk(i,j) * dmassdz;
        newthk(i,j) = 0.0;
      }
      if (M.ocean(i,j) && (newthk(i,j) > 0.0)) {
        my_oceanlost += newthk(i,j) * dmassdz;
        newthk(i,j) = 0.0;
      }
    }
  }
  ierr = newthk.end_access(); CHKERRQ(ierr);
  ierr = cellarea->end_access(); CHKERRQ(ierr);
  ierr = mask->end_access(); CHKERRQ(ierr);

  // make global over all proc domains (i.e. whole glacier/ice sheet)
  ierr = PISMGlobalSum(&my_icefreelost, &icefreelost, grid.com); CHKERRQ(ierr);
  ierr = PISMGlobalSum(&my_oceanlost, &oceanlost, grid.com); CHKERRQ(ierr);
  ierr = PISMGlobalSum(&my_negativegain, &negativegain, grid.com); CHKERRQ(ierr);

  if (stripwidth <= 0.0) {
    nullstriplost = 0.0;
    return 0;
  }
  ierr = newthk.begin_access(); CHKERRQ(ierr);
  ierr = cellarea->begin_access(); CHKERRQ(ierr);
  PetscReal my_nullstriplost = 0.0;
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      const PetscReal dmassdz = (*cellarea)(i,j) * fresh_water_density; // kg m-1
      if (in_null_strip(i,j)) {
        my_nullstriplost += newthk(i,j) * dmassdz;
        newthk(i,j) = 0.0;
      }
    }
  }
  ierr = newthk.end_access(); CHKERRQ(ierr);
  ierr = cellarea->end_access(); CHKERRQ(ierr);

  ierr = PISMGlobalSum(&my_nullstriplost, &nullstriplost, grid.com); CHKERRQ(ierr);

  return 0;
}


//! Copies the W variable, the modeled transportable water layer thickness.
PetscErrorCode PISMRoutingHydrology::subglacial_water_thickness(IceModelVec2S &result) {
  PetscErrorCode ierr = W.copy_to(result); CHKERRQ(ierr);
  return 0;
}


//! Returns the (trivial) overburden pressure as the pressure of the transportable water, because this is the model.
PetscErrorCode PISMRoutingHydrology::subglacial_water_pressure(IceModelVec2S &result) {
  PetscErrorCode ierr = overburden_pressure(result); CHKERRQ(ierr);
  return 0;
}


//! Get the hydraulic potential from bedrock topography and current state variables.
/*!
Computes \f$\psi = P + \rho_w g (b + W)\f$ except where floating, where \f$\psi = P_o\f$.
Calls subglacial_water_pressure() method to get water pressure.
 */
PetscErrorCode PISMRoutingHydrology::subglacial_hydraulic_potential(IceModelVec2S &result) {
  PetscErrorCode ierr;

  const PetscReal
    rg = config.get("fresh_water_density") * config.get("standard_gravity");

  ierr = subglacial_water_pressure(result); CHKERRQ(ierr);
  ierr = result.add(rg, (*bed)); CHKERRQ(ierr); // result  <-- P + rhow g b
  ierr = result.add(rg, W); CHKERRQ(ierr);      // result  <-- result + rhow g (b + W)

  // now mask: psi = P_o if ocean
  MaskQuery M(*mask);
  ierr = overburden_pressure(Pover); CHKERRQ(ierr);
  ierr = Pover.begin_access(); CHKERRQ(ierr);
  ierr = mask->begin_access(); CHKERRQ(ierr);
  ierr = result.begin_access(); CHKERRQ(ierr);
  for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
      if (M.ocean(i,j))
        result(i,j) = Pover(i,j);
    }
  }
  ierr = Pover.end_access(); CHKERRQ(ierr);
  ierr = mask->end_access(); CHKERRQ(ierr);
  ierr = result.end_access(); CHKERRQ(ierr);
  return 0;
}


//! Average the regular grid water thickness to values at the center of cell edges.
/*! Uses mask values to avoid averaging using water thickness values from
either ice-free or floating areas. */
PetscErrorCode PISMRoutingHydrology::water_thickness_staggered(IceModelVec2Stag &result) {
  PetscErrorCode ierr;

  MaskQuery M(*mask);
  ierr = mask->begin_access(); CHKERRQ(ierr);
  ierr = W.begin_access(); CHKERRQ(ierr);
  ierr = result.begin_access(); CHKERRQ(ierr);
  for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
      // east
      if (M.grounded_ice(i,j))
        if (M.grounded_ice(i+1,j))
          result(i,j,0) = 0.5 * (W(i,j) + W(i+1,j));
        else
          result(i,j,0) = W(i,j);
      else
        if (M.grounded_ice(i+1,j))
          result(i,j,0) = W(i+1,j);
        else
          result(i,j,0) = 0.0;
      // north
      if (M.grounded_ice(i,j))
        if (M.grounded_ice(i,j+1))
          result(i,j,1) = 0.5 * (W(i,j) + W(i,j+1));
        else
          result(i,j,1) = W(i,j);
      else
        if (M.grounded_ice(i,j+1))
          result(i,j,1) = W(i,j+1);
        else
          result(i,j,1) = 0.0;
    }
  }
  ierr = mask->end_access(); CHKERRQ(ierr);
  ierr = W.end_access(); CHKERRQ(ierr);
  ierr = result.end_access(); CHKERRQ(ierr);
  return 0;
}


//! Compute the nonlinear conductivity at the center of cell edges.
/*!
Computes
    \f[ K = K(W,\nabla P, \nabla b) = k W^{\alpha-1} |\nabla(P+\rho_w g b)|^{\beta-2} \f]
on the staggered grid.  We denote \f$R = P+\rho_w g b\f$ and
\f$\Pi = |\nabla R|^2\f$ internally; the latter is computed on a staggered grid
by a [\ref Mahaffy] -like scheme.  This requires \f$R\f$ to be defined on a box
stencil of width 1.

Also returns the maximum over all staggered points of \f$ K W \f$.
 */
PetscErrorCode PISMRoutingHydrology::conductivity_staggered(
                       IceModelVec2Stag &result, PetscReal &maxKW) {
  PetscErrorCode ierr;
  const PetscReal
    k     = config.get("hydrology_hydraulic_conductivity"),
    alpha = config.get("hydrology_thickness_power_in_flux"),
    beta  = config.get("hydrology_gradient_power_in_flux"),
    rg    = config.get("standard_gravity") * config.get("fresh_water_density");
  if (alpha < 1.0) {
    PetscPrintf(grid.com,
           "PISM ERROR: alpha = %f < 1 which is not allowed\n"
           "ENDING ... \n\n", alpha);
    PISMEnd();
  }

  // the following calculation is bypassed if beta == 2.0 exactly; it puts
  // the squared norm of the gradient of the simplified hydrolic potential
  // temporarily in "result"
  if (beta != 2.0) {
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

  PetscReal betapow = (beta-2.0)/2.0, mymaxKW = 0.0;
  ierr = Wstag.begin_access(); CHKERRQ(ierr);
  ierr = result.begin_access(); CHKERRQ(ierr);
  for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
      for (PetscInt o = 0; o < 2; ++o) {
        PetscReal Ktmp = k * pow(Wstag(i,j,o),alpha-1.0);
        if (beta < 2.0) {
          // regularize negative power |\grad psi|^{beta-2} by adding eps because
          //   large head gradient might be 10^7 Pa per 10^4 m or 10^3 Pa/m
          const PetscReal eps = 1.0;   // Pa m-1
          result(i,j,o) = Ktmp * pow(result(i,j,o) + eps * eps,betapow);
        } else if (beta > 2.0) {
          result(i,j,o) = Ktmp * pow(result(i,j,o),betapow);
        } else { // beta == 2.0
          result(i,j,o) = Ktmp;
        }
        mymaxKW = PetscMax( mymaxKW, result(i,j,o) * Wstag(i,j,o) );
      }
    }
  }
  ierr = Wstag.end_access(); CHKERRQ(ierr);
  ierr = result.end_access(); CHKERRQ(ierr);

  ierr = PISMGlobalMax(&mymaxKW, &maxKW, grid.com); CHKERRQ(ierr);

  return 0;
}


//! Compute the wall melt rate which comes from (turbulent) dissipation of flow energy.
/*!
This code fills `result` with
    \f[ \frac{m_{wall}}{\rho_w} = - \frac{1}{L \rho_w} \mathbf{q} \cdot \nabla \psi = \left(\frac{k}{L \rho_w}\right) W^\alpha |\nabla R|^\beta \f]
where \f$R = P+\rho_w g b\f$.

Note that conductivity_staggered() computes the related quantity
\f$K = k W^{\alpha-1} |\nabla R|^{\beta-2}\f$ on the staggered grid, but
contriving to reuse that code would be inefficient because of the
staggered-versus-regular change.

At the current state of the code, this is a diagnostic calculation only.
 */
PetscErrorCode PISMRoutingHydrology::wall_melt(IceModelVec2S &result) {
  PetscErrorCode ierr;

  const PetscReal
    k     = config.get("hydrology_hydraulic_conductivity"),
    L     = config.get("water_latent_heat_fusion"),
    alpha = config.get("hydrology_thickness_power_in_flux"),
    beta  = config.get("hydrology_gradient_power_in_flux"),
    rhow  = config.get("standard_gravity"),
    g     = config.get("fresh_water_density"),
    rg    = rhow * g,
    CC    = k / (L * rhow);

  // FIXME:  could be scaled with overall factor hydrology_coefficient_wall_melt ?
  if (alpha < 1.0) {
    PetscPrintf(grid.com,
           "PISM ERROR: alpha = %f < 1 which is not allowed\n"
           "ENDING ... \n\n", alpha);
    PISMEnd();
  }

  ierr = subglacial_water_pressure(R); CHKERRQ(ierr);  // yes, it updates ghosts
  ierr = R.add(rg, (*bed)); CHKERRQ(ierr); // R  <-- P + rhow g b
  ierr = R.update_ghosts(); CHKERRQ(ierr);

  PetscReal dRdx, dRdy;
  ierr = R.begin_access(); CHKERRQ(ierr);
  ierr = W.begin_access(); CHKERRQ(ierr);
  ierr = result.begin_access(); CHKERRQ(ierr);
  for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
      if (W(i,j) > 0.0) {
        dRdx = 0.0;
        if (W(i+1,j) > 0.0) {
          dRdx =  ( R(i+1,j) - R(i,j) ) / (2.0 * grid.dx);
        }
        if (W(i-1,j) > 0.0) {
          dRdx += ( R(i,j) - R(i-1,j) ) / (2.0 * grid.dx);
        }
        dRdy = 0.0;
        if (W(i,j+1) > 0.0) {
          dRdy =  ( R(i,j+1) - R(i,j) ) / (2.0 * grid.dy);
        }
        if (W(i,j-1) > 0.0) {
          dRdy += ( R(i,j) - R(i,j-1) ) / (2.0 * grid.dy);
        }
        result(i,j) = CC * pow(W(i,j),alpha) * pow(dRdx * dRdx + dRdy * dRdy, beta/2.0);
      } else
        result(i,j) = 0.0;
    }
  }
  ierr = R.end_access(); CHKERRQ(ierr);
  ierr = W.end_access(); CHKERRQ(ierr);
  ierr = result.end_access(); CHKERRQ(ierr);

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
    dtmax = config.get("hydrology_maximum_time_step_years",
                       "years", "seconds"),
    rg    = config.get("standard_gravity") * config.get("fresh_water_density");
  PetscReal tmp[2];
  ierr = V.absmaxcomponents(tmp); CHKERRQ(ierr); // V could be zero if P is constant and bed is flat
  maxV_result = sqrt(tmp[0]*tmp[0] + tmp[1]*tmp[1]);
  maxD_result = rg * maxKW;
  dtCFL_result = 0.5 / (tmp[0]/grid.dx + tmp[1]/grid.dy); // FIXME: is regularization needed?
  dtDIFFW_result = 1.0/(grid.dx*grid.dx) + 1.0/(grid.dy*grid.dy);
  dtDIFFW_result = 0.25 / (maxD_result * dtDIFFW_result);
  // dt = min { te-t, dtmax, dtCFL, dtDIFFW }
  dt_result = PetscMin(t_end - t_current, dtmax);
  dt_result = PetscMin(dt_result, dtCFL_result);
  dt_result = PetscMin(dt_result, dtDIFFW_result);
  return 0;
}


//! The computation of Wtilnew, called by update().
/*!
Does an implicit (backward Euler) step of the integration
  \f[ \frac{\partial W_{til}}{\partial t} = \mu \left(\min\{\omega W,W_{til}^{max} - W_{til}\right)\f]
where \f$\mu=\f$`hydrology_tillwat_rate`, \f$\omega=\f$`hydrology_tillwat_transfer_proportion`,
and \f$W_{til}^{max}\f$=`hydrology_tillwat_max`.  There is no time-step restriction.
The solution satisfies the inequalities
  \f[ 0 \le W_{til} \le W_{til}^{max}.\f]
 */
PetscErrorCode PISMRoutingHydrology::raw_update_Wtil(PetscReal hdt) {
    PetscErrorCode ierr;
    PetscReal tmpmin;
    const PetscReal Wtilmax  = config.get("hydrology_tillwat_max"),
                    mu       = config.get("hydrology_tillwat_rate"),
                    omega    = config.get("hydrology_tillwat_transfer_proportion"),
                    denom    = 1.0 + mu * hdt;
    if ((Wtilmax < 0.0) || (mu < 0.0) || (omega < 0.0)) {
      PetscPrintf(grid.com,
         "PISMRoutingHydrology ERROR: scalar config parameter is negative in raw_update_Wtil();\n"
         "            this is not allowed ... ENDING ... \n\n");
      PISMEnd();
    }
    ierr = W.begin_access(); CHKERRQ(ierr);
    ierr = Wtil.begin_access(); CHKERRQ(ierr);
    ierr = Wtilnew.begin_access(); CHKERRQ(ierr);
    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
        tmpmin = PetscMin(omega * W(i,j), Wtilmax);
        Wtilnew(i,j) = (Wtil(i,j) + mu * hdt * tmpmin) / denom;
      }
    }
    ierr = W.end_access(); CHKERRQ(ierr);
    ierr = Wtil.end_access(); CHKERRQ(ierr);
    ierr = Wtilnew.end_access(); CHKERRQ(ierr);

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
    ierr = Wtil.begin_access(); CHKERRQ(ierr);
    ierr = Wtilnew.begin_access(); CHKERRQ(ierr);
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
        Wnew(i,j) = W(i,j) - Wtilnew(i,j) + Wtil(i,j)
                    + hdt * (- divadflux + diffW + total_input(i,j));
      }
    }
    ierr = W.end_access(); CHKERRQ(ierr);
    ierr = Wtil.end_access(); CHKERRQ(ierr);
    ierr = Wtilnew.end_access(); CHKERRQ(ierr);
    ierr = Wstag.end_access(); CHKERRQ(ierr);
    ierr = Kstag.end_access(); CHKERRQ(ierr);
    ierr = Qstag.end_access(); CHKERRQ(ierr);
    ierr = total_input.end_access(); CHKERRQ(ierr);
    ierr = Wnew.end_access(); CHKERRQ(ierr);

    return 0;
}


//! Update the model state variables W and Wtil by applying the subglacial hydrology model equations.
/*!
Runs the hydrology model from time icet to time icet + icedt.  Here [icet,icedt]
is generally on the order of months to years.  This hydrology model will take its
own shorter time steps, perhaps hours to weeks.

For updating W = `bwat`, calls raw_update_W().  For updating Wtil = `tillwat`,
calls raw_update_Wtil().
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

#if (PISM_DEBUG==1)
    ierr = check_water_thickness_nonnegative(W); CHKERRQ(ierr);
    ierr = check_Wtil_bounds(); CHKERRQ(ierr);
#endif

    // update Wtilnew (the actual step) from W and Wtil
    ierr = raw_update_Wtil(hdt); CHKERRQ(ierr);
    ierr = boundary_mass_changes(Wtilnew, delta_icefree, delta_ocean,
                                 delta_neggain, delta_nullstrip); CHKERRQ(ierr);
    icefreelost  += delta_icefree;
    oceanlost    += delta_ocean;
    negativegain += delta_neggain;
    nullstriplost+= delta_nullstrip;

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

    // update Wnew (the actual step) from W, Wtil, Wtilnew, Wstag, Qstag, total_input
    ierr = raw_update_W(hdt); CHKERRQ(ierr);
    ierr = boundary_mass_changes(Wnew, delta_icefree, delta_ocean,
                                 delta_neggain, delta_nullstrip); CHKERRQ(ierr);
    icefreelost  += delta_icefree;
    oceanlost    += delta_ocean;
    negativegain += delta_neggain;
    nullstriplost+= delta_nullstrip;

    // transfer new into old
    ierr = Wnew.update_ghosts(W); CHKERRQ(ierr);
    ierr = Wtilnew.copy_to(Wtil); CHKERRQ(ierr);

    ht += hdt;
  } // end of hydrology model time-stepping loop

  if (report_mass_accounting) {
    ierr = verbPrintf(2, grid.com,
                      " 'routing' hydrology summary:\n"
                      "     %d hydrology sub-steps with average dt = %.6f years = %.2f s\n"
                      "        (max |V| = %.2e m s-1; max D = %.2e m^2 s-1)\n"
                      "     ice free land lost = %.3e kg, ocean lost = %.3e kg\n"
                      "     negative bmelt gain = %.3e kg, null strip lost = %.3e kg\n",
                      hydrocount, grid.convert(dt/hydrocount, "seconds", "years"), dt/hydrocount,
                      maxV, maxD,
                      icefreelost, oceanlost,
                      negativegain, nullstriplost); CHKERRQ(ierr);
  }
  return 0;
}


PISMRoutingHydrology_bwatvel::PISMRoutingHydrology_bwatvel(PISMRoutingHydrology *m, IceGrid &g, PISMVars &my_vars)
    : PISMDiag<PISMRoutingHydrology>(m, g, my_vars) {

  // set metadata:
  dof = 2;
  vars.resize(dof, NCSpatialVariable(g.get_unit_system()));
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
