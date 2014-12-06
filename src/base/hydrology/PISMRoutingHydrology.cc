// Copyright (C) 2012-2014 PISM Authors
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
#include "error_handling.hh"

namespace pism {

RoutingHydrology::RoutingHydrology(IceGrid &g, const Config &conf)
    : Hydrology(g, conf)
{
  stripwidth = config.get("hydrology_null_strip_width");

  // model state variables; need ghosts
  W.create(grid, "bwat", WITH_GHOSTS, 1);
  W.set_attrs("model_state",
              "thickness of transportable subglacial water layer",
              "m", "");
  W.metadata().set_double("valid_min", 0.0);

  // auxiliary variables which NEED ghosts
  Wstag.create(grid, "W_staggered", WITH_GHOSTS, 1);
  Wstag.set_attrs("internal",
                  "cell face-centered (staggered) values of water layer thickness",
                  "m", "");
  Wstag.metadata().set_double("valid_min", 0.0);
  Kstag.create(grid, "K_staggered", WITH_GHOSTS, 1);
  Kstag.set_attrs("internal",
                  "cell face-centered (staggered) values of nonlinear conductivity",
                  "", "");
  Kstag.metadata().set_double("valid_min", 0.0);
  Qstag.create(grid, "advection_flux", WITH_GHOSTS, 1);
  Qstag.set_attrs("internal",
                  "cell face-centered (staggered) components of advective subglacial water flux",
                  "m2 s-1", "");
  R.create(grid, "potential_workspace", WITH_GHOSTS, 1); // box stencil used
  R.set_attrs("internal",
              "work space for modeled subglacial water hydraulic potential",
              "Pa", "");

  // auxiliary variables which do not need ghosts
  Pover.create(grid, "overburden_pressure_internal", WITHOUT_GHOSTS);
  Pover.set_attrs("internal",
                  "overburden pressure",
                  "Pa", "");
  Pover.metadata().set_double("valid_min", 0.0);
  V.create(grid, "water_velocity", WITHOUT_GHOSTS);
  V.set_attrs("internal",
              "cell face-centered (staggered) components of water velocity in subglacial water layer",
              "m s-1", "");

  // temporaries during update; do not need ghosts
  Wnew.create(grid, "Wnew_internal", WITHOUT_GHOSTS);
  Wnew.set_attrs("internal",
                 "new thickness of transportable subglacial water layer during update",
                 "m", "");
  Wnew.metadata().set_double("valid_min", 0.0);
  Wtilnew.create(grid, "Wtilnew_internal", WITHOUT_GHOSTS);
  Wtilnew.set_attrs("internal",
                    "new thickness of till (subglacial) water layer during update",
                    "m", "");
  Wtilnew.metadata().set_double("valid_min", 0.0);
}

RoutingHydrology::~RoutingHydrology() {
  // empty
}


void RoutingHydrology::init(Vars &vars) {
  verbPrintf(2, grid.com,
             "* Initializing the routing subglacial hydrology model ...\n");
  // initialize water layer thickness from the context if present,
  //   otherwise from -i or -boot_file, otherwise with constant value
  bool stripset;
  {
    OptionsIsSet("-report_mass_accounting",
                 "Report to stdout on mass accounting in hydrology models",
                 report_mass_accounting);

    stripwidth = grid.convert(stripwidth, "m", "km");
    OptionsReal("-hydrology_null_strip",
                "set the width, in km, of the strip around the edge"
                " of the computational domain in which hydrology is inactivated",
                stripwidth, stripset);
    if (stripset) {
      stripwidth = grid.convert(stripwidth, "km", "m");
    }
  }

  Hydrology::init(vars);

  init_bwat(vars);
}


void RoutingHydrology::init_bwat(Vars &vars) {

  // initialize water layer thickness from the context if present,
  //   otherwise from -i or -boot_file, otherwise with constant value
  bool i, bootstrap;
  {
    OptionsIsSet("-i", "PISM input file", i);
    OptionsIsSet("-boot_file", "PISM bootstrapping file",
                 bootstrap);
  }

  const PetscReal bwatdefault = config.get("bootstrapping_bwat_value_no_var");
  IceModelVec2S *W_input = NULL;
  try {
    // FIXME: this is not an "exceptional" situation...
    // a variable called "bwat" is already in context
    W_input = vars.get_2d_scalar("bwat");
    W.copy_from(*W_input);
  } catch (RuntimeError) {
    if (i || bootstrap) {
      std::string filename;
      int start;
      find_pism_input(filename, bootstrap, start);
      if (i) {
        PIO nc(grid, "guess_mode");
        nc.open(filename, PISM_READONLY);
        bool bwat_exists = nc.inq_var("bwat");
        if (bwat_exists == true) {
          W.read(filename, start);
        } else {
          verbPrintf(2, grid.com,
                     "PISM WARNING: bwat for hydrology model not found in '%s'."
                     "  Setting it to %.2f ...\n",
                     filename.c_str(), bwatdefault);
          W.set(bwatdefault);
        }
        nc.close();
      } else {
        W.regrid(filename, OPTIONAL, bwatdefault);
      }
    } else { // not sure if this case can be reached, but ...
      W.set(bwatdefault);
    }
  }

  // however we initialized it, we could be asked to regrid from file
  regrid("RoutingHydrology", &W);
}


void RoutingHydrology::add_vars_to_output(const std::string &keyword, std::set<std::string> &result) {
  Hydrology::add_vars_to_output(keyword, result);
  result.insert("bwat");
}


void RoutingHydrology::define_variables(const std::set<std::string> &vars, const PIO &nc,
                                                 IO_Type nctype) {
  Hydrology::define_variables(vars, nc, nctype);
  if (set_contains(vars, "bwat")) {
    W.define(nc, nctype);
  }
}


void RoutingHydrology::write_variables(const std::set<std::string> &vars, const PIO &nc) {
  Hydrology::write_variables(vars, nc);
  if (set_contains(vars, "bwat")) {
    W.write(nc);
  }
}


void RoutingHydrology::get_diagnostics(std::map<std::string, Diagnostic*> &dict,
                                           std::map<std::string, TSDiagnostic*> &/*ts_dict*/) {
  // bwat is state
  dict["bwp"] = new Hydrology_bwp(this, grid, *variables);
  dict["bwprel"] = new Hydrology_bwprel(this, grid, *variables);
  dict["effbwp"] = new Hydrology_effbwp(this, grid, *variables);
  dict["hydrobmelt"] = new Hydrology_hydrobmelt(this, grid, *variables);
  dict["hydroinput"] = new Hydrology_hydroinput(this, grid, *variables);
  dict["wallmelt"] = new Hydrology_wallmelt(this, grid, *variables);
  // add diagnostic that only makes sense if transport is modeled
  dict["bwatvel"] = new RoutingHydrology_bwatvel(this, grid, *variables);
}


//! Check thk >= 0 and fails with message if not satisfied.
void RoutingHydrology::check_water_thickness_nonnegative(IceModelVec2S &waterthk) {
  IceModelVec::AccessList list;
  list.add(waterthk);

  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (waterthk(i,j) < 0.0) {
      throw RuntimeError::formatted("RoutingHydrology: disallowed negative water layer thickness\n"
                                    "waterthk(i,j) = %.6f m at (i,j)=(%d,%d)",
                                    waterthk(i,j),i,j);
    }
  }
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
void RoutingHydrology::boundary_mass_changes(IceModelVec2S &newthk,
                                             double &icefreelost, double &oceanlost,
                                             double &negativegain, double &nullstriplost) {
  double fresh_water_density = config.get("fresh_water_density");
  double my_icefreelost = 0.0, my_oceanlost = 0.0, my_negativegain = 0.0;
  MaskQuery M(*mask);

  IceModelVec::AccessList list;
  list.add(newthk);
  list.add(*cellarea);
  list.add(*mask);

  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    const double dmassdz = (*cellarea)(i,j) * fresh_water_density; // kg m-1
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

  // make global over all proc domains (i.e. whole glacier/ice sheet)
  GlobalSum(grid.com, &my_icefreelost,  &icefreelost);
  GlobalSum(grid.com, &my_oceanlost,  &oceanlost);
  GlobalSum(grid.com, &my_negativegain,  &negativegain);

  if (stripwidth <= 0.0) {
    nullstriplost = 0.0;
    return;
  }

  double my_nullstriplost = 0.0;
  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    const double dmassdz = (*cellarea)(i,j) * fresh_water_density; // kg m-1
    if (grid.in_null_strip(i, j, stripwidth)) {
      my_nullstriplost += newthk(i,j) * dmassdz;
      newthk(i,j) = 0.0;
    }
  }

  GlobalSum(grid.com, &my_nullstriplost,  &nullstriplost);
}


//! Copies the W variable, the modeled transportable water layer thickness.
void RoutingHydrology::subglacial_water_thickness(IceModelVec2S &result) {
  W.copy_to(result);
}


//! Returns the (trivial) overburden pressure as the pressure of the transportable water, because this is the model.
void RoutingHydrology::subglacial_water_pressure(IceModelVec2S &result) {
  overburden_pressure(result);
}


//! Get the hydraulic potential from bedrock topography and current state variables.
/*!
Computes \f$\psi = P + \rho_w g (b + W)\f$ except where floating, where \f$\psi = P_o\f$.
Calls subglacial_water_pressure() method to get water pressure.
 */
void RoutingHydrology::subglacial_hydraulic_potential(IceModelVec2S &result) {

  const double
    rg = config.get("fresh_water_density") * config.get("standard_gravity");

  subglacial_water_pressure(result);
  result.add(rg, (*bed)); // result  <-- P + rhow g b
  result.add(rg, W);      // result  <-- result + rhow g (b + W)

  // now mask: psi = P_o if ocean
  MaskQuery M(*mask);
  overburden_pressure(Pover);

  IceModelVec::AccessList list;
  list.add(Pover);
  list.add(*mask);
  list.add(result);

  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (M.ocean(i,j)) {
      result(i,j) = Pover(i,j);
    }
  }
}


//! Average the regular grid water thickness to values at the center of cell edges.
/*! Uses mask values to avoid averaging using water thickness values from
either ice-free or floating areas. */
void RoutingHydrology::water_thickness_staggered(IceModelVec2Stag &result) {
  MaskQuery M(*mask);

  IceModelVec::AccessList list;
  list.add(*mask);
  list.add(W);
  list.add(result);

  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    // east
    if (M.grounded_ice(i,j)) {
      if (M.grounded_ice(i+1,j)) {
        result(i,j,0) = 0.5 * (W(i,j) + W(i+1,j));
      } else {
        result(i,j,0) = W(i,j);
      }
    } else {
      if (M.grounded_ice(i+1,j)) {
        result(i,j,0) = W(i+1,j);
      } else {
        result(i,j,0) = 0.0;
      }
    }
    // north
    if (M.grounded_ice(i,j)) {
      if (M.grounded_ice(i,j+1)) {
        result(i,j,1) = 0.5 * (W(i,j) + W(i,j+1));
      } else {
        result(i,j,1) = W(i,j);
      }
    } else {
      if (M.grounded_ice(i,j+1)) {
        result(i,j,1) = W(i,j+1);
      } else {
        result(i,j,1) = 0.0;
      }
    }
  }
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
void RoutingHydrology::conductivity_staggered(IceModelVec2Stag &result,
                                                        double &maxKW) {
  const double
    k     = config.get("hydrology_hydraulic_conductivity"),
    alpha = config.get("hydrology_thickness_power_in_flux"),
    beta  = config.get("hydrology_gradient_power_in_flux"),
    rg    = config.get("standard_gravity") * config.get("fresh_water_density");
  if (alpha < 1.0) {
    throw RuntimeError::formatted("alpha = %f < 1 which is not allowed", alpha);
  }

  IceModelVec::AccessList list(result);

  // the following calculation is bypassed if beta == 2.0 exactly; it puts
  // the squared norm of the gradient of the simplified hydrolic potential
  // temporarily in "result"
  if (beta != 2.0) {
    subglacial_water_pressure(R);  // yes, it updates ghosts
    R.add(rg, (*bed)); // R  <-- P + rhow g b
    R.update_ghosts();

    list.add(R);
    for (Points p(grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      double dRdx, dRdy;
      dRdx = (R(i+1,j) - R(i,j)) / grid.dx();
      dRdy = (R(i+1,j+1) + R(i,j+1) - R(i+1,j-1) - R(i,j-1)) / (4.0 * grid.dy());
      result(i,j,0) = dRdx * dRdx + dRdy * dRdy;
      dRdx = (R(i+1,j+1) + R(i+1,j) - R(i-1,j+1) - R(i-1,j)) / (4.0 * grid.dx());
      dRdy = (R(i,j+1) - R(i,j)) / grid.dy();
      result(i,j,1) = dRdx * dRdx + dRdy * dRdy;
    }
  }

  double betapow = (beta-2.0)/2.0, mymaxKW = 0.0;

  list.add(Wstag);
  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    for (int o = 0; o < 2; ++o) {
      double Ktmp = k * pow(Wstag(i,j,o),alpha-1.0);
      if (beta < 2.0) {
        // regularize negative power |\grad psi|^{beta-2} by adding eps because
        //   large head gradient might be 10^7 Pa per 10^4 m or 10^3 Pa/m
        const double eps = 1.0;   // Pa m-1
        result(i,j,o) = Ktmp * pow(result(i,j,o) + eps * eps,betapow);
      } else if (beta > 2.0) {
        result(i,j,o) = Ktmp * pow(result(i,j,o),betapow);
      } else { // beta == 2.0
        result(i,j,o) = Ktmp;
      }
      mymaxKW = std::max(mymaxKW, result(i,j,o) * Wstag(i,j,o));
    }
  }

  GlobalMax(grid.com, &mymaxKW,  &maxKW);
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
void RoutingHydrology::wall_melt(IceModelVec2S &result) {

  const double
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
    throw RuntimeError::formatted("alpha = %f < 1 which is not allowed", alpha);
  }

  subglacial_water_pressure(R);  // yes, it updates ghosts
  R.add(rg, (*bed)); // R  <-- P + rhow g b
  R.update_ghosts();

  IceModelVec::AccessList list;
  list.add(R);
  list.add(W);
  list.add(result);

  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();
    double dRdx, dRdy;

    if (W(i,j) > 0.0) {
      dRdx = 0.0;
      if (W(i+1,j) > 0.0) {
        dRdx =  (R(i+1,j) - R(i,j)) / (2.0 * grid.dx());
      }
      if (W(i-1,j) > 0.0) {
        dRdx += (R(i,j) - R(i-1,j)) / (2.0 * grid.dx());
      }
      dRdy = 0.0;
      if (W(i,j+1) > 0.0) {
        dRdy =  (R(i,j+1) - R(i,j)) / (2.0 * grid.dy());
      }
      if (W(i,j-1) > 0.0) {
        dRdy += (R(i,j) - R(i,j-1)) / (2.0 * grid.dy());
      }
      result(i,j) = CC * pow(W(i,j),alpha) * pow(dRdx * dRdx + dRdy * dRdy, beta/2.0);
    } else {
      result(i,j) = 0.0;
    }
  }
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
void RoutingHydrology::velocity_staggered(IceModelVec2Stag &result) {
  const double  rg = config.get("standard_gravity") * config.get("fresh_water_density");
  double dbdx, dbdy, dPdx, dPdy;

  subglacial_water_pressure(R);  // R=P; yes, it updates ghosts

  IceModelVec::AccessList list;
  list.add(R);
  list.add(Wstag);
  list.add(Kstag);
  list.add(*bed);
  list.add(result);

  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (Wstag(i,j,0) > 0.0) {
      dPdx = (R(i+1,j) - R(i,j)) / grid.dx();
      dbdx = ((*bed)(i+1,j) - (*bed)(i,j)) / grid.dx();
      result(i,j,0) = - Kstag(i,j,0) * (dPdx + rg * dbdx);
    } else {
      result(i,j,0) = 0.0;
    }
    
    if (Wstag(i,j,1) > 0.0) {
      dPdy = (R(i,j+1) - R(i,j)) / grid.dy();
      dbdy = ((*bed)(i,j+1) - (*bed)(i,j)) / grid.dy();
      result(i,j,1) = - Kstag(i,j,1) * (dPdy + rg * dbdy);
    } else {
      result(i,j,1) = 0.0;
    }

    if (grid.in_null_strip(i,j, stripwidth) || grid.in_null_strip(i+1,j, stripwidth)) {
      result(i,j,0) = 0.0;
    }

    if (grid.in_null_strip(i,j, stripwidth) || grid.in_null_strip(i,j+1, stripwidth)) {
      result(i,j,1) = 0.0;
    }
  }
}


//! Compute Q = V W at edge-centers (staggered grid) by first-order upwinding.
/*!
The field W must have valid ghost values, but V does not need them.

FIXME:  This could be re-implemented using the Koren (1993) flux-limiter.
 */
void RoutingHydrology::advective_fluxes(IceModelVec2Stag &result) {
  IceModelVec::AccessList list;
  list.add(W);
  list.add(V);
  list.add(result);

  assert(W.get_stencil_width() >= 1);

  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    result(i,j,0) = (V(i,j,0) >= 0.0) ? V(i,j,0) * W(i,j) :  V(i,j,0) * W(i+1,j);
    result(i,j,1) = (V(i,j,1) >= 0.0) ? V(i,j,1) * W(i,j) :  V(i,j,1) * W(i,  j+1);
  }
}


//! Compute the adaptive time step for evolution of W.
void RoutingHydrology::adaptive_for_W_evolution(double t_current, double t_end, double maxKW,
                                                double &dt_result,
                                                double &maxV_result, double &maxD_result,
                                                double &dtCFL_result, double &dtDIFFW_result) {
  const double
    dtmax = config.get("hydrology_maximum_time_step_years",
                       "years", "seconds"),
    rg    = config.get("standard_gravity") * config.get("fresh_water_density");
  double tmp[2];
  V.absmaxcomponents(tmp); // V could be zero if P is constant and bed is flat
  maxV_result = sqrt(tmp[0]*tmp[0] + tmp[1]*tmp[1]);
  maxD_result = rg * maxKW;
  dtCFL_result = 0.5 / (tmp[0]/grid.dx() + tmp[1]/grid.dy()); // FIXME: is regularization needed?
  dtDIFFW_result = 1.0/(grid.dx()*grid.dx()) + 1.0/(grid.dy()*grid.dy());
  dtDIFFW_result = 0.25 / (maxD_result * dtDIFFW_result);
  // dt = min { te-t, dtmax, dtCFL, dtDIFFW }
  dt_result = std::min(t_end - t_current, dtmax);
  dt_result = std::min(dt_result, dtCFL_result);
  dt_result = std::min(dt_result, dtDIFFW_result);
}


//! The computation of Wtilnew, called by update().
/*!
Does a step of the trivial integration
  \f[ \frac{\partial W_{til}}{\partial t} = \frac{m}{\rho_w} - C\f]
where \f$C=\f$`hydrology_tillwat_decay_rate`.  Enforces bounds
\f$0 \le W_{til} \le W_{til}^{max}\f$ where the upper bound is
`hydrology_tillwat_max`.  Here \f$m/\rho_w\f$ is `total_input`.

Compare NullTransportHydrology::update().  The current code is not quite "code
duplication" because the code here: (1) computes `Wtilnew` instead of updating
`Wtil` in place; (2) uses time steps determined by the rest of the
RoutingHydrology model; (3) does not check mask because the boundary_mass_changes()
call addresses that.  Otherwise this is the same physical model with the
same configurable parameters.
 */
void RoutingHydrology::raw_update_Wtil(double hdt) {
  const double tillwat_max = config.get("hydrology_tillwat_max"),
               C           = config.get("hydrology_tillwat_decay_rate");

  IceModelVec::AccessList list;
  list.add(Wtil);
  list.add(Wtilnew);
  list.add(total_input);

  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    Wtilnew(i,j) = Wtil(i,j) + hdt * (total_input(i,j) - C);
    Wtilnew(i,j) = std::min(std::max(0.0, Wtilnew(i,j)), tillwat_max);
  }
}


//! The computation of Wnew, called by update().
void RoutingHydrology::raw_update_W(double hdt) {
  const double
    wux  = 1.0 / (grid.dx() * grid.dx()),
    wuy  = 1.0 / (grid.dy() * grid.dy()),
    rg   = config.get("standard_gravity") * config.get("fresh_water_density");
  double divadflux, diffW;

  IceModelVec::AccessList list;
  list.add(W);
  list.add(Wtil);
  list.add(Wtilnew);
  list.add(Wstag);
  list.add(Kstag);
  list.add(Qstag);
  list.add(total_input);
  list.add(Wnew);

  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    divadflux =   (Qstag(i,j,0) - Qstag(i-1,j  ,0)) / grid.dx()
      + (Qstag(i,j,1) - Qstag(i,  j-1,1)) / grid.dy();
    const double  De = rg * Kstag(i,  j,0) * Wstag(i,  j,0),
      Dw = rg * Kstag(i-1,j,0) * Wstag(i-1,j,0),
      Dn = rg * Kstag(i,j  ,1) * Wstag(i,j  ,1),
      Ds = rg * Kstag(i,j-1,1) * Wstag(i,j-1,1);
    diffW =   wux * (De * (W(i+1,j) - W(i,j)) - Dw * (W(i,j) - W(i-1,j)))
      + wuy * (Dn * (W(i,j+1) - W(i,j)) - Ds * (W(i,j) - W(i,j-1)));
    Wnew(i,j) = W(i,j) - Wtilnew(i,j) + Wtil(i,j)
      + hdt * (- divadflux + diffW + total_input(i,j));
  }
}


//! Update the model state variables W and Wtil by applying the subglacial hydrology model equations.
/*!
Runs the hydrology model from time icet to time icet + icedt.  Here [icet,icedt]
is generally on the order of months to years.  This hydrology model will take its
own shorter time steps, perhaps hours to weeks.

For updating W = `bwat`, calls raw_update_W().  For updating Wtil = `tillwat`,
calls raw_update_Wtil().
 */
void RoutingHydrology::update(double icet, double icedt) {

  // if asked for the identical time interval versus last time, then
  //   do nothing; otherwise assume that [my_t,my_t+my_dt] is the time
  //   interval on which we are solving
  if ((fabs(icet - m_t) < 1e-12) && (fabs(icedt - m_dt) < 1e-12)) {
    return;
  }
  // update Component times: t = current time, t+dt = target time
  m_t = icet;
  m_dt = icedt;

  if (config.get("hydrology_tillwat_max") < 0.0) {
    throw RuntimeError("RoutingHydrology: hydrology_tillwat_max is negative.\n"
                       "This is not allowed.");
  }

  // make sure W has valid ghosts before starting hydrology steps
  W.update_ghosts();

  MaskQuery M(*mask);
  double ht = m_t, hdt, // hydrology model time and time step
            maxKW, maxV, maxD, dtCFL, dtDIFFW;
  double icefreelost = 0, oceanlost = 0, negativegain = 0, nullstriplost = 0,
            delta_icefree, delta_ocean, delta_neggain, delta_nullstrip;
  int hydrocount = 0; // count hydrology time steps

  while (ht < m_t + m_dt) {
    hydrocount++;

#if (PISM_DEBUG==1)
    check_water_thickness_nonnegative(W);
    check_Wtil_bounds();
#endif

    water_thickness_staggered(Wstag);
    Wstag.update_ghosts();

    conductivity_staggered(Kstag,maxKW);
    Kstag.update_ghosts();

    velocity_staggered(V);

    // to get Qstag, W needs valid ghosts
    advective_fluxes(Qstag);
    Qstag.update_ghosts();

    adaptive_for_W_evolution(ht, m_t+m_dt, maxKW,
                             hdt, maxV, maxD, dtCFL, dtDIFFW);

    if ((inputtobed != NULL) || (hydrocount==1)) {
      get_input_rate(ht,hdt,total_input);
    }

    // update Wtilnew from Wtil
    raw_update_Wtil(hdt);
    boundary_mass_changes(Wtilnew, delta_icefree, delta_ocean,
                          delta_neggain, delta_nullstrip);
    icefreelost  += delta_icefree;
    oceanlost    += delta_ocean;
    negativegain += delta_neggain;
    nullstriplost+= delta_nullstrip;

    // update Wnew from W, Wtil, Wtilnew, Wstag, Qstag, total_input
    raw_update_W(hdt);
    boundary_mass_changes(Wnew, delta_icefree, delta_ocean,
                          delta_neggain, delta_nullstrip);
    icefreelost  += delta_icefree;
    oceanlost    += delta_ocean;
    negativegain += delta_neggain;
    nullstriplost+= delta_nullstrip;

    // transfer new into old
    Wnew.update_ghosts(W);
    Wtilnew.copy_to(Wtil);

    ht += hdt;
  } // end of hydrology model time-stepping loop

  if (report_mass_accounting) {
    verbPrintf(2, grid.com,
               " 'routing' hydrology summary:\n"
               "     %d hydrology sub-steps with average dt = %.6f years = %.2f s\n"
               "        (max |V| = %.2e m s-1; max D = %.2e m^2 s-1)\n"
               "     ice free land loss = %.3e kg, ocean loss = %.3e kg\n"
               "     negative bmelt gain = %.3e kg, null strip loss = %.3e kg\n",
               hydrocount, grid.convert(m_dt/hydrocount, "seconds", "years"), m_dt/hydrocount,
               maxV, maxD,
               icefreelost, oceanlost,
               negativegain, nullstriplost);
  }
}


RoutingHydrology_bwatvel::RoutingHydrology_bwatvel(RoutingHydrology *m, IceGrid &g, Vars &my_vars)
    : Diag<RoutingHydrology>(m, g, my_vars) {

  // set metadata:
  dof = 2;
  vars.push_back(NCSpatialVariable(grid.get_unit_system(), "bwatvel[0]", grid));
  vars.push_back(NCSpatialVariable(grid.get_unit_system(), "bwatvel[1]", grid));

  set_attrs("velocity of water in subglacial layer, i-offset", "",
            "m s-1", "m year-1", 0);
  set_attrs("velocity of water in subglacial layer, j-offset", "",
            "m s-1", "m year-1", 1);
}


void RoutingHydrology_bwatvel::compute(IceModelVec* &output) {

  IceModelVec2Stag *result = new IceModelVec2Stag;
  result->create(grid, "bwatvel", WITHOUT_GHOSTS);
  result->metadata(0) = vars[0];
  result->metadata(1) = vars[1];
  result->write_in_glaciological_units = true;

  model->velocity_staggered(*result);

  output = result;
}

} // end of namespace pism
