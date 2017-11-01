// Copyright (C) 2012-2017 PISM Authors
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

#include <cassert>

#include "Hydrology.hh"
#include "pism/util/Mask.hh"
#include "pism/util/Vars.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/io/PIO.hh"
#include "pism/util/pism_options.hh"
#include "hydrology_diagnostics.hh"
#include "pism/util/MaxTimestep.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/util/IceModelVec2CellType.hh"

namespace pism {
namespace hydrology {

Routing::BoundaryAccounting::BoundaryAccounting() {
  reset();
}

void Routing::BoundaryAccounting::reset() {
  ice_free_land_loss      = 0.0;
  ocean_loss              = 0.0;
  null_strip_loss         = 0.0;
  negative_thickness_gain = 0.0;
}

Routing::BoundaryAccounting Routing::boundary_mass_accounting() const {
  return m_boundary_accounting;
}

Routing::Routing(IceGrid::ConstPtr g)
  : Hydrology(g), m_dx(g->dx()), m_dy(g->dy())
{
  m_stripwidth = m_config->get_double("hydrology.null_strip_width");

  // these variables are also set to zero every time init() is called
  m_boundary_accounting.reset();

  // model state variables; need ghosts
  m_W.create(m_grid, "bwat", WITH_GHOSTS, 1);
  m_W.set_attrs("model_state",
              "thickness of transportable subglacial water layer",
              "m", "");
  m_W.metadata().set_double("valid_min", 0.0);

  // auxiliary variables which NEED ghosts
  m_Wstag.create(m_grid, "W_staggered", WITH_GHOSTS, 1);
  m_Wstag.set_attrs("internal",
                  "cell face-centered (staggered) values of water layer thickness",
                  "m", "");
  m_Wstag.metadata().set_double("valid_min", 0.0);
  m_K.create(m_grid, "K_staggered", WITH_GHOSTS, 1);
  m_K.set_attrs("internal",
                  "cell face-centered (staggered) values of nonlinear conductivity",
                  "", "");
  m_K.metadata().set_double("valid_min", 0.0);
  m_Q.create(m_grid, "advection_flux", WITH_GHOSTS, 1);
  m_Q.set_attrs("internal",
                  "cell face-centered (staggered) components of advective subglacial water flux",
                  "m2 s-1", "");
  m_R.create(m_grid, "potential_workspace", WITH_GHOSTS, 1); // box stencil used
  m_R.set_attrs("internal",
              "work space for modeled subglacial water hydraulic potential",
              "Pa", "");

  // auxiliary variables which do not need ghosts

  m_V.create(m_grid, "water_velocity", WITHOUT_GHOSTS);
  m_V.set_attrs("internal",
              "cell face-centered (staggered) components of water velocity in subglacial water layer",
              "m s-1", "");

  // temporaries during update; do not need ghosts
  m_Wnew.create(m_grid, "Wnew_internal", WITHOUT_GHOSTS);
  m_Wnew.set_attrs("internal",
                 "new thickness of transportable subglacial water layer during update",
                 "m", "");
  m_Wnew.metadata().set_double("valid_min", 0.0);
  m_Wtilnew.create(m_grid, "Wtilnew_internal", WITHOUT_GHOSTS);
  m_Wtilnew.set_attrs("internal",
                    "new thickness of till (subglacial) water layer during update",
                    "m", "");
  m_Wtilnew.metadata().set_double("valid_min", 0.0);

  {
    double alpha = m_config->get_double("hydrology.thickness_power_in_flux");
    if (alpha < 1.0) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                    "alpha = %f < 1 which is not allowed", alpha);
    }
  }
}

Routing::~Routing() {
  // empty
}


void Routing::init() {
  m_log->message(2,
             "* Initializing the routing subglacial hydrology model ...\n");
  // initialize water layer thickness from the context if present,
  //   otherwise from -i file, otherwise with constant value

  options::Real hydrology_null_strip("-hydrology_null_strip",
                                     "set the width, in km, of the strip around the edge"
                                     " of the computational domain in which hydrology is inactivated",
                                     units::convert(m_sys, m_stripwidth, "m", "km"));
  m_stripwidth = units::convert(m_sys, hydrology_null_strip, "km", "m");

  Hydrology::init();

  init_bwat();

  m_boundary_accounting.reset();
}

MaxTimestep Routing::max_timestep_impl(double t) const {
  (void) t;
  return MaxTimestep("routing hydrology");
}


void Routing::init_bwat() {

  // initialize water layer thickness from the context if present, otherwise from -i file,
  // otherwise with constant value

  InputOptions opts = process_input_options(m_grid->com);

  const PetscReal bwatdefault = m_config->get_double("bootstrapping.defaults.bwat");

  switch (opts.type) {
  case INIT_RESTART:
    m_W.read(opts.filename, opts.record);
    break;
  case INIT_BOOTSTRAP:
    m_W.regrid(opts.filename, OPTIONAL, bwatdefault);
    break;
  default:
    m_W.set(bwatdefault);
  }

  // however we initialized it, we could be asked to regrid from file
  regrid("hydrology::Routing", m_W);
}


void Routing::define_model_state_impl(const PIO &output) const {
  Hydrology::define_model_state_impl(output);
  m_W.define(output);
}

void Routing::write_model_state_impl(const PIO &output) const {
  Hydrology::write_model_state_impl(output);
  m_W.write(output);
}

std::map<std::string, Diagnostic::Ptr> Routing::diagnostics_impl() const {
  std::map<std::string, Diagnostic::Ptr> result = {
    {"bwat",       Diagnostic::Ptr(new BasalWaterThickness(this))},
    {"bwatvel",    Diagnostic::Ptr(new BasalWaterVelocity(this))},
    {"bwp",        Diagnostic::Ptr(new BasalWaterPressure(this))},
    {"bwprel",     Diagnostic::Ptr(new RelativeBasalWaterPressure(this))},
    {"effbwp",     Diagnostic::Ptr(new EffectiveBasalWaterPressure(this))},
    {"wallmelt",   Diagnostic::Ptr(new WallMelt(this))},
  };
  return combine(result, Hydrology::diagnostics_impl());
}

std::map<std::string, TSDiagnostic::Ptr> Routing::ts_diagnostics_impl() const {
  std::map<std::string, TSDiagnostic::Ptr> result = {
    {"hydro_ice_free_land_loss",      TSDiagnostic::Ptr(new MCHydrology_ice_free_land_loss(this))},
    {"hydro_ocean_loss",              TSDiagnostic::Ptr(new MCHydrology_ocean_loss(this))},
    {"hydro_negative_thickness_gain", TSDiagnostic::Ptr(new MCHydrology_negative_thickness_gain(this))},
    {"hydro_null_strip_loss",         TSDiagnostic::Ptr(new MCHydrology_null_strip_loss(this))}
  };
  return result;
}


//! Check thk >= 0 and fails with message if not satisfied.
void Routing::check_water_thickness_nonnegative(const IceModelVec2S &waterthk) {
  IceModelVec::AccessList list{&waterthk};

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      if (waterthk(i, j) < 0.0) {
        throw RuntimeError::formatted(PISM_ERROR_LOCATION, "hydrology::Routing: disallowed negative water layer thickness\n"
                                      "waterthk(i, j) = %.6f m at (i, j)=(%d, %d)",
                                      waterthk(i, j), i, j);
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

}


//! Correct the new water thickness based on boundary requirements.  Do mass accounting.
/*!
At ice free locations and ocean locations we require that water thicknesses
(i.e. both the transportable water thickness \f$W\f$ and the till water
thickness \f$W_{til}\f$) be zero at the end of each time step.  Also we require
that any negative water thicknesses be set to zero (i.e. we do projection to
enforce lower bound).  This method does not enforce any upper bounds.

This method should be called once for each thickness field which needs to be
processed.  This method alters the field water_thickness in-place and sums
the boundary adjustments.
 */
void Routing::boundary_mass_changes(IceModelVec2S &water_thickness,
                                    double &icefreelost, double &oceanlost,
                                    double &negativegain, double &nullstriplost) {
  const double fresh_water_density = m_config->get_double("constants.fresh_water.density");

  double
    my_icefreelost  = 0.0,
    my_oceanlost    = 0.0,
    my_negativegain = 0.0;

  const IceModelVec2S        &cell_area = *m_grid->variables().get_2d_scalar("cell_area");
  const IceModelVec2CellType &mask      = *m_grid->variables().get_2d_cell_type("mask");

  IceModelVec::AccessList list{&water_thickness, &cell_area, &mask};

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    const double dmassdz = cell_area(i, j) * fresh_water_density; // kg m-1
    if (water_thickness(i, j) < 0.0) {
      my_negativegain += -water_thickness(i, j) * dmassdz;
      water_thickness(i, j) = 0.0;
    }
    if (mask.ice_free_land(i, j) and (water_thickness(i, j) > 0.0)) {
      my_icefreelost += water_thickness(i, j) * dmassdz;
      water_thickness(i, j) = 0.0;
    }
    if (mask.ocean(i, j) and (water_thickness(i, j) > 0.0)) {
      my_oceanlost += water_thickness(i, j) * dmassdz;
      water_thickness(i, j) = 0.0;
    }
  }

  // make global over all processor domains (i.e. whole glacier/ice sheet)
  icefreelost  = GlobalSum(m_grid->com, my_icefreelost);
  oceanlost    = GlobalSum(m_grid->com, my_oceanlost);
  negativegain = GlobalSum(m_grid->com, my_negativegain);

  if (m_stripwidth <= 0.0) {
    nullstriplost = 0.0;
    return;
  }

  double my_nullstriplost = 0.0;
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    const double dmassdz = cell_area(i, j) * fresh_water_density; // kg m-1
    if (in_null_strip(*m_grid, i, j, m_stripwidth)) {
      my_nullstriplost += water_thickness(i, j) * dmassdz;
      water_thickness(i, j) = 0.0;
    }
  }

  nullstriplost = GlobalSum(m_grid->com, my_nullstriplost);
}


//! Copies the W variable, the modeled transportable water layer thickness.
const IceModelVec2S& Routing::subglacial_water_thickness() const {
  return m_W;
}


//! Returns the (trivial) overburden pressure as the pressure of the transportable water,
//! because this is the model.
const IceModelVec2S& Routing::subglacial_water_pressure() const {
  return m_Pover;
}


//! Get the hydraulic potential from bedrock topography and current state variables.
/*!
  Computes \f$\psi = P + \rho_w g (b + W)\f$ except where floating, where \f$\psi = P_o\f$.
*/
void Routing::compute_hydraulic_potential(const IceModelVec2S &W,
                                          const IceModelVec2S &P,
                                          const IceModelVec2S &P_overburden,
                                          const IceModelVec2S &bed,
                                          const IceModelVec2CellType &mask,
                                          IceModelVec2S &result) const {

  const double
    rho = m_config->get_double("constants.fresh_water.density"),
    g   = m_config->get_double("constants.standard_gravity");

  IceModelVec::AccessList list{&P, &P_overburden, &W, &mask, &bed, &result};

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (mask.ocean(i, j)) {
      result(i, j) = P_overburden(i, j);
    } else {
      result(i, j) = P(i, j) + rho * g * (bed(i, j) + W(i, j));
    }
  }
}


//! Average the regular grid water thickness to values at the center of cell edges.
/*! Uses mask values to avoid averaging using water thickness values from
  either ice-free or floating areas. */
void Routing::water_thickness_staggered(const IceModelVec2S &W,
                                        const IceModelVec2CellType &mask,
                                        IceModelVec2Stag &result) {

  IceModelVec::AccessList list{&mask, &W, &result};

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    // east
    if (mask.grounded_ice(i, j)) {
      if (mask.grounded_ice(i+1, j)) {
        result(i, j, 0) = 0.5 * (W(i, j) + W(i+1, j));
      } else {
        result(i, j, 0) = W(i, j);
      }
    } else {
      if (mask.grounded_ice(i+1, j)) {
        result(i, j, 0) = W(i+1, j);
      } else {
        result(i, j, 0) = 0.0;
      }
    }
    // north
    if (mask.grounded_ice(i, j)) {
      if (mask.grounded_ice(i, j+1)) {
        result(i, j, 1) = 0.5 * (W(i, j) + W(i, j+1));
      } else {
        result(i, j, 1) = W(i, j);
      }
    } else {
      if (mask.grounded_ice(i, j+1)) {
        result(i, j, 1) = W(i, j+1);
      } else {
        result(i, j, 1) = 0.0;
      }
    }
  }

  result.update_ghosts();
}


//! Compute the nonlinear conductivity at the center of cell edges.
/*!
  Computes
  \f[ K = K(W, \nabla P, \nabla b) = k W^{\alpha-1} |\nabla R|^{\beta-2} \f]
  on the staggered grid, where \f$R = P+\rho_w g b\f$.  We denote
  \f$\Pi = |\nabla R|^2\f$ internally; this is computed on a staggered grid
  by a [\ref Mahaffy] -like scheme.  This requires \f$R\f$ to be defined on a box
  stencil of width 1.

  Also returns the maximum over all staggered points of \f$ K W \f$.
*/
void Routing::compute_conductivity(const IceModelVec2Stag &W,
                                   const IceModelVec2S &P,
                                   const IceModelVec2S &bed,
                                   IceModelVec2Stag &result,
                                   double &maxKW) const {
  const double
    k     = m_config->get_double("hydrology.hydraulic_conductivity"),
    alpha = m_config->get_double("hydrology.thickness_power_in_flux"),
    beta  = m_config->get_double("hydrology.gradient_power_in_flux"),
    rg    = (m_config->get_double("constants.standard_gravity") *
             m_config->get_double("constants.fresh_water.density"));

  IceModelVec::AccessList list(result);

  // the following calculation is bypassed if beta == 2.0 exactly; it puts
  // the squared norm of the gradient of the simplified hydrolic potential
  // temporarily in "result"
  if (beta != 2.0) {
    m_R.copy_from(P);  // yes, it updates ghosts
    m_R.add(rg, bed); // R  <-- P + rhow g b
    m_R.update_ghosts();

    list.add(m_R);
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      double dRdx, dRdy;
      dRdx = (m_R(i + 1, j) - m_R(i, j)) / m_dx;
      dRdy = (m_R(i + 1, j + 1) + m_R(i, j + 1) - m_R(i + 1, j - 1) - m_R(i, j - 1)) / (4.0 * m_dy);
      result(i, j, 0) = dRdx * dRdx + dRdy * dRdy;
      dRdx = (m_R(i + 1, j + 1) + m_R(i + 1, j) - m_R(i - 1, j + 1) - m_R(i - 1, j)) / (4.0 * m_dx);
      dRdy = (m_R(i, j + 1) - m_R(i, j)) / m_dy;
      result(i, j, 1) = dRdx * dRdx + dRdy * dRdy;
    }
  }

  double betapow = (beta-2.0)/2.0, mymaxKW = 0.0;

  list.add(W);

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    for (int o = 0; o < 2; ++o) {
      double Ktmp = k * pow(W(i, j, o), alpha-1.0);
      if (beta < 2.0) {
        // regularize negative power |\grad psi|^{beta-2} by adding eps because
        //   large head gradient might be 10^7 Pa per 10^4 m or 10^3 Pa/m
        const double eps = 1.0;   // Pa m-1
        result(i, j, o) = Ktmp * pow(result(i, j, o) + eps * eps, betapow);
      } else if (beta > 2.0) {
        result(i, j, o) = Ktmp * pow(result(i, j, o), betapow);
      } else { // beta == 2.0
        result(i, j, o) = Ktmp;
      }
      mymaxKW = std::max(mymaxKW, result(i, j, o) * W(i, j, o));
    }
  }

  maxKW = GlobalMax(m_grid->com, mymaxKW);

  result.update_ghosts();
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
void wall_melt(const Routing &model, IceModelVec2S &result) {

  IceGrid::ConstPtr grid = result.get_grid();

  Config::ConstPtr config = grid->ctx()->config();

  const double
    k     = config->get_double("hydrology.hydraulic_conductivity"),
    L     = config->get_double("constants.fresh_water.latent_heat_of_fusion"),
    alpha = config->get_double("hydrology.thickness_power_in_flux"),
    beta  = config->get_double("hydrology.gradient_power_in_flux"),
    rhow  = config->get_double("constants.standard_gravity"),
    g     = config->get_double("constants.fresh_water.density"),
    rg    = rhow * g,
    CC    = k / (L * rhow);

  const IceModelVec2S *bed = grid->variables().get_2d_scalar("bedrock_altitude");

  // FIXME:  could be scaled with overall factor hydrology_coefficient_wall_melt ?
  if (alpha < 1.0) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "alpha = %f < 1 which is not allowed", alpha);
  }

  IceModelVec2S R;
  R.create(grid, "R", WITH_GHOSTS);

  R.copy_from(model.subglacial_water_pressure());  // yes, it updates ghosts
  R.add(rg, *bed); // R  <-- P + rhow g b

  IceModelVec2S W;
  W.create(grid, "W", WITH_GHOSTS);
  W.copy_from(model.subglacial_water_thickness());

  IceModelVec::AccessList list{&R, &W, &result};

  double dx = grid->dx();
  double dy = grid->dy();

  for (Points p(*grid); p; p.next()) {
    const int i = p.i(), j = p.j();
    double dRdx, dRdy;

    if (W(i, j) > 0.0) {
      dRdx = 0.0;
      if (W(i + 1, j) > 0.0) {
        dRdx = (R(i + 1, j) - R(i, j)) / (2.0 * dx);
      }
      if (W(i - 1, j) > 0.0) {
        dRdx += (R(i, j) - R(i - 1, j)) / (2.0 * dx);
      }
      dRdy = 0.0;
      if (W(i, j + 1) > 0.0) {
        dRdy = (R(i, j + 1) - R(i, j)) / (2.0 * dy);
      }
      if (W(i, j - 1) > 0.0) {
        dRdy += (R(i, j) - R(i, j - 1)) / (2.0 * dy);
      }
      result(i, j) = CC * pow(W(i, j), alpha) * pow(dRdx * dRdx + dRdy * dRdy, beta/2.0);
    } else {
      result(i, j) = 0.0;
    }
  }
}


//! Get the advection velocity V at the center of cell edges.
/*!
Computes the advection velocity @f$\mathbf{V}@f$ on the staggered
(edge-centered) grid.  If V = (u, v) in components then we have
<code> result(i, j, 0) = u(i+1/2, j) </code> and
<code> result(i, j, 1) = v(i, j+1/2) </code>

The advection velocity is given by the formula
  @f[ \mathbf{V} = - K \left(\nabla P + \rho_w g \nabla b\right) @f]
where @f$\mathbf{V}@f$ is the water velocity, @f$P@f$ is the water
pressure, and @f$b@f$ is the bedrock elevation.

If the corresponding staggered grid value of the water thickness is zero then
that component of V is set to zero.  This does not change the flux value (which
would be zero anyway) but it does provide the correct max velocity in the
CFL calculation.  We assume Wstag and K are up-to-date.  We assume P and b
have valid ghosts.
 */
void Routing::compute_velocity(const IceModelVec2Stag &W,
                               const IceModelVec2S &P,
                               const IceModelVec2S &bed,
                               const IceModelVec2Stag &K,
                               IceModelVec2Stag &result) const {
  const double  rg = (m_config->get_double("constants.standard_gravity") *
                      m_config->get_double("constants.fresh_water.density"));
  double dbdx, dbdy, dPdx, dPdy;

  IceModelVec2S &pressure = m_R;
  pressure.copy_from(P);  // yes, it updates ghosts

  IceModelVec::AccessList list{&pressure, &W, &K, &bed, &result};

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (W(i, j, 0) > 0.0) {
      dPdx = (pressure(i + 1, j) - pressure(i, j)) / m_dx;
      dbdx = (bed(i + 1, j) - bed(i, j)) / m_dx;
      result(i, j, 0) =  - K(i, j, 0) * (dPdx + rg * dbdx);
    } else {
      result(i, j, 0) = 0.0;
    }

    if (W(i, j, 1) > 0.0) {
      dPdy = (pressure(i, j + 1) - pressure(i, j)) / m_dy;
      dbdy = (bed(i, j + 1) - bed(i, j)) / m_dy;
      result(i, j, 1) =  - K(i, j, 1) * (dPdy + rg * dbdy);
    } else {
      result(i, j, 1) = 0.0;
    }

    if (in_null_strip(*m_grid, i, j, m_stripwidth) or
        in_null_strip(*m_grid, i + 1, j, m_stripwidth)) {
      result(i, j, 0) = 0.0;
    }

    if (in_null_strip(*m_grid, i, j, m_stripwidth) or
        in_null_strip(*m_grid, i, j + 1, m_stripwidth)) {
      result(i, j, 1) = 0.0;
    }
  }
}


//! Compute Q = V W at edge-centers (staggered grid) by first-order upwinding.
/*!
The field W must have valid ghost values, but V does not need them.

FIXME:  This could be re-implemented using the Koren (1993) flux-limiter.
 */
void Routing::advective_fluxes(const IceModelVec2Stag &V,
                               const IceModelVec2S &W,
                               IceModelVec2Stag &result) const {
  IceModelVec::AccessList list{&W, &V, &result};

  assert(W.get_stencil_width() >= 1);

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    result(i, j, 0) = V(i, j, 0) * (V(i, j, 0) >= 0.0 ? W(i, j) :  W(i+1, j));
    result(i, j, 1) = V(i, j, 1) * (V(i, j, 1) >= 0.0 ? W(i, j) :  W(i, j+1));
  }

  result.update_ghosts();
}


//! Compute the adaptive time step for evolution of W.
void Routing::adaptive_for_W_evolution(double t_current, double t_end, double maxKW,
                                       double &dt_result,
                                       double &maxV_result, double &maxD_result,
                                       double &dtCFL_result, double &dtDIFFW_result) {
  const double
    dtmax = m_config->get_double("hydrology.maximum_time_step", "seconds"),
    rg    = m_config->get_double("constants.standard_gravity") * m_config->get_double("constants.fresh_water.density");

  // V could be zero if P is constant and bed is flat
  std::vector<double> tmp = m_V.absmaxcomponents();
  maxV_result = sqrt(tmp[0]*tmp[0] + tmp[1]*tmp[1]);
  maxD_result = rg * maxKW;
  dtCFL_result = 0.5 / (tmp[0]/m_dx + tmp[1]/m_dy); // FIXME: is regularization needed?
  dtDIFFW_result = 1.0/(m_dx*m_dx) + 1.0/(m_dy*m_dy);
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

Compare hydrology::NullTransport::update().  The current code is not quite "code
duplication" because the code here: (1) computes `Wtilnew` instead of updating
`Wtil` in place; (2) uses time steps determined by the rest of the
hydrology::Routing model; (3) does not check mask because the boundary_mass_changes()
call addresses that.  Otherwise this is the same physical model with the
same configurable parameters.
 */
void Routing::raw_update_Wtil(double hdt) {
  const double
    tillwat_max = m_config->get_double("hydrology.tillwat_max"),
    C           = m_config->get_double("hydrology.tillwat_decay_rate");

  IceModelVec::AccessList list{&m_Wtil, &m_Wtilnew, &m_total_input};

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    m_Wtilnew(i, j) = m_Wtil(i, j) + hdt * (m_total_input(i, j) - C);
    m_Wtilnew(i, j) = std::min(std::max(0.0, m_Wtilnew(i, j)), tillwat_max);
  }
}


//! The computation of Wnew, called by update().
void Routing::raw_update_W(double hdt) {
  const double
    wux = 1.0 / (m_dx * m_dx),
    wuy = 1.0 / (m_dy * m_dy),
    rg  = m_config->get_double("constants.standard_gravity") * m_config->get_double("constants.fresh_water.density");

  IceModelVec::AccessList list{&m_W, &m_Wtil, &m_Wtilnew, &m_Wstag, &m_K, &m_Q,
      &m_total_input, &m_Wnew};

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    const double divadflux =
      (m_Q(i, j, 0) - m_Q(i - 1, j, 0)) / m_dx +
      (m_Q(i, j, 1) - m_Q(i, j - 1, 1)) / m_dy;

    const double
      De = rg * m_K(i,     j,     0) * m_Wstag(i,     j,     0),
      Dw = rg * m_K(i - 1, j,     0) * m_Wstag(i - 1, j,     0),
      Dn = rg * m_K(i,     j,     1) * m_Wstag(i,     j,     1),
      Ds = rg * m_K(i,     j - 1, 1) * m_Wstag(i,     j - 1, 1);

    const double diffW =
      wux * (De * (m_W(i + 1, j) - m_W(i, j)) - Dw * (m_W(i, j) - m_W(i - 1, j))) +
      wuy * (Dn * (m_W(i, j + 1) - m_W(i, j)) - Ds * (m_W(i, j) - m_W(i, j - 1)));

    m_Wnew(i, j) = m_W(i, j) - m_Wtilnew(i, j) + m_Wtil(i, j) + hdt * ( - divadflux + diffW + m_total_input(i, j));
  }
}


//! Update the model state variables W and Wtil by applying the subglacial hydrology model equations.
/*!
Runs the hydrology model from time icet to time icet + icedt.  Here [icet, icedt]
is generally on the order of months to years.  This hydrology model will take its
own shorter time steps, perhaps hours to weeks.

To update W = `bwat` we call raw_update_W(), and to update Wtil = `tillwat` we
call raw_update_Wtil().
 */
void Routing::update_impl(double icet, double icedt) {

  // if asked for the identical time interval versus last time, then
  //   do nothing; otherwise assume that [my_t, my_t+my_dt] is the time
  //   interval on which we are solving
  if ((fabs(icet - m_t) < 1e-12) && (fabs(icedt - m_dt) < 1e-12)) {
    return;
  }
  // update Component times: t = current time, t+dt = target time
  m_t = icet;
  m_dt = icedt;

  m_boundary_accounting.reset();

  if (m_config->get_double("hydrology.tillwat_max") < 0.0) {
    throw RuntimeError(PISM_ERROR_LOCATION, "hydrology::Routing: hydrology.tillwat_max is negative.\n"
                       "This is not allowed.");
  }

  // make sure W has valid ghosts before starting hydrology steps
  m_W.update_ghosts();

  double ht = m_t, hdt = 0.0, // hydrology model time and time step
    maxKW = 0.0, maxV = 0.0, maxD = 0.0, dtCFL = 0.0, dtDIFFW = 0.0;
  double icefreelost = 0.0, oceanlost = 0.0, negativegain = 0.0, nullstriplost = 0.0,
    delta_icefree = 0.0, delta_ocean = 0.0, delta_neggain = 0.0, delta_nullstrip = 0.0;

  const IceModelVec2CellType &mask = *m_grid->variables().get_2d_cell_type("mask");
  const IceModelVec2S &bed = *m_grid->variables().get_2d_scalar("bedrock_altitude");

  unsigned int step_counter = 0;
  while (ht < m_t + m_dt) {
    step_counter++;

#if (PISM_DEBUG==1)
    check_water_thickness_nonnegative(m_W);
    check_Wtil_bounds();
#endif

    water_thickness_staggered(m_W,
                              mask,
                              m_Wstag);

    compute_conductivity(m_Wstag,
                         subglacial_water_pressure(),
                         bed,
                         m_K, maxKW);

    compute_velocity(m_Wstag,
                     subglacial_water_pressure(),
                     bed,
                     m_K,
                     m_V);

    // to get Q, W needs valid ghosts
    advective_fluxes(m_V, m_W, m_Q);

    adaptive_for_W_evolution(ht, m_t+m_dt, maxKW,
                             hdt, maxV, maxD, dtCFL, dtDIFFW);

    if ((m_inputtobed != NULL) || (step_counter==1)) {
      get_input_rate(ht, hdt, m_total_input);
    }

    // update Wtilnew from Wtil
    raw_update_Wtil(hdt);

    boundary_mass_changes(m_Wtilnew, delta_icefree, delta_ocean,
                          delta_neggain, delta_nullstrip);
    icefreelost  += delta_icefree;
    oceanlost    += delta_ocean;
    negativegain += delta_neggain;
    nullstriplost+= delta_nullstrip;

    // update Wnew from W, Wtil, Wtilnew, Wstag, Q, total_input
    raw_update_W(hdt);
    boundary_mass_changes(m_Wnew, delta_icefree, delta_ocean,
                          delta_neggain, delta_nullstrip);
    icefreelost  += delta_icefree;
    oceanlost    += delta_ocean;
    negativegain += delta_neggain;
    nullstriplost+= delta_nullstrip;

    // transfer new into old
    m_Wnew.update_ghosts(m_W);
    m_Wtil.copy_from(m_Wtilnew);

    ht += hdt;
  } // end of hydrology model time-stepping loop

  m_log->message(2,
             "  'routing' hydrology took %d hydrology sub-steps with average dt = %.6f years\n",
             step_counter, units::convert(m_sys, m_dt/step_counter, "seconds", "years"));

  m_log->message(3,
             "  (hydrology info: dt = %.2f s, max |V| = %.2e m s-1, max D = %.2e m^2 s-1)\n",
             m_dt/step_counter, maxV, maxD);

  m_boundary_accounting.ice_free_land_loss      += icefreelost;
  m_boundary_accounting.ocean_loss              += oceanlost;
  m_boundary_accounting.negative_thickness_gain += negativegain;
  m_boundary_accounting.null_strip_loss         += nullstriplost;
}

const IceModelVec2Stag& Routing::velocity_staggered() const {
  return m_V;
}

BasalWaterVelocity::BasalWaterVelocity(const Routing *m)
  : Diag<Routing>(m) {

  // set metadata:
  m_dof = 2;
  m_vars = {SpatialVariableMetadata(m_sys, "bwatvel[0]"),
            SpatialVariableMetadata(m_sys, "bwatvel[1]")};

  set_attrs("velocity of water in subglacial layer, i-offset", "",
            "m s-1", "m year-1", 0);
  set_attrs("velocity of water in subglacial layer, j-offset", "",
            "m s-1", "m year-1", 1);
}

IceModelVec::Ptr BasalWaterVelocity::compute_impl() const {
  IceModelVec2Stag::Ptr result(new IceModelVec2Stag);
  result->create(m_grid, "bwatvel", WITHOUT_GHOSTS);
  result->metadata(0) = m_vars[0];
  result->metadata(1) = m_vars[1];

  result->copy_from(model->velocity_staggered());

  return result;
}

} // end of namespace hydrology
} // end of namespace pism
