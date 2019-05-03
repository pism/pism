// Copyright (C) 2012-2018 PISM Authors
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

#include "Routing.hh"
#include "pism/util/IceModelVec2CellType.hh"
#include "pism/util/Mask.hh"
#include "pism/util/MaxTimestep.hh"

#include "pism/util/error_handling.hh"

#include "pism/util/pism_options.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/util/Vars.hh"

namespace pism {
namespace hydrology {

namespace diagnostics {

//! \brief Reports the pressure of the transportable water in the subglacial layer.
class BasalWaterPressure : public Diag<Routing>
{
public:
  BasalWaterPressure(const Routing *m)
    : Diag<Routing>(m) {
    m_vars = {SpatialVariableMetadata(m_sys, "bwp")};
    set_attrs("pressure of transportable water in subglacial layer", "", "Pa", "Pa", 0);
  }

protected:
  virtual IceModelVec::Ptr compute_impl() const {
    IceModelVec2S::Ptr result(new IceModelVec2S(m_grid, "bwp", WITHOUT_GHOSTS));
    result->metadata() = m_vars[0];
    result->copy_from(model->subglacial_water_pressure());
    return result;
  }
};


//! \brief Reports the pressure of the transportable water in the subglacial layer as a
//! fraction of the overburden pressure.
class RelativeBasalWaterPressure : public Diag<Routing>
{
public:
  RelativeBasalWaterPressure(const Routing *m)
    : Diag<Routing>(m) {
    m_vars = {SpatialVariableMetadata(m_sys, "bwprel")};
    set_attrs("pressure of transportable water in subglacial layer"
              " as fraction of the overburden pressure", "",
              "", "", 0);
    m_vars[0].set_double("_FillValue", m_fill_value);
  }

protected:
  virtual IceModelVec::Ptr compute_impl() const {
    double fill_value = m_fill_value;

    IceModelVec2S::Ptr result(new IceModelVec2S(m_grid, "bwprel", WITHOUT_GHOSTS));
    result->metadata(0) = m_vars[0];

    const IceModelVec2S
      &P  = model->subglacial_water_pressure(),
      &Po = model->overburden_pressure();

    IceModelVec::AccessList list{result.get(), &Po, &P};
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      if (Po(i,j) > 0.0) {
        (*result)(i,j) = P(i, j) / Po(i,j);
      } else {
        (*result)(i,j) = fill_value;
      }
    }

    return result;
  }
};


//! \brief Reports the effective pressure of the transportable water in the subglacial
//! layer, that is, the overburden pressure minus the pressure.
class EffectiveBasalWaterPressure : public Diag<Routing>
{
public:
  EffectiveBasalWaterPressure(const Routing *m)
    : Diag<Routing>(m) {
    m_vars = {SpatialVariableMetadata(m_sys, "effbwp")};
    set_attrs("effective pressure of transportable water in subglacial layer"
              " (overburden pressure minus water pressure)",
              "", "Pa", "Pa", 0);
  }

protected:
  virtual IceModelVec::Ptr compute_impl() const {

    IceModelVec2S::Ptr result(new IceModelVec2S(m_grid, "effbwp", WITHOUT_GHOSTS));
    result->metadata() = m_vars[0];

    const IceModelVec2S
      &P  = model->subglacial_water_pressure(),
      &Po = model->overburden_pressure();

    IceModelVec::AccessList list{&Po, &P, result.get()};

    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      (*result)(i, j) = Po(i, j) - P(i, j);
    }

    return result;
  }
};


//! \brief Report the wall melt rate from dissipation of the potential energy of the
//! transportable water.
class WallMelt : public Diag<Routing>
{
public:
  WallMelt(const Routing *m)
    : Diag<Routing>(m) {
    m_vars = {SpatialVariableMetadata(m_sys, "wallmelt")};
    set_attrs("wall melt into subglacial hydrology layer"
              " from (turbulent) dissipation of energy in transportable water",
              "", "m s-1", "m year-1", 0);
  }

protected:
  virtual IceModelVec::Ptr compute_impl() const {
    IceModelVec2S::Ptr result(new IceModelVec2S(m_grid, "wallmelt", WITHOUT_GHOSTS));
    result->metadata() = m_vars[0];

    const IceModelVec2S &bed_elevation = *m_grid->variables().get_2d_scalar("bedrock_altitude");

    wall_melt(*model, bed_elevation, *result);
    return result;
  }
};

//! @brief Diagnostically reports the staggered-grid components of the velocity of the
//! water in the subglacial layer.
class BasalWaterVelocity : public Diag<Routing>
{
public:
  BasalWaterVelocity(const Routing *m)
    : Diag<Routing>(m) {

    // set metadata:
    m_vars = {SpatialVariableMetadata(m_sys, "bwatvel[0]"),
              SpatialVariableMetadata(m_sys, "bwatvel[1]")};

    set_attrs("velocity of water in subglacial layer, i-offset", "",
              "m s-1", "m year-1", 0);
    set_attrs("velocity of water in subglacial layer, j-offset", "",
              "m s-1", "m year-1", 1);
  }
protected:
  virtual IceModelVec::Ptr compute_impl() const {
    IceModelVec2Stag::Ptr result(new IceModelVec2Stag);
    result->create(m_grid, "bwatvel", WITHOUT_GHOSTS);
    result->metadata(0) = m_vars[0];
    result->metadata(1) = m_vars[1];

    result->copy_from(model->velocity_staggered());

    return result;
  }
};

} // end of namespace diagnostics

Routing::Routing(IceGrid::ConstPtr g)
  : Hydrology(g), m_dx(g->dx()), m_dy(g->dy()) {

  m_W.metadata().set_string("pism_intent", "model_state");

  m_rg    = (m_config->get_double("constants.fresh_water.density") *
             m_config->get_double("constants.standard_gravity"));

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
                "cell face-centered (staggered) components of water velocity"
                " in subglacial water layer",
                "m s-1", "");

  // temporaries during update; do not need ghosts
  m_Wnew.create(m_grid, "W_new", WITHOUT_GHOSTS);
  m_Wnew.set_attrs("internal",
                   "new thickness of transportable subglacial water layer during update",
                   "m", "");
  m_Wnew.metadata().set_double("valid_min", 0.0);

  m_Wtillnew.create(m_grid, "Wtill_new", WITHOUT_GHOSTS);
  m_Wtillnew.set_attrs("internal",
                       "new thickness of till (subglacial) water layer during update",
                       "m", "");
  m_Wtillnew.metadata().set_double("valid_min", 0.0);

  {
    double alpha = m_config->get_double("hydrology.thickness_power_in_flux");
    if (alpha < 1.0) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                    "alpha = %f < 1 which is not allowed", alpha);
    }

    if (m_config->get_double("hydrology.tillwat_max") < 0.0) {
      throw RuntimeError(PISM_ERROR_LOCATION,
                         "hydrology::Routing: hydrology.tillwat_max is negative.\n"
                         "This is not allowed.");
    }
  }
}

Routing::~Routing() {
  // empty
}

void Routing::initialization_message() const {
  m_log->message(2,
                 "* Initializing the routing subglacial hydrology model ...\n");
}

void Routing::restart_impl(const PIO &input_file, int record) {
  Hydrology::restart_impl(input_file, record);

  m_W.read(input_file, record);

  regrid("Hydrology", m_W);
}

void Routing::bootstrap_impl(const PIO &input_file,
                             const IceModelVec2S &ice_thickness) {
  Hydrology::bootstrap_impl(input_file, ice_thickness);

  double bwat_default = m_config->get_double("bootstrapping.defaults.bwat");
  m_W.regrid(input_file, OPTIONAL, bwat_default);

  regrid("Hydrology", m_W);
}

void Routing::initialize_impl(const IceModelVec2S &W_till,
                              const IceModelVec2S &W,
                              const IceModelVec2S &P) {
  Hydrology::initialize_impl(W_till, W, P);

  m_W.copy_from(W);
}

void Routing::define_model_state_impl(const PIO &output) const {
  Hydrology::define_model_state_impl(output);
  m_W.define(output);
}

void Routing::write_model_state_impl(const PIO &output) const {
  Hydrology::write_model_state_impl(output);
  m_W.write(output);
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

  IceModelVec::AccessList list{&P, &P_overburden, &W, &mask, &bed, &result};

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (mask.ocean(i, j)) {
      result(i, j) = P_overburden(i, j);
    } else {
      result(i, j) = P(i, j) + m_rg * (bed(i, j) + W(i, j));
    }
  }
}


//! Average the regular grid water thickness to values at the center of cell edges.
/*! Uses mask values to avoid averaging using water thickness values from
  either ice-free or floating areas. */
void Routing::water_thickness_staggered(const IceModelVec2S &W,
                                        const IceModelVec2CellType &mask,
                                        IceModelVec2Stag &result) {

  IceModelVec::AccessList list{ &mask, &W, &result };

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    // east
    if (mask.grounded_ice(i, j)) {
      if (mask.grounded_ice(i + 1, j)) {
        result(i, j, 0) = 0.5 * (W(i, j) + W(i + 1, j));
      } else {
        result(i, j, 0) = W(i, j);
      }
    } else {
      if (mask.grounded_ice(i + 1, j)) {
        result(i, j, 0) = W(i + 1, j);
      } else {
        result(i, j, 0) = 0.0;
      }
    }
    // north
    if (mask.grounded_ice(i, j)) {
      if (mask.grounded_ice(i, j + 1)) {
        result(i, j, 1) = 0.5 * (W(i, j) + W(i, j + 1));
      } else {
        result(i, j, 1) = W(i, j);
      }
    } else {
      if (mask.grounded_ice(i, j + 1)) {
        result(i, j, 1) = W(i, j + 1);
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

  \f[ \Pi = |\nabla R|^2 \f]

  internally; this is computed on a staggered grid by a Mahaffy-like ([@ref Mahaffy])
  scheme. This requires \f$R\f$ to be defined on a box stencil of width 1.

  Also returns the maximum over all staggered points of \f$ K W \f$.
*/
void Routing::compute_conductivity(const IceModelVec2Stag &W,
                                   const IceModelVec2S &P,
                                   const IceModelVec2S &bed_elevation,
                                   IceModelVec2Stag &result,
                                   double &KW_max) const {
  const double
    k     = m_config->get_double("hydrology.hydraulic_conductivity"),
    alpha = m_config->get_double("hydrology.thickness_power_in_flux"),
    beta  = m_config->get_double("hydrology.gradient_power_in_flux"),
    betapow = (beta - 2.0) / 2.0;

  IceModelVec::AccessList list({&result, &W});

  KW_max = 0.0;

  if (beta != 2.0) {
    // Put the squared norm of the gradient of the simplified hydrolic potential (Pi) in
    // "result"
    //
    // FIXME: we don't need to re-compute this during every hydrology time step: the
    // simplified hydrolic potential does not depend on the water amount and can be
    // computed *once* in update_impl(), before entering the time-stepping loop
    {
      // R  <-- P + rhow g b
      P.add(m_rg, bed_elevation, m_R);  // yes, it updates ghosts

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

    // We regularize negative power |\grad psi|^{beta-2} by adding eps because large
    // head gradient might be 10^7 Pa per 10^4 m or 10^3 Pa/m.
    const double eps = beta < 2.0 ? 1.0 : 0.0;

    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      for (int o = 0; o < 2; ++o) {
        const double Pi = result(i, j, o);

        // FIXME: same as Pi above: we don't need to re-compute this each time we make a
        // short hydrology time step
        double B = pow(Pi + eps * eps, betapow);

        result(i, j, o) = k * pow(W(i, j, o), alpha - 1.0) * B;

        KW_max = std::max(KW_max, result(i, j, o) * W(i, j, o));
      }
    }
  } else {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      for (int o = 0; o < 2; ++o) {
        result(i, j, o) = k * pow(W(i, j, o), alpha - 1.0);

        KW_max = std::max(KW_max, result(i, j, o) * W(i, j, o));
      }
    }
  }

  KW_max = GlobalMax(m_grid->com, KW_max);

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
void wall_melt(const Routing &model,
               const IceModelVec2S &bed_elevation,
               IceModelVec2S &result) {

  IceGrid::ConstPtr grid = result.grid();

  Config::ConstPtr config = grid->ctx()->config();

  const double
    k     = config->get_double("hydrology.hydraulic_conductivity"),
    L     = config->get_double("constants.fresh_water.latent_heat_of_fusion"),
    alpha = config->get_double("hydrology.thickness_power_in_flux"),
    beta  = config->get_double("hydrology.gradient_power_in_flux"),
    g     = config->get_double("constants.standard_gravity"),
    rhow  = config->get_double("constants.fresh_water.density"),
    rg    = rhow * g,
    CC    = k / (L * rhow);

  // FIXME:  could be scaled with overall factor hydrology_coefficient_wall_melt ?
  if (alpha < 1.0) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "alpha = %f < 1 which is not allowed", alpha);
  }

  IceModelVec2S R;
  R.create(grid, "R", WITH_GHOSTS);

  // R  <-- P + rhow g b
  model.subglacial_water_pressure().add(rg, bed_elevation, R);
  // yes, it updates ghosts

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

  If the corresponding staggered grid value of the water thickness is zero then that
  component of V is set to zero. This does not change the flux value (which would be zero
  anyway) but it does provide the correct max velocity in the CFL calculation. We assume
  bed has valid ghosts.
*/
void Routing::compute_velocity(const IceModelVec2Stag &W,
                               const IceModelVec2S &P,
                               const IceModelVec2S &bed,
                               const IceModelVec2Stag &K,
                               const IceModelVec2Int *no_model_mask,
                               IceModelVec2Stag &result) const {
  double dbdx, dbdy, dPdx, dPdy;

  IceModelVec2S &pressure = m_R;
  pressure.copy_from(P);  // yes, it updates ghosts

  IceModelVec::AccessList list{&pressure, &W, &K, &bed, &result};

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (W(i, j, 0) > 0.0) {
      dPdx = (pressure(i + 1, j) - pressure(i, j)) / m_dx;
      dbdx = (bed(i + 1, j) - bed(i, j)) / m_dx;
      result(i, j, 0) =  - K(i, j, 0) * (dPdx + m_rg * dbdx);
    } else {
      result(i, j, 0) = 0.0;
    }

    if (W(i, j, 1) > 0.0) {
      dPdy = (pressure(i, j + 1) - pressure(i, j)) / m_dy;
      dbdy = (bed(i, j + 1) - bed(i, j)) / m_dy;
      result(i, j, 1) =  - K(i, j, 1) * (dPdy + m_rg * dbdy);
    } else {
      result(i, j, 1) = 0.0;
    }
  }

  if (no_model_mask) {
    const IceModelVec2Int &M = *no_model_mask;
    list.add(M);

    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      if (M(i, j) or M(i + 1, j)) {
        result(i, j, 0) = 0.0;
      }

      if (M(i, j) or M(i, j + 1)) {
        result(i, j, 1) = 0.0;
      }
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

  assert(W.stencil_width() >= 1);

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    result(i, j, 0) = V(i, j, 0) * (V(i, j, 0) >= 0.0 ? W(i, j) :  W(i + 1, j));
    result(i, j, 1) = V(i, j, 1) * (V(i, j, 1) >= 0.0 ? W(i, j) :  W(i, j + 1));
  }

  result.update_ghosts();
}

/*!
 * See equation (51) in Bueler and van Pelt.
 */
double Routing::max_timestep_W_diff(double KW_max) const {
  double D_max = m_rg * KW_max;
  double result = 1.0/(m_dx*m_dx) + 1.0/(m_dy*m_dy);
  return 0.25 / (D_max * result);
}

/*!
 * See equation (50) in Bueler and van Pelt.
 */
double Routing::max_timestep_W_cfl() const {
  // V could be zero if P is constant and bed is flat
  std::vector<double> tmp = m_V.absmaxcomponents();

  return 0.5 / (tmp[0]/m_dx + tmp[1]/m_dy); // FIXME: is regularization needed?
}


//! The computation of Wtillnew, called by update().
/*!
  Does a step of the trivial integration
  \f[ \frac{\partial W_{till}}{\partial t} = \frac{m}{\rho_w} - C\f]

  where \f$C=\f$`hydrology_tillwat_decay_rate`.  Enforces bounds
  \f$0 \le W_{till} \le W_{till}^{max}\f$ where the upper bound is
  `hydrology_tillwat_max`.  Here \f$m/\rho_w\f$ is `total_input`.

  Compare hydrology::NullTransport::update_impl().

  The current code is not quite "code duplication" because the code here:

  1. computes `Wtill_new` instead of updating `Wtill` in place;
  2. uses time steps determined by the rest of the hydrology::Routing model;
  3. does not check mask because the enforce_bounds() call addresses that.

  Otherwise this is the same physical model with the same configurable parameters.
*/
void Routing::update_Wtill(double dt,
                           const IceModelVec2S &Wtill,
                           const IceModelVec2S &input_rate,
                           IceModelVec2S &Wtill_new) {
  const double
    tillwat_max = m_config->get_double("hydrology.tillwat_max"),
    C           = m_config->get_double("hydrology.tillwat_decay_rate", "m / second");

  IceModelVec::AccessList list{&Wtill, &Wtill_new, &input_rate};

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    Wtill_new(i, j) = clip(Wtill(i, j) + dt * (input_rate(i, j) - C),
                           0, tillwat_max);
  }
}

void Routing::W_change_due_to_flow(double dt,
                                   const IceModelVec2S    &W,
                                   const IceModelVec2Stag &Wstag,
                                   const IceModelVec2Stag &K,
                                   const IceModelVec2Stag &Q,
                                   IceModelVec2S &result) {
  const double
    wux = 1.0 / (m_dx * m_dx),
    wuy = 1.0 / (m_dy * m_dy);

  IceModelVec::AccessList list{&W, &Wstag, &K, &Q, &result};

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    auto q = Q.star(i, j);
    const double divQ = (q.e - q.w) / m_dx + (q.n - q.s) / m_dy;

    auto k  = K.star(i, j);
    auto ws = Wstag.star(i, j);

    const double
      De = m_rg * k.e * ws.e,
      Dw = m_rg * k.w * ws.w,
      Dn = m_rg * k.n * ws.n,
      Ds = m_rg * k.s * ws.s;

    auto w = W.star(i, j);
    const double diffW = (wux * (De * (w.e - w.ij) - Dw * (w.ij - w.w)) +
                          wuy * (Dn * (w.n - w.ij) - Ds * (w.ij - w.s)));

    result(i, j) = dt * (- divQ + diffW);
  }
}


//! The computation of Wnew, called by update().
void Routing::update_W(double dt,
                       const IceModelVec2S    &input_rate,
                       const IceModelVec2S    &W,
                       const IceModelVec2Stag &Wstag,
                       const IceModelVec2S    &Wtill,
                       const IceModelVec2S    &Wtill_new,
                       const IceModelVec2Stag &K,
                       const IceModelVec2Stag &Q,
                       IceModelVec2S &W_new) {

  W_change_due_to_flow(dt, W, Wstag, K, Q, m_flow_change_incremental);

  IceModelVec::AccessList list{&W, &Wtill, &Wtill_new, &input_rate,
      &m_flow_change_incremental, &W_new};

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    double Wtill_change = Wtill_new(i, j) - Wtill(i, j);
    W_new(i, j) = W(i, j) - Wtill_change + dt * input_rate(i, j) + m_flow_change_incremental(i, j);
  }

  m_flow_change.add(1.0, m_flow_change_incremental);
  m_input_change.add(dt, input_rate);
}


//! Update the model state variables W and Wtill by applying the subglacial hydrology model equations.
/*!
  Runs the hydrology model from time t to time t + dt.  Here [t, dt]
  is generally on the order of months to years.  This hydrology model will take its
  own shorter time steps, perhaps hours to weeks.

  To update W = `bwat` we call update_W(), and to update Wtill = `tillwat` we
  call update_Wtill().
*/
void Routing::update_impl(double t, double dt, const Inputs& inputs) {

  double
    ht  = t,
    hdt = 0.0;

  const double
    t_final = t + dt,
    dt_max  = m_config->get_double("hydrology.maximum_time_step", "seconds");

  // make sure W has valid ghosts before starting hydrology steps
  m_W.update_ghosts();

  unsigned int step_counter = 0;
  for (; ht < t_final; ht += hdt) {
    step_counter++;

#if (PISM_DEBUG==1)
    double huge_number = 1e6;
    check_bounds(m_W, huge_number);

    check_bounds(m_Wtill, m_config->get_double("hydrology.tillwat_max"));
#endif

    water_thickness_staggered(m_W,
                              *inputs.cell_type,
                              m_Wstag);

    double maxKW = 0.0;
    compute_conductivity(m_Wstag,
                         subglacial_water_pressure(),
                         *inputs.bed_elevation,
                         m_K, maxKW);

    compute_velocity(m_Wstag,
                     subglacial_water_pressure(),
                     *inputs.bed_elevation,
                     m_K,
                     inputs.no_model_mask,
                     m_V);

    // to get Q, W needs valid ghosts
    advective_fluxes(m_V, m_W, m_Q);

    {
      const double
        dt_cfl    = max_timestep_W_cfl(),
        dt_diff_w = max_timestep_W_diff(maxKW);

      hdt = std::min(t_final - ht, dt_max);
      hdt = std::min(hdt, dt_cfl);
      hdt = std::min(hdt, dt_diff_w);
    }

    m_log->message(3, "  hydrology step %05d, dt = %f s\n", step_counter, hdt);

    // update Wtillnew from Wtill and input_rate
    update_Wtill(hdt,
                 m_Wtill,
                 m_input_rate,
                 m_Wtillnew);
    // remove water in ice-free areas and account for changes
    enforce_bounds(*inputs.cell_type,
                   inputs.no_model_mask,
                   0.0,        // do not limit maximum thickness
                   m_Wtillnew,
                   m_grounded_margin_change,
                   m_grounding_line_change,
                   m_conservation_error_change,
                   m_no_model_mask_change);

    // update Wnew from W, Wtill, Wtillnew, Wstag, Q, input_rate
    update_W(hdt,
             m_input_rate,
             m_W, m_Wstag,
             m_Wtill, m_Wtillnew,
             m_K, m_Q,
             m_Wnew);
    // remove water in ice-free areas and account for changes
    enforce_bounds(*inputs.cell_type,
                   inputs.no_model_mask,
                   0.0,        // do not limit maximum thickness
                   m_Wnew,
                   m_grounded_margin_change,
                   m_grounding_line_change,
                   m_conservation_error_change,
                   m_no_model_mask_change);

    // transfer new into old
    m_W.copy_from(m_Wnew);
    m_Wtill.copy_from(m_Wtillnew);
  } // end of the time-stepping loop

  m_log->message(2,
                 "  took %d hydrology sub-steps with average dt = %.6f years (%.3f s or %.3f hours)\n",
                 step_counter,
                 units::convert(m_sys, dt / step_counter, "seconds", "years"),
                 dt / step_counter,
                 (dt / step_counter) / 3600.0);
}

const IceModelVec2Stag& Routing::velocity_staggered() const {
  return m_V;
}

std::map<std::string, Diagnostic::Ptr> Routing::diagnostics_impl() const {
  using namespace diagnostics;

  DiagnosticList result = {
    {"bwatvel",  Diagnostic::Ptr(new BasalWaterVelocity(this))},
    {"bwp",      Diagnostic::Ptr(new BasalWaterPressure(this))},
    {"bwprel",   Diagnostic::Ptr(new RelativeBasalWaterPressure(this))},
    {"effbwp",   Diagnostic::Ptr(new EffectiveBasalWaterPressure(this))},
    {"wallmelt", Diagnostic::Ptr(new WallMelt(this))},
  };
  return combine(result, Hydrology::diagnostics_impl());
}

std::map<std::string, TSDiagnostic::Ptr> Routing::ts_diagnostics_impl() const {
  std::map<std::string, TSDiagnostic::Ptr> result = {
    // FIXME: add mass-conservation diagnostics
  };
  return result;
}

} // end of namespace hydrology
} // end of namespace pism
