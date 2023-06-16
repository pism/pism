// Copyright (C) 2012-2023 PISM Authors
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
#include "pism/util/array/CellType.hh"
#include "pism/util/Mask.hh"
#include "pism/util/MaxTimestep.hh"

#include "pism/util/error_handling.hh"

#include "pism/util/pism_options.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/util/Vars.hh"
#include "pism/geometry/Geometry.hh"
#include "pism/util/Profiling.hh"
#include "pism/util/Context.hh"

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
  virtual array::Array::Ptr compute_impl() const {
    array::Scalar::Ptr result(new array::Scalar(m_grid, "bwp"));
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
    m_vars[0]["_FillValue"] = {m_fill_value};
  }

protected:
  virtual array::Array::Ptr compute_impl() const {
    double fill_value = m_fill_value;

    array::Scalar::Ptr result(new array::Scalar(m_grid, "bwprel"));
    result->metadata(0) = m_vars[0];

    const array::Scalar
      &P  = model->subglacial_water_pressure(),
      &Po = model->overburden_pressure();

    array::AccessScope list{result.get(), &Po, &P};
    for (auto p = m_grid->points(); p; p.next()) {
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
  virtual array::Array::Ptr compute_impl() const {

    array::Scalar::Ptr result(new array::Scalar(m_grid, "effbwp"));
    result->metadata() = m_vars[0];

    const array::Scalar
      &P  = model->subglacial_water_pressure(),
      &Po = model->overburden_pressure();

    array::AccessScope list{&Po, &P, result.get()};

    for (auto p = m_grid->points(); p; p.next()) {
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
  virtual array::Array::Ptr compute_impl() const {
    array::Scalar::Ptr result(new array::Scalar(m_grid, "wallmelt"));
    result->metadata() = m_vars[0];

    const array::Scalar &bed_elevation = *m_grid->variables().get_2d_scalar("bedrock_altitude");

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
              "m s-1", "m s-1", 0);
    set_attrs("velocity of water in subglacial layer, j-offset", "",
              "m s-1", "m s-1", 1);
  }
protected:
  virtual array::Array::Ptr compute_impl() const {
    auto result = std::make_shared<array::Staggered>(m_grid, "bwatvel");
    result->metadata(0) = m_vars[0];
    result->metadata(1) = m_vars[1];

    result->copy_from(model->velocity_staggered());

    return result;
  }
};

//! Compute the hydraulic potential.
/*!
  Computes \f$\psi = P + \rho_w g (b + W)\f$.
*/
void hydraulic_potential(const array::Scalar &W,
                         const array::Scalar &P,
                         const array::Scalar &sea_level,
                         const array::Scalar &bed,
                         const array::Scalar &ice_thickness,
                         array::Scalar &result) {

  auto grid = result.grid();

  Config::ConstPtr config = grid->ctx()->config();

  double
    ice_density       = config->get_number("constants.ice.density"),
    sea_water_density = config->get_number("constants.sea_water.density"),
    C                 = ice_density / sea_water_density,
    rg                = (config->get_number("constants.fresh_water.density") *
                         config->get_number("constants.standard_gravity"));

  array::AccessScope list{&P, &W, &sea_level, &ice_thickness, &bed, &result};

  for (auto p = grid->points(); p; p.next()) {
    const int i = p.i(), j = p.j();

    double b = std::max(bed(i, j), sea_level(i, j) - C * ice_thickness(i, j));

    result(i, j) = P(i, j) + rg * (b + W(i, j));
  }
}

/*! @brief Report hydraulic potential in the subglacial hydrology system */
class HydraulicPotential : public Diag<Routing>
{
public:
  HydraulicPotential(const Routing *m)
    : Diag<Routing>(m) {

    m_vars = {SpatialVariableMetadata(m_sys, "hydraulic_potential")};

    set_attrs("hydraulic potential in the subglacial hydrology system", "",
              "Pa", "Pa", 0);
  }

protected:
  array::Array::Ptr compute_impl() const {

    array::Scalar::Ptr result(new array::Scalar(m_grid, "hydraulic_potential"));
    result->metadata(0) = m_vars[0];

    const array::Scalar        &sea_level     = *m_grid->variables().get_2d_scalar("sea_level");
    const array::Scalar        &bed_elevation = *m_grid->variables().get_2d_scalar("bedrock_altitude");
    const array::Scalar        &ice_thickness = *m_grid->variables().get_2d_scalar("land_ice_thickness");

    hydraulic_potential(model->subglacial_water_thickness(),
                        model->subglacial_water_pressure(),
                        sea_level,
                        bed_elevation,
                        ice_thickness,
                        *result);

    return result;
  }
};

} // end of namespace diagnostics

Routing::Routing(std::shared_ptr<const IceGrid> grid)
  : Hydrology(grid),
    m_Qstag(grid, "advection_flux"),
    m_Qstag_average(grid, "cumulative_advection_flux"),
    m_Vstag(grid, "water_velocity"),
    m_Wstag(grid, "W_staggered"),
    m_Kstag(grid, "K_staggered"),
    m_Wnew(grid, "W_new"),
    m_Wtillnew(grid, "Wtill_new"),
    m_R(grid, "potential_workspace"), /* box stencil used */
    m_dx(grid->dx()),
    m_dy(grid->dy()),
    m_bottom_surface(grid, "ice_bottom_surface_elevation") {

  m_W.metadata()["pism_intent"] = "model_state";

  m_rg = (m_config->get_number("constants.fresh_water.density") *
          m_config->get_number("constants.standard_gravity"));

  m_Qstag.set_attrs("internal",
                    "cell face-centered (staggered) components of advective subglacial water flux",
                    "m2 s-1", "m2 s-1", "", 0);

  m_Qstag_average.set_attrs("internal",
                            "average (over time) advection flux on the staggered grid",
                            "m2 s-1", "m2 s-1", "", 0);

  m_Vstag.set_attrs("internal",
                    "cell face-centered (staggered) components of water velocity"
                    " in subglacial water layer",
                    "m s-1", "m s-1", "", 0);

  // auxiliary variables which NEED ghosts
  m_Wstag.set_attrs("internal",
                    "cell face-centered (staggered) values of water layer thickness",
                    "m", "m", "", 0);
  m_Wstag.metadata()["valid_min"] = {0.0};

  m_Kstag.set_attrs("internal",
                    "cell face-centered (staggered) values of nonlinear conductivity",
                    "", "", "", 0);
  m_Kstag.metadata()["valid_min"] = {0.0};

  m_R.set_attrs("internal",
                "work space for modeled subglacial water hydraulic potential",
                "Pa", "Pa", "", 0);

  // temporaries during update; do not need ghosts
  m_Wnew.set_attrs("internal",
                   "new thickness of transportable subglacial water layer during update",
                   "m", "m", "", 0);
  m_Wnew.metadata()["valid_min"] = {0.0};

  m_Wtillnew.set_attrs("internal",
                       "new thickness of till (subglacial) water layer during update",
                       "m", "m", "", 0);
  m_Wtillnew.metadata()["valid_min"] = {0.0};

  {
    double alpha = m_config->get_number("hydrology.thickness_power_in_flux");
    if (alpha < 1.0) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                    "alpha = %f < 1 which is not allowed", alpha);
    }

    if (m_config->get_number("hydrology.tillwat_max") < 0.0) {
      throw RuntimeError(PISM_ERROR_LOCATION,
                         "hydrology::Routing: hydrology.tillwat_max is negative.\n"
                         "This is not allowed.");
    }
  }
}

void Routing::initialization_message() const {
  m_log->message(2,
                 "* Initializing the routing subglacial hydrology model ...\n");

  if (m_config->get_flag("hydrology.routing.include_floating_ice")) {
    m_log->message(2, "  ... routing subglacial water under grounded and floating ice.\n");
  } else {
    m_log->message(2, "  ... routing subglacial water under grounded ice only.\n");
  }
}

void Routing::restart_impl(const File &input_file, int record) {
  Hydrology::restart_impl(input_file, record);

  m_W.read(input_file, record);

  regrid("Hydrology", m_W);
}

void Routing::bootstrap_impl(const File &input_file,
                             const array::Scalar &ice_thickness) {
  Hydrology::bootstrap_impl(input_file, ice_thickness);

  double bwat_default = m_config->get_number("bootstrapping.defaults.bwat");
  m_W.regrid(input_file, OPTIONAL, bwat_default);

  regrid("Hydrology", m_W);
}

void Routing::init_impl(const array::Scalar &W_till,
                              const array::Scalar &W,
                              const array::Scalar &P) {
  Hydrology::init_impl(W_till, W, P);

  m_W.copy_from(W);
}

void Routing::define_model_state_impl(const File &output) const {
  Hydrology::define_model_state_impl(output);
  m_W.define(output);
}

void Routing::write_model_state_impl(const File &output) const {
  Hydrology::write_model_state_impl(output);
  m_W.write(output);
}

//! Returns the (trivial) overburden pressure as the pressure of the transportable water,
//! because this is the model.
const array::Scalar& Routing::subglacial_water_pressure() const {
  return m_Pover;
}

const array::Staggered& Routing::velocity_staggered() const {
  return m_Vstag;
}


//! Average the regular grid water thickness to values at the center of cell edges.
/*! Uses mask values to avoid averaging using water thickness values from
  either ice-free or floating areas. */
void Routing::water_thickness_staggered(const array::Scalar &W,
                                        const array::CellType1 &mask,
                                        array::Staggered &result) {

  bool include_floating = m_config->get_flag("hydrology.routing.include_floating_ice");

  array::AccessScope list{ &mask, &W, &result };

  for (auto p = m_grid->points(); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (include_floating) {
      // east
      if (mask.icy(i, j)) {
        result(i, j, 0) = mask.icy(i + 1, j) ? 0.5 * (W(i, j) + W(i + 1, j)) : W(i, j);
      } else {
        result(i, j, 0) = mask.icy(i + 1, j) ? W(i + 1, j) : 0.0;
      }
      // north
      if (mask.icy(i, j)) {
        result(i, j, 1) = mask.icy(i, j + 1) ? 0.5 * (W(i, j) + W(i, j + 1)) : W(i, j);
      } else {
        result(i, j, 1) = mask.icy(i, j + 1) ? W(i, j + 1) : 0.0;
      }
    } else {
      // east
      if (mask.grounded_ice(i, j)) {
        result(i, j, 0) = mask.grounded_ice(i + 1, j) ? 0.5 * (W(i, j) + W(i + 1, j)) : W(i, j);
      } else {
        result(i, j, 0) = mask.grounded_ice(i + 1, j) ? W(i + 1, j) : 0.0;
      }
      // north
      if (mask.grounded_ice(i, j)) {
        result(i, j, 1) = mask.grounded_ice(i, j + 1) ? 0.5 * (W(i, j) + W(i, j + 1)) : W(i, j);
      } else {
        result(i, j, 1) = mask.grounded_ice(i, j + 1) ? W(i, j + 1) : 0.0;
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
void Routing::compute_conductivity(const array::Staggered &W,
                                   const array::Scalar &P,
                                   const array::Scalar &bed_elevation,
                                   array::Staggered &result,
                                   double &KW_max) const {
  const double
    k     = m_config->get_number("hydrology.hydraulic_conductivity"),
    alpha = m_config->get_number("hydrology.thickness_power_in_flux"),
    beta  = m_config->get_number("hydrology.gradient_power_in_flux"),
    betapow = (beta - 2.0) / 2.0;

  array::AccessScope list({&result, &W});

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
      for (auto p = m_grid->points(); p; p.next()) {
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

    for (auto p = m_grid->points(); p; p.next()) {
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
    for (auto p = m_grid->points(); p; p.next()) {
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
               const array::Scalar &bed_elevation,
               array::Scalar &result) {

  auto grid = result.grid();

  Config::ConstPtr config = grid->ctx()->config();

  const double
    k     = config->get_number("hydrology.hydraulic_conductivity"),
    L     = config->get_number("constants.fresh_water.latent_heat_of_fusion"),
    alpha = config->get_number("hydrology.thickness_power_in_flux"),
    beta  = config->get_number("hydrology.gradient_power_in_flux"),
    g     = config->get_number("constants.standard_gravity"),
    rhow  = config->get_number("constants.fresh_water.density"),
    rg    = rhow * g,
    CC    = k / (L * rhow);

  // FIXME:  could be scaled with overall factor hydrology_coefficient_wall_melt ?
  if (alpha < 1.0) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "alpha = %f < 1 which is not allowed", alpha);
  }

  array::Scalar1 R(grid, "R");

  // R  <-- P + rhow g b
  model.subglacial_water_pressure().add(rg, bed_elevation, R);
  // yes, it updates ghosts

  array::Scalar1 W(grid, "W");
  W.copy_from(model.subglacial_water_thickness());

  array::AccessScope list{&R, &W, &result};

  double dx = grid->dx();
  double dy = grid->dy();

  for (auto p = grid->points(); p; p.next()) {
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
void Routing::compute_velocity(const array::Staggered &W,
                               const array::Scalar &pressure,
                               const array::Scalar &bed,
                               const array::Staggered &K,
                               const array::Scalar1 *no_model_mask,
                               array::Staggered &result) const {
  array::Scalar &P = m_R;
  P.copy_from(pressure);  // yes, it updates ghosts

  array::AccessScope list{&P, &W, &K, &bed, &result};

  for (auto p = m_grid->points(); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (W(i, j, 0) > 0.0) {
      double
        P_x = (P(i + 1, j) - P(i, j)) / m_dx,
        b_x = (bed(i + 1, j) - bed(i, j)) / m_dx;
      result(i, j, 0) = - K(i, j, 0) * (P_x + m_rg * b_x);
    } else {
      result(i, j, 0) = 0.0;
    }

    if (W(i, j, 1) > 0.0) {
      double
        P_y = (P(i, j + 1) - P(i, j)) / m_dy,
        b_y = (bed(i, j + 1) - bed(i, j)) / m_dy;
      result(i, j, 1) = - K(i, j, 1) * (P_y + m_rg * b_y);
    } else {
      result(i, j, 1) = 0.0;
    }
  }

  if (no_model_mask) {
    list.add(*no_model_mask);

    for (auto p = m_grid->points(); p; p.next()) {
      const int i = p.i(), j = p.j();

      auto M = no_model_mask->star(i, j);

      if (M.c or M.e) {
        result(i, j, 0) = 0.0;
      }

      if (M.c or M.n) {
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
void Routing::advective_fluxes(const array::Staggered &V,
                               const array::Scalar &W,
                               array::Staggered &result) const {
  array::AccessScope list{&W, &V, &result};

  assert(W.stencil_width() >= 1);

  for (auto p = m_grid->points(); p; p.next()) {
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
  double result = 1.0 / (m_dx * m_dx) + 1.0 / (m_dy * m_dy);
  return 0.25 / (D_max * result);
}

/*!
 * See equation (50) in Bueler and van Pelt.
 */
double Routing::max_timestep_W_cfl() const {
  // V could be zero if P is constant and bed is flat
  auto tmp = absmax(m_Vstag);

  // add a safety margin
  double alpha = 0.95;
  double eps = 1e-6;

  return alpha * 0.5 / (tmp[0]/m_dx + tmp[1]/m_dy + eps);
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
                           const array::Scalar &Wtill,
                           const array::Scalar &surface_input_rate,
                           const array::Scalar &basal_melt_rate,
                           array::Scalar &Wtill_new) {
  const double
    tillwat_max = m_config->get_number("hydrology.tillwat_max"),
    C           = m_config->get_number("hydrology.tillwat_decay_rate", "m / second");

  array::AccessScope list{&Wtill, &Wtill_new, &basal_melt_rate};

  bool add_surface_input = m_config->get_flag("hydrology.add_water_input_to_till_storage");
  if (add_surface_input) {
    list.add(surface_input_rate);
  }

  for (auto p = m_grid->points(); p; p.next()) {
    const int i = p.i(), j = p.j();

    double input_rate = basal_melt_rate(i, j);
    if (add_surface_input) {
      input_rate += surface_input_rate(i, j);
    }

    Wtill_new(i, j) = clip(Wtill(i, j) + dt * (input_rate - C),
                           0.0, tillwat_max);
  }
}

void Routing::W_change_due_to_flow(double dt,
                                   const array::Scalar1    &W,
                                   const array::Staggered1 &Wstag,
                                   const array::Staggered1 &K,
                                   const array::Staggered1 &Q,
                                   array::Scalar &result) {
  const double
    wux = 1.0 / (m_dx * m_dx),
    wuy = 1.0 / (m_dy * m_dy);

  array::AccessScope list{&W, &Wstag, &K, &Q, &result};

  for (auto p = m_grid->points(); p; p.next()) {
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
    const double diffW = (wux * (De * (w.e - w.c) - Dw * (w.c - w.w)) +
                          wuy * (Dn * (w.n - w.c) - Ds * (w.c - w.s)));

    result(i, j) = dt * (- divQ + diffW);
  }
}


//! The computation of Wnew, called by update().
void Routing::update_W(double dt,
                       const array::Scalar     &surface_input_rate,
                       const array::Scalar     &basal_melt_rate,
                       const array::Scalar1    &W,
                       const array::Staggered1 &Wstag,
                       const array::Scalar     &Wtill,
                       const array::Scalar     &Wtill_new,
                       const array::Staggered1 &K,
                       const array::Staggered1 &Q,
                       array::Scalar &W_new) {

  W_change_due_to_flow(dt, W, Wstag, K, Q, m_flow_change_incremental);

  array::AccessScope list{&W, &Wtill, &Wtill_new, &surface_input_rate,
                               &basal_melt_rate, &m_flow_change_incremental, &W_new};

  for (auto p = m_grid->points(); p; p.next()) {
    const int i = p.i(), j = p.j();

    double input_rate = surface_input_rate(i, j) + basal_melt_rate(i, j);

    double Wtill_change = Wtill_new(i, j) - Wtill(i, j);
    W_new(i, j) = (W(i, j) + (dt * input_rate - Wtill_change) + m_flow_change_incremental(i, j));
  }

  m_flow_change.add(1.0, m_flow_change_incremental);
  m_input_change.add(dt, surface_input_rate);
  m_input_change.add(dt, basal_melt_rate);
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

  ice_bottom_surface(*inputs.geometry, m_bottom_surface);

  double
    ht  = t,
    hdt = 0.0;

  const double
    t_final     = t + dt,
    dt_max      = m_config->get_number("hydrology.maximum_time_step", "seconds"),
    tillwat_max = m_config->get_number("hydrology.tillwat_max");

  m_Qstag_average.set(0.0);

  // make sure W has valid ghosts before starting hydrology steps
  m_W.update_ghosts();

  unsigned int step_counter = 0;
  for (; ht < t_final; ht += hdt) {
    step_counter++;

#if (Pism_DEBUG==1)
    double huge_number = 1e6;
    check_bounds(m_W, huge_number);
    check_bounds(m_Wtill, tillwat_max);
#endif

    // updates ghosts of m_Wstag
    water_thickness_staggered(m_W,
                              inputs.geometry->cell_type,
                              m_Wstag);

    double maxKW = 0.0;
    // updates ghosts of m_Kstag
    profiling().begin("routing_conductivity");
    compute_conductivity(m_Wstag,
                         subglacial_water_pressure(),
                         m_bottom_surface,
                         m_Kstag, maxKW);
    profiling().end("routing_conductivity");

    // ghosts of m_Vstag are not updated
    profiling().begin("routing_velocity");
    compute_velocity(m_Wstag,
                     subglacial_water_pressure(),
                     m_bottom_surface,
                     m_Kstag,
                     inputs.no_model_mask,
                     m_Vstag);
    profiling().end("routing_velocity");

    // to get Q, W needs valid ghosts (ghosts of m_Vstag are not used)
    // updates ghosts of m_Qstag
    profiling().begin("routing_flux");
    advective_fluxes(m_Vstag, m_W, m_Qstag);
    profiling().end("routing_flux");

    m_Qstag_average.add(hdt, m_Qstag);

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
    {
      profiling().begin("routing_Wtill");
      update_Wtill(hdt,
                   m_Wtill,
                   m_surface_input_rate,
                   m_basal_melt_rate,
                   m_Wtillnew);
      // remove water in ice-free areas and account for changes
      enforce_bounds(inputs.geometry->cell_type,
                     inputs.no_model_mask,
                     0.0,           // do not limit maximum thickness
                     tillwat_max,   // till water thickness under the ocean
                     m_Wtillnew,
                     m_grounded_margin_change,
                     m_grounding_line_change,
                     m_conservation_error_change,
                     m_no_model_mask_change);
      profiling().end("routing_Wtill");
    }

    // update Wnew from W, Wtill, Wtillnew, Wstag, Q, input_rate
    // uses ghosts of m_W, m_Wstag, m_Qstag, m_Kstag
    {
      profiling().begin("routing_W");
      update_W(hdt,
               m_surface_input_rate,
               m_basal_melt_rate,
               m_W, m_Wstag,
               m_Wtill, m_Wtillnew,
               m_Kstag, m_Qstag,
               m_Wnew);
      // remove water in ice-free areas and account for changes
      enforce_bounds(inputs.geometry->cell_type,
                     inputs.no_model_mask,
                     0.0,        // do not limit maximum thickness
                     0.0,        // transportable water thickness under the ocean
                     m_Wnew,
                     m_grounded_margin_change,
                     m_grounding_line_change,
                     m_conservation_error_change,
                     m_no_model_mask_change);

      // transfer new into old (updates ghosts of m_W)
      m_W.copy_from(m_Wnew);
      profiling().end("routing_W");
    }

    // m_Wtill has no ghosts
    m_Wtill.copy_from(m_Wtillnew);
  } // end of the time-stepping loop

  staggered_to_regular(inputs.geometry->cell_type, m_Qstag_average,
                       m_config->get_flag("hydrology.routing.include_floating_ice"),
                       m_Q);
  m_Q.scale(1.0 / dt);

  m_log->message(2,
                 "  took %d hydrology sub-steps with average dt = %.6f years (%.3f s or %.3f hours)\n",
                 step_counter,
                 units::convert(m_sys, dt / step_counter, "seconds", "years"),
                 dt / step_counter,
                 (dt / step_counter) / 3600.0);
}

std::map<std::string, Diagnostic::Ptr> Routing::diagnostics_impl() const {
  using namespace diagnostics;

  DiagnosticList result = {
    {"bwatvel",             Diagnostic::Ptr(new BasalWaterVelocity(this))},
    {"bwp",                 Diagnostic::Ptr(new BasalWaterPressure(this))},
    {"bwprel",              Diagnostic::Ptr(new RelativeBasalWaterPressure(this))},
    {"effbwp",              Diagnostic::Ptr(new EffectiveBasalWaterPressure(this))},
    {"wallmelt",            Diagnostic::Ptr(new WallMelt(this))},
    {"hydraulic_potential", Diagnostic::Ptr(new HydraulicPotential(this))},
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
