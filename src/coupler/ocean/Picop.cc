// Copyright (C) 2012-2019, 2021, 2022, 2023, 2024, 2025 Constantine Khrulev, Ricarda Winkelmann, Ronja Reese, Torsten
// Albrecht, and Matthias Mengel
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
//
// Please cite this model as:
// 1.
// Antarctic sub-shelf melt rates via PICO
// R. Reese, T. Albrecht, M. Mengel, X. Asay-Davis and R. Winkelmann
// The Cryosphere, 12, 1969-1985, (2018)
// DOI: 10.5194/tc-12-1969-2018
//
// 2.
// A box model of circulation and melting in ice shelf caverns
// D. Olbers & H. Hellmer
// Ocean Dynamics (2010), Volume 60, Issue 1, pp 141–153
// DOI: 10.1007/s10236-009-0252-z
//
// 3.
// PICOP, a new ocean melt parameterization under ice shelves
// combining PICO and a plume model.
// T. Pelle, M. Morlighem, J.H. Bondzio
// The Cryosphere, 13, 1043-49, (2019)
// DOI: 10.5194/tc-13-1043-2019

#include <cmath>
#include <stdexcept>

#include "pism/coupler/util/options.hh"
#include "pism/geometry/Geometry.hh"
#include "pism/util/Config.hh"
#include "pism/util/Grid.hh"
#include "pism/stressbalance/StressBalance.hh"

#include "pism/coupler/ocean/Picop.hh"
#include "pism/coupler/ocean/PicopPhysics.hh"
#include "pism/util/Logger.hh"
#include "pism/util/pism_utilities.hh"

namespace pism {

namespace ocean {

Picop::Picop(std::shared_ptr<const Grid> grid)
  : CompleteOceanModel(grid, std::shared_ptr<OceanModel>()),
    m_pico(std::make_shared<Pico>(grid)),
    m_basal_melt_rate(m_grid, "picop_basal_melt_rate"),
    m_grounding_line_elevation(grid, "picop_grounding_line_elevation"),
    m_grounding_line_slope(grid, "picop_grounding_line_slope"),
    m_theta_ocean(m_pico->get_temperature()),
    m_salinity_ocean(m_pico->get_salinity()),
    m_flow_direction(grid, "ice_flow_direction"),
    m_work(grid, "temporary_storage")
{

  ForcingOptions opt(*m_grid->ctx(), "ocean.picop");

  m_grounding_line_elevation.metadata(0).long_name("grounding line elevation").units("m");
  m_grounding_line_elevation.metadata()["_FillValue"] = { 0.0 };
  m_grounding_line_elevation.set(0.0);

  m_grounding_line_slope.metadata(0).long_name("grounding line slope").units("rad");
  m_grounding_line_slope.metadata()["_FillValue"] = { 0.0 };
  m_grounding_line_slope.set(0.0);
  
}

void Picop::init_impl(const Geometry &geometry) {
  (void) geometry;

  m_pico->init(geometry);
  m_log->message(2, "* Initializing the Plume extension of PICO for the ocean ...\n");

  PicopPhysics picop_physics(*m_config);

  double
    ice_density   = m_config->get_number("constants.ice.density"),
    water_density = m_config->get_number("constants.sea_water.density"),
    g             = m_config->get_number("constants.standard_gravity");

  compute_average_water_column_pressure(geometry, ice_density, water_density, g,
                                        *m_water_column_pressure);
}

void Picop::define_model_state_impl(const File &output) const {

  m_pico->define_model_state(output);
  OceanModel::define_model_state_impl(output);
}

void Picop::write_model_state_impl(const File &output) const {

  m_pico->write_model_state(output);
  OceanModel::write_model_state_impl(output);
}

void Picop::update_impl(const Inputs &inputs, double t, double dt) {

  m_pico->update(inputs, t, dt);

  const auto &geometry = *inputs.geometry;

  double
    ice_density   = m_config->get_number("constants.ice.density"),
    water_density = m_config->get_number("constants.sea_water.density"),
    g             = m_config->get_number("constants.standard_gravity");

  compute_average_water_column_pressure(geometry, ice_density, water_density, g,
                                        *m_water_column_pressure);

  if (inputs.stress_balance == nullptr) {
    // Use outputs from PICO if the stress balance is not available
    m_shelf_base_temperature->copy_from(m_pico->shelf_base_temperature());
    m_shelf_base_mass_flux->copy_from(m_pico->shelf_base_mass_flux());
    return;
  }
  
  PicopPhysics picop_physics(*m_config);
  
  compute_grounding_line_elevation(inputs, m_grounding_line_elevation);

  {
    // FIXME: remove this once the rest of PICOP is ready to test
    m_shelf_base_temperature->copy_from(m_pico->shelf_base_temperature());
    m_shelf_base_mass_flux->copy_from(m_pico->shelf_base_mass_flux());
    return;
  }

  compute_grounding_line_slope(inputs, m_grounding_line_slope);
  
  compute_melt_rate(inputs, picop_physics, m_theta_ocean, m_salinity_ocean,
                    m_grounding_line_elevation,m_grounding_line_slope,
                    m_basal_melt_rate);

  m_shelf_base_mass_flux->copy_from(m_basal_melt_rate);
  m_shelf_base_mass_flux->scale(ice_density);

}


MaxTimestep Picop::max_timestep_impl(double t) const {
  (void) t;

  auto pico_dt_max = m_pico->max_timestep(t);
  if (pico_dt_max.finite()) {
    return { pico_dt_max.value(), "ocean picop" };
  }

  return MaxTimestep("ocean picop");
}

void Picop::compute_melt_rate(const Inputs &inputs,
                              const PicopPhysics &physics,
                              const array::Scalar &T_a,
                              const array::Scalar &S_a,
                              array::Scalar1 &grounding_line_elevation,
                              array::Scalar1 &grounding_line_slope,
                              array::Scalar1 &melt_rate) const {

  const auto &ice_surface_elevation = inputs.geometry->ice_surface_elevation;
  const auto &ice_thickness = inputs.geometry->ice_thickness;
  
  array::AccessScope scope{&T_a, &S_a,
                           &ice_surface_elevation, &ice_thickness,
                           &grounding_line_slope, &grounding_line_elevation,
                           &melt_rate};

  for (auto p = m_grid->points(); p; p.next()) {
    int i = p.i(), j = p.j();

    const double z_b = ice_surface_elevation(i, j) - ice_thickness(i, j);
    const double z_gl = grounding_line_elevation(i, j);
    const double alpha = 0.0;  // FIX
    const double s_a = S_a(i, j);
    const double t_a = T_a(i, j);
      
    const double t_f_gl = physics.characteristic_freezing_poing(s_a, z_b);
    const double Gamma_TS = physics.effective_heat_exchange_coefficient(t_a, t_f_gl, alpha);
    const double l = physics.length_scaling(t_a, t_f_gl, alpha);
    const double g_alpha = physics.geometric_scaling(Gamma_TS, alpha);
    const double M = physics.melt_function(t_a, s_a, z_gl, g_alpha);
    const double X_hat = physics.dimensionless_coordinate(z_b, z_gl, l);
    
    melt_rate(i, j) = M * physics.dimensionless_melt_curve(X_hat);
  }
  
}


/*! Use bilinear interpolation to estimate the value in a cell with a,b,c,d at corners.
 *
 * `alpha` and `beta` are interpolation weights.
 *
 * d--------c
 * |        |
 * |        |
 * |        |
 * a--------b
 */
static double interpolate(double a, double b, double c, double d,
                          double alpha, double beta) {
  return a * (1 - alpha) * (1 - beta) + b * alpha * (1 - beta) + c * alpha * beta +
         d * (1 - alpha) * (beta);
}

/*! Use bilinear interpolation to estimate the value in a "box" stencil at the location (x,y).
 *
 * Each side of the box has length 2, from -1 to 1 in both X and Y directions. The point
 * "c" is at (0,0).
 *
 * nw-----n-----ne
 *  |     |     |
 *  |     |     |
 *  w-----c-----e
 *  |     |     |
 *  |     |     |
 * sw-----s-----se
 */
static double interpolate(const stencils::Box<double> &B, double x, double y) {
  if (x <= 0 and y <= 0) {
    /*!
     *  w--------c
     *  |        |
     *  |        |
     *  |        |
     * sw--------s
     */
    return interpolate(B.sw, B.s, B.c, B.w, 1 - (-x), 1 - (-y));
  }

  if (x > 0 and y <= 0) {
    /*!
     *  c--------e
     *  |        |
     *  |        |
     *  |        |
     *  s--------se
     */
    return interpolate(B.s, B.se, B.e, B.c, x, 1 - (-y));
  }

  if (x <= 0 and y > 0) {
    /*!
     * nw--------n
     *  |        |
     *  |        |
     *  |        |
     *  w--------c
     */
    return interpolate(B.w, B.c, B.n, B.nw, 1 - (-x), y);
  }

  if (x > 0 and y > 0) {
    /*!
     *  n--------ne
     *  |        |
     *  |        |
     *  |        |
     *  c--------e
     */
    return interpolate(B.c, B.e, B.ne, B.n, x, y);
  }

  throw std::runtime_error("logic failed: one of inputs might be nan");
}

/*!
 * Perform one iteration of the naive semi-Lagrangian transport algorithm, transporting
 * `U_old` over the area covered by floating ice in the direction defined by
 * `flow_direction` (normalized flow velocity), storing result in `U_new`.
 *
 * Values at locations where `cell_type` is other than floating_ice are held constant.
 */
void transport_step(const array::Scalar1 &U_old, const array::CellType &cell_type,
               const array::Vector &flow_direction, array::Scalar &U_new) {

  auto grid = U_new.grid();

  array::AccessScope scope{ &U_old, &cell_type, &U_new, &flow_direction };

  for (Points p(*grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (cell_type.floating_ice(i, j)) {
      auto U = U_old.box(i, j);
      auto D = flow_direction(i, j);

      U_new(i, j) = interpolate(U, -D.u, -D.v);
    } else {
      U_new(i, j) = U_old(i, j);
    }
  }
}

void Picop::compute_grounding_line_elevation(const Inputs &inputs,
                                             array::Scalar1 &result) {

  const auto &cell_type = inputs.geometry->cell_type;
  const auto &adv_vel   = inputs.stress_balance->advective_velocity();

  array::AccessScope scope{ &adv_vel, &m_flow_direction };

  // Step 1: Initialize zgl0 at grounding line: bed elevation
  //         Normalize velocities
  for (auto p = m_grid->points(); p; p.next()) {
    int i = p.i(), j = p.j();
    double flow_speed = adv_vel(i, j).magnitude();
    if (flow_speed > 0.0) {
      m_flow_direction(i, j) = adv_vel(i, j) / flow_speed;
    } else {
      m_flow_direction(i, j) = 0.0;
    }
  }

  // FIXME: this is the right way to initialize it *if* we don't have a better guess, but
  // we may benefit from re-using the result of this computation from the previous time
  // step, especially if the bed elevation did not change
  result.copy_from(inputs.geometry->bed_elevation);
  
  double tol = 1.0;             // meters
  double max_iter = 500;
  
  const auto &ice_surface_elevation = inputs.geometry->ice_surface_elevation;
  const auto &ice_thickness = inputs.geometry->ice_thickness;

  auto &result_old = result;
  auto &result_new = m_work;

  scope.add({ &result_old, &result_new, &cell_type, &ice_surface_elevation, &ice_thickness });

  for (int iter = 0; iter < max_iter; ++iter) {

    transport_step(result_old, cell_type, m_flow_direction, result_new);

    // elevation of a plume origin that reached (x,y) cannot be above the elevation of the
    // bottom of the ice at (x,y)
    for (auto p = m_grid->points(); p; p.next()) {
      int i = p.i(), j = p.j();

      if (cell_type.floating_ice(i, j)) {
        double ice_bottom_elevation = ice_surface_elevation(i, j) - ice_thickness(i, j);

        if (result_new(i, j) > ice_bottom_elevation) {
          result_new(i, j) = ice_bottom_elevation;
        }
      }
    }

    double residual = 0.0;
    for (auto p = m_grid->points(); p; p.next()) {
      int i = p.i(), j = p.j();
      residual = std::max(residual, std::abs(result_new(i, j) - result_old(i, j)));
    }

    residual = GlobalMax(m_grid->com, residual);

    // copy into `result`, updating ghosts for the next iteration
    result.copy_from(result_new);

    m_log->message(2, "picop iteration %03d, max change = %f m\n", iter, residual);

    if (residual < tol) {
      break;
    }
  }
}

void Picop::compute_grounding_line_slope(const Inputs &inputs,
                                             array::Scalar1 &result) {

  const auto &cell_type = inputs.geometry->cell_type;
  const auto &adv_vel   = inputs.stress_balance->advective_velocity();

  array::AccessScope scope{ &adv_vel, &m_flow_direction };

  // Step 1: Initialize zgl0 at grounding line: bed elevation
  //         Normalize velocities
  for (auto p = m_grid->points(); p; p.next()) {
    int i = p.i(), j = p.j();
    double flow_speed = adv_vel(i, j).magnitude();
    if (flow_speed > 0.0) {
      m_flow_direction(i, j) = adv_vel(i, j) / flow_speed;
    } else {
      m_flow_direction(i, j) = 0.0;
    }
  }

  // FIXME: this is the right way to initialize it *if* we don't have a better guess, but
  // we may benefit from re-using the result of this computation from the previous time
  // step, especially if the bed elevation did not change
  result.copy_from(inputs.geometry->bed_elevation);
  
  double tol = 1.0;             // meters
  double max_iter = 500;
  
  const auto &ice_surface_elevation = inputs.geometry->ice_surface_elevation;
  const auto &ice_thickness = inputs.geometry->ice_thickness;

  auto &result_old = result;
  auto &result_new = m_work;

  scope.add({ &result_old, &result_new, &cell_type, &ice_surface_elevation, &ice_thickness });

  for (int iter = 0; iter < max_iter; ++iter) {

    transport_step(result_old, cell_type, m_flow_direction, result_new);

    // elevation of a plume origin that reached (x,y) cannot be above the elevation of the
    // bottom of the ice at (x,y)
    for (auto p = m_grid->points(); p; p.next()) {
      int i = p.i(), j = p.j();

      if (cell_type.floating_ice(i, j)) {
        double ice_bottom_elevation = ice_surface_elevation(i, j) - ice_thickness(i, j);

        if (result_new(i, j) > ice_bottom_elevation) {
          result_new(i, j) = ice_bottom_elevation;
        }
      }
    }

    double residual = 0.0;
    for (auto p = m_grid->points(); p; p.next()) {
      int i = p.i(), j = p.j();
      residual = std::max(residual, std::abs(result_new(i, j) - result_old(i, j)));
    }

    residual = GlobalMax(m_grid->com, residual);

    // copy into `result`, updating ghosts for the next iteration
    result.copy_from(result_new);

    m_log->message(2, "picop iteration %03d, max change = %f m\n", iter, residual);

    if (residual < tol) {
      break;
    }
  }
}

// Write diagnostic variables to extra files if requested
DiagnosticList Picop::diagnostics_impl() const {

  DiagnosticList result = {
    { "picop_grounding_line_elevation", Diagnostic::wrap(m_grounding_line_elevation) },
    { "picop_grounding_line_slope", Diagnostic::wrap(m_grounding_line_slope) },
  };

  return combine(result, OceanModel::diagnostics_impl());
}


} // end of namespace ocean
} // end of namespace pism
