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
#include "pism/util/Profiling.hh"
#include "pism/util/pism_utilities.hh"

namespace pism {

namespace ocean {

Picop::Picop(std::shared_ptr<const Grid> grid)
  : CompleteOceanModel(grid, std::shared_ptr<OceanModel>()),
    m_pico(std::make_shared<Pico>(grid)),
    m_basal_melt_rate(m_grid, "picop_basal_melt_rate"),
    m_grounding_line_elevation(grid, "picop_grounding_line_elevation"),
    m_shelf_base_elevation(grid, "picop_shelf_base_elevation"),
    m_grounding_line_slope(grid, "picop_grounding_line_slope"),
    m_geometric_scale(grid, "picop_geometric_scale"),
    m_length_scale(grid, "picop_length_scale"),
    m_theta_ocean(m_pico->get_temperature()),
    m_salinity_ocean(m_pico->get_salinity()),
    m_flow_direction(grid, "ice_flow_direction"),
    m_work(grid, "temporary_storage"),
    m_zb_x(grid, "staggered slope x"),
    m_zb_y(grid, "staggered slope y") {

  ForcingOptions opt(*m_grid->ctx(), "ocean.picop");

  m_basal_melt_rate.metadata(0)
      .long_name("PICOP sub-shelf melt rate")
      .units("m s^-1")
      .output_units("m year^-1");
  m_basal_melt_rate.metadata()["_FillValue"] = {0.0};
  
  m_grounding_line_elevation.metadata(0).long_name("grounding line elevation").units("m");
  m_grounding_line_elevation.metadata()["_FillValue"] = { 0.0 };
  m_grounding_line_elevation.set(0.0);

  m_shelf_base_elevation.metadata(0).long_name("shelf base elevation").units("m");
  m_shelf_base_elevation.metadata()["_FillValue"] = { 0.0 };
  m_shelf_base_elevation.set(0.0);
  
  m_grounding_line_slope.metadata(0).long_name("shelf base slope").units("rad").output_units("degree");
  m_grounding_line_slope.metadata()["_FillValue"] = { 0.0 };
  m_grounding_line_slope.set(0.0);
      
  m_length_scale.metadata(0).long_name("length scale").units("m");
  m_length_scale.metadata()["_FillValue"] = { 0.0 };
  m_length_scale.set(0.0);

  m_geometric_scale.metadata(0).long_name("length scale").units("");
  m_geometric_scale.metadata()["_FillValue"] = { 0.0 };
  m_geometric_scale.set(0.0);

  m_shelf_base_temperature->metadata()["_FillValue"] = {0.0};
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

// CODE Duplication: should this be a public member of PICO?
/*!
* Extend basal melt rates to grounded and ocean neighbors for consitency with subgl_melt.
* Note that melt rates are then simply interpolated into partially floating cells, they
* are not included in the calculations of PICO.
*/
static void extend_basal_melt_rates(const array::CellType1 &cell_type,
                                    array::Scalar1 &basal_melt_rate) {

  auto grid = basal_melt_rate.grid();

  // update ghosts of the basal melt rate so that we can use basal_melt_rate.box(i,j)
  // below
  basal_melt_rate.update_ghosts();

  array::AccessScope list{&cell_type, &basal_melt_rate};

  for (auto p = grid->points(); p; p.next()) {

    const int i = p.i(), j = p.j();

    auto M = cell_type.box(i, j);

    bool potential_partially_filled_cell =
      ((M.c  == MASK_GROUNDED or M.c  == MASK_ICE_FREE_OCEAN) and
       (M.w  == MASK_FLOATING or M.e  == MASK_FLOATING or
        M.s  == MASK_FLOATING or M.n  == MASK_FLOATING or
        M.sw == MASK_FLOATING or M.nw == MASK_FLOATING or
        M.se == MASK_FLOATING or M.ne == MASK_FLOATING));

    if (potential_partially_filled_cell) {
      auto BMR = basal_melt_rate.box(i, j);

      int N = 0;
      double melt_sum = 0.0;

      melt_sum += M.nw == MASK_FLOATING ? (++N, BMR.nw) : 0.0;
      melt_sum += M.n  == MASK_FLOATING ? (++N, BMR.n)  : 0.0;
      melt_sum += M.ne == MASK_FLOATING ? (++N, BMR.ne) : 0.0;
      melt_sum += M.e  == MASK_FLOATING ? (++N, BMR.e)  : 0.0;
      melt_sum += M.se == MASK_FLOATING ? (++N, BMR.se) : 0.0;
      melt_sum += M.s  == MASK_FLOATING ? (++N, BMR.s)  : 0.0;
      melt_sum += M.sw == MASK_FLOATING ? (++N, BMR.sw) : 0.0;
      melt_sum += M.w  == MASK_FLOATING ? (++N, BMR.w)  : 0.0;

      if (N != 0) { // If there are floating neigbors, return average melt rates
        basal_melt_rate(i, j) = melt_sum / N;
      }
    }
  } // end of the loop over grid points
}

void Picop::update_impl(const Inputs &inputs, double t, double dt) {

  m_pico->update(inputs, t, dt);

  const auto &cell_type = inputs.geometry->cell_type;

  if (inputs.stress_balance == nullptr) {
    // Use outputs from PICO if the stress balance is not available
    m_shelf_base_temperature->copy_from(m_pico->shelf_base_temperature());
    m_shelf_base_mass_flux->copy_from(m_pico->shelf_base_mass_flux());
    return;
  }
  
  PicopPhysics picop_physics(*m_config);

  
  compute_shelf_base_elevation(inputs, m_shelf_base_elevation);

  profiling().begin("ocean.compute_grounding_line_elevation");
  compute_grounding_line_elevation(inputs, m_grounding_line_elevation);
  profiling().end("ocean.compute_grounding_line_elevation");
  

  profiling().begin("ocean.compute_grounding_line_slope");
  compute_grounding_line_slope(inputs, m_grounding_line_slope);
  profiling().end("ocean.compute_grounding_line_slope");
  
  compute_melt_rate(inputs, picop_physics, m_theta_ocean, m_salinity_ocean,
                    m_basal_melt_rate);

  extend_basal_melt_rates(cell_type, m_basal_melt_rate);
    
  double
    ice_density   = m_config->get_number("constants.ice.density"),
    water_density = m_config->get_number("constants.sea_water.density"),
    g             = m_config->get_number("constants.standard_gravity");
  
  m_shelf_base_temperature->copy_from(m_pico->shelf_base_temperature());
  m_shelf_base_mass_flux->copy_from(m_basal_melt_rate);
  m_shelf_base_mass_flux->scale(ice_density);

  compute_average_water_column_pressure(*inputs.geometry, ice_density, water_density, g,
                                        *m_water_column_pressure);
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
                              array::Scalar1 &result)  {

  const auto &ice_surface_elevation = inputs.geometry->ice_surface_elevation;
  const auto &ice_thickness = inputs.geometry->ice_thickness;
  const auto &cell_type = inputs.geometry->cell_type;
  
  array::AccessScope scope{&T_a, &S_a, &cell_type,
                           &ice_surface_elevation, &ice_thickness,
                           &m_geometric_scale, &m_length_scale,
                           &m_grounding_line_slope, &m_grounding_line_elevation,
                           &result};

  for (auto p = m_grid->points(); p; p.next()) {
    int i = p.i(), j = p.j();

    if (cell_type.floating_ice(i, j)) {
      const double z_b = ice_surface_elevation(i, j) - ice_thickness(i, j);
      const double z_gl = m_grounding_line_elevation(i, j);
      const double alpha = m_grounding_line_slope(i, j);
            
      const double s_a = S_a(i, j);
      const double t_min = physics.characteristic_freezing_point(s_a, 0.0);
      double t_a = T_a(i, j);
      /* Low bound for Toc to ensure X_hat is between 0 and 1 */
      if (t_a < t_min) {
        t_a = t_min;
      }
      
      const double t_f_gl = physics.characteristic_freezing_point(s_a, z_gl);
      const double Gamma_TS = physics.effective_heat_exchange_coefficient(t_a, t_f_gl, alpha);
      const double l = physics.length_scaling(t_a, t_f_gl, Gamma_TS, alpha);
      const double g_alpha = physics.geometric_scaling(Gamma_TS, alpha);
      m_geometric_scale(i, j) = g_alpha;
      m_length_scale(i, j) = l;
      double X_hat = physics.dimensionless_coordinate(z_b, z_gl, l);
      const double M = physics.melt_function(t_a, t_f_gl, g_alpha);
      double m =  M * physics.dimensionless_melt_curve(X_hat);

      m_log->message(2, "(%i, %i) s_a=%f, t_a=%f, t_f_gl = %f, t_min=%f, zb=%f, zgl=%f, Gamma_TS=%f, alpha=%f, l=%f, X_hat=%f\n", i, j, s_a, t_a, t_f_gl, t_min, z_b, z_gl, Gamma_TS, alpha, l, X_hat);
      result(i, j) = m;
    }    
  }
}

/*!
 * Compute the weight used to determine if the difference between locations `i,j` and `n`
 * (neighbor) should be used in the computation of the surface gradient in
 * SSA::compute_driving_stress().
 *
 * We avoid differencing across
 *
 * - ice margins if stress boundary condition at ice margins (CFBC) is active
 * - grounding lines
 * - ice margins next to ice free locations above the surface elevation of the ice (fjord
 *   walls, nunataks, headwalls)
 */
static int weight(bool margin_bc, int M_ij, int M_n, double h_ij, double h_n) {
  using mask::floating_ice;
  using mask::grounded;
  using mask::ice_free;
  using mask::ice_free_ocean;
  using mask::icy;

  // grounding lines and calving fronts
  if ((grounded(M_ij) and floating_ice(M_n)) or (floating_ice(M_ij) and grounded(M_n)) or
      (floating_ice(M_ij) and ice_free_ocean(M_n))) {
    return 0;
  }

  // fjord walls, nunataks, headwalls
  if ((icy(M_ij) and ice_free(M_n) and h_n > h_ij) or
      (ice_free(M_ij) and icy(M_n) and h_ij > h_n)) {
    return 0;
  }

  // This condition has to match the one used to implement the calving front stress
  // boundary condition in assemble_rhs().
  if (margin_bc and ((icy(M_ij) and ice_free(M_n)) or (ice_free(M_ij) and icy(M_n)))) {
    return 0;
  }

  return 1;
}

// Use "upwinded" (-ish) finite difference to approximate the surface elevation
// difference.
static double diff_uphill(double L, double C, double R) {
  double dL = C - L, dR = R - C;

  if (dR * dL > 0) {
    // dL and dR have the same sign
    //
    // If dL < 0 then L > C > R and "dL = C - L" is the "uphill" difference.
    //
    // If dL > 0 then L < C < R and "dR = R - C" is the "uphill" difference.
    return dL < 0.0 ? dL : dR;
  }

  // centered
  return 0.5 * (dL + dR);
}

static double diff_centered(double L, double /* unused */, double R) {
  return 0.5 * (R - L);
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

  array::AccessScope scope{ &adv_vel, &m_flow_direction, &result };

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
  
  const double rtol = 0.001;             // meters
  const int max_iter = 500;
  
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
      double denom = std::max(std::abs(result_old(i, j)), 1e-8);  // avoid divide-by-zero
      double rel_change = std::abs(result_new(i, j) - result_old(i, j)) / denom;
      residual = std::max(residual, rel_change);
    }

    residual = GlobalMax(m_grid->com, residual);

    // copy into `result`, updating ghosts for the next iteration
    result.copy_from(result_new);


    if (residual < rtol) {
      m_log->message(2, "grounding line elevation converged iteration %03d, max rel. change = %f\n", iter, residual);
      break;
    }
    if (iter == max_iter) {
      m_log->message(2, "grounding line elevation maximum number of iterations reached %03d, max rel. change = %f\n", max_iter, residual);
    }
  }
}

void Picop::compute_shelf_base_elevation(const Inputs &inputs,
                                         array::Scalar1 &result) {

  const auto &cell_type = inputs.geometry->cell_type;
  const auto &ice_surface_elevation = inputs.geometry->ice_surface_elevation;
  const auto &ice_thickness = inputs.geometry->ice_thickness;

  array::AccessScope scope{&cell_type,
                           &ice_surface_elevation, &ice_thickness, &result };

  for (auto p = m_grid->points(); p; p.next()) {
    int i = p.i(), j = p.j();      
      result(i, j) = ice_surface_elevation(i, j) - ice_thickness(i, j);
  }
  result.update_ghosts();
}


void Picop::compute_grounding_line_slope(const Inputs &inputs,
                                         array::Scalar1 &result) {

  const auto &cell_type = inputs.geometry->cell_type;
  const auto &ice_surface_elevation = inputs.geometry->ice_surface_elevation;
  const auto &ice_thickness = inputs.geometry->ice_thickness;

  const auto &zb = m_shelf_base_elevation;
  
  array::AccessScope scope{&zb, &ice_surface_elevation, &ice_thickness, &cell_type, &result };
  
  const double
    dx = m_grid->dx(),
    dy = m_grid->dy();

  using mask::floating_ice;
  using mask::ice_free_ocean;

  bool cfbc = m_config->get_flag("stress_balance.calving_front_stress_bc");

  auto diff_grounded = diff_centered;
  if (m_config->get_flag("stress_balance.ssa.fd.upstream_surface_slope_approximation")) {
    diff_grounded = diff_uphill;
  }
  
  for (auto p = m_grid->points(); p; p.next()) {
    const int i = p.i(), j = p.j(); {

    // To compute the x-derivative we use
    //
    // * away from the grounding line, ice margins, and no_model mask transitions -- 2nd
    //   order centered difference
    //
    // * at the grounded cell near the grounding line -- 1st order
    //   one-sided difference using the grounded neighbor
    //
    // * at the floating cell near the grounding line -- 1st order
    //   one-sided difference using the floating neighbor
    //
    // All these cases can be combined by writing h_x as the weighted
    // average of one-sided differences, with weights of 0 if a finite
    // difference is not used and 1 if it is.
    //
    // The y derivative is handled the same way.

    auto M = cell_type.star_int(i, j);
    auto h = zb.star(i, j);

    // x-derivative
    double h_x = 0.0;
    {
      int west = weight(cfbc, M.c, M.w, h.c, h.w),
          east = weight(cfbc, M.c, M.e, h.c, h.e);

      if (east + west == 2 and mask::grounded_ice(M.c)) {
        // interior of the ice blob: use the "uphill-biased" difference
        h_x = diff_grounded(h.w, h.c, h.e) / dx;
      } else if (east + west > 0) {
        h_x = 1.0 / ((west + east) * dx) * (west * (h.c - h.w) + east * (h.e - h.c));
        if (floating_ice(M.c) and (ice_free_ocean(M.e) or ice_free_ocean(M.w))) {
          // at the ice front: use constant extrapolation to approximate the value outside
          // the ice extent (see the notes in the manual)
          h_x /= 2.0;
        }
      } else {
        h_x = 0.0;
      }
    }

    // y-derivative
    double h_y = 0.0;
    {
      int south = weight(cfbc, M.c, M.s, h.c, h.s),
          north = weight(cfbc, M.c, M.n, h.c, h.n);

      if (north + south == 2 and mask::grounded_ice(M.c)) {
        // interior of the ice blob: use the "uphill-biased" difference
        h_y = diff_grounded(h.s, h.c, h.n) / dy;
      } else if (north + south > 0) {
        h_y = 1.0 / ((south + north) * dy) * (south * (h.c - h.s) + north * (h.n - h.c));
        if (floating_ice(M.c) and (ice_free_ocean(M.s) or ice_free_ocean(M.n))) {
          // at the ice front: use constant extrapolation to approximate the value outside
          // the ice extent
          h_y /= 2.0;
        }
      } else {
        h_y = 0.0;
      }
    }
    
    double slope =  atan(sqrt(h_x*h_x + h_y*h_y));
    if (slope >= M_PI) {
      slope = M_PI - 0.001;
    }
    result(i, j) = slope;
    }
  }
  
  const double rtol = 0.001;
  const int max_iter = 500;
  
  auto &result_old = result;
  auto &result_new = m_work;

  const auto &adv_vel   = inputs.stress_balance->advective_velocity();

  scope.add({ &result_old, &result_new, &adv_vel });

  for (int iter = 0; iter < max_iter; ++iter) {

    transport_step(result_old, cell_type, m_flow_direction, result_new);
    
    double residual = 0.0;
    for (auto p = m_grid->points(); p; p.next()) {
      int i = p.i(), j = p.j();
      double denom = std::max(std::abs(result_old(i, j)), 1e-8);  // avoid divide-by-zero
      double rel_change = std::abs(result_new(i, j) - result_old(i, j)) / denom;
      residual = std::max(residual, rel_change);
    }

    residual = GlobalMax(m_grid->com, residual);

    // copy into `result`, updating ghosts for the next iteration
    result.copy_from(result_new);

    if (residual < rtol) {
      m_log->message(2, "grounding line slope converged iteration %03d, max rel. change = %f\n", iter, residual);
      break;
    }
    if (iter == max_iter) {
      m_log->message(2, "grounding line slope maximum number of iterations reached %03d, max rel .change = %f\n", max_iter, residual);
    }
  
  }

}
// Write diagnostic variables to extra files if requested
DiagnosticList Picop::diagnostics_impl() const {

  DiagnosticList result = {
    { "picop_grounding_line_elevation", Diagnostic::wrap(m_grounding_line_elevation) },
    { "picop_basal_melt_rate", Diagnostic::wrap(m_basal_melt_rate) },
    { "picop_temperature", Diagnostic::wrap(m_theta_ocean) },
    { "picop_salinity", Diagnostic::wrap(m_salinity_ocean) },
    { "picop_shelf_base_elevation", Diagnostic::wrap(m_shelf_base_elevation) },
    { "picop_grounding_line_slope", Diagnostic::wrap(m_grounding_line_slope) },
    { "picop_geometric_scale", Diagnostic::wrap(m_geometric_scale) },
    { "picop_length_scale", Diagnostic::wrap(m_length_scale) },
  };

  return combine(result, OceanModel::diagnostics_impl());
}


} // end of namespace ocean
} // end of namespace pism
