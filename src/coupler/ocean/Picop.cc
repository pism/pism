// Copyright (C) 2012-2019, 2021, 2022, 2023, 2024, 2025, 2026 Constantine Khrulev, Ricarda Winkelmann, Ronja Reese, Torsten
// Albrecht, Matthias Mengel, and Andy Aschwanden
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
#include <vector>
#include <algorithm>
#include <mpi.h>

#include "pism/coupler/util/options.hh"
#include "pism/geometry/Geometry.hh"
#include "pism/util/Config.hh"
#include "pism/util/Grid.hh"
#include "pism/util/pism_utilities.hh"  // GlobalSum
#include "pism/stressbalance/StressBalance.hh"
#include "pism/hydrology/Hydrology.hh"

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
    m_local_slope(grid, "picop_local_slope"),
    m_fresh_water_melt_rate(grid, "picop_fresh_water_melt_rate"),
    m_discharge_flux(grid, "picop_discharge_flux"),
    m_disch_q0(grid, "picop_disch_q0"),
    m_disch_L5(grid, "picop_disch_L5"),
    m_disch_s(grid, "picop_disch_s"),
    m_theta_ocean(m_pico->get_temperature()),
    m_salinity_ocean(m_pico->get_salinity()),
    m_flow_direction(grid, "ice_flow_direction"),
    m_work(grid, "temporary_storage"),
    m_zb_x(grid, "staggered slope x"),
    m_zb_y(grid, "staggered slope y") {

  ForcingOptions opt(*m_grid->ctx(), "ocean.picop");

  m_add_fresh_water_melt = m_config->get_flag("ocean.picop.add_fresh_water_melt");
  {
    const auto method = m_config->get_string("ocean.picop.discharge_method");
    if (method == "along_flow") {
      m_discharge_method = DISCHARGE_ALONG_FLOW;
    } else if (method == "downstream_gate") {
      m_discharge_method = DISCHARGE_DOWNSTREAM_GATE;
    } else {
      m_discharge_method = DISCHARGE_ISOTROPIC;
    }
  }

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
  
  m_local_slope.metadata(0).long_name("shelf base slope").units("rad").output_units("degree");
  m_local_slope.metadata()["_FillValue"] = { 0.0 };
  m_local_slope.set(0.0);
      
  m_fresh_water_melt_rate.metadata(0)
      .long_name("PICOP fresh water melt rate")
      .units("m s^-1")
      .output_units("m year^-1");
  m_fresh_water_melt_rate.metadata()["_FillValue"] = {0.0};
  m_fresh_water_melt_rate.set(0.0);

  m_discharge_flux.metadata(0)
      .long_name("PICOP subglacial discharge flux on floating ice")
      .units("m^2 s^-1");
  m_discharge_flux.metadata()["_FillValue"] = {0.0};
  m_discharge_flux.set(0.0);

  // Internal along-flow transport tracers, exposed as diagnostics for debugging.
  m_disch_q0.metadata(0).long_name("PICOP transported source discharge q_sg0").units("m^2 s^-1");
  m_disch_L5.metadata(0).long_name("PICOP transported governing length scale 5L'").units("m");
  m_disch_s.metadata(0).long_name("PICOP along-flow distance from discharge outflow").units("m");
  m_disch_q0.set(0.0);
  m_disch_L5.set(0.0);
  m_disch_s.set(0.0);

  m_shelf_base_temperature->metadata()["_FillValue"] = {0.0};
}

void Picop::init_impl(const Geometry &geometry) {
  (void) geometry;

  m_pico->init(geometry);
  m_log->message(2, "* Initializing the Plume extension of PICO (PICOP) for the ocean ...\n");
  m_log->message(2, "  Note: PICOP requires stress balance computation to be enabled.\n");

  PicopPhysics picop_physics(*m_config);

  double
    ice_density   = m_config->get_number("constants.ice.density"),
    water_density = m_config->get_number("constants.sea_water.density"),
    g             = m_config->get_number("constants.standard_gravity");

  compute_average_water_column_pressure(geometry, ice_density, water_density, g,
                                        *m_water_column_pressure);
}

std::set<VariableMetadata> Picop::state_impl() const {
  return m_pico->state();
}

void Picop::write_state_impl(const OutputFile &output) const {
  m_pico->write_state(output);
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

  for (auto p : grid->points()) {

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
    m_log->message(3,
                   "WARNING: PICOP requires stress balance for plume transport calculations.\n"
                   "         Stress balance not available - falling back to PICO melt rates.\n"
                   "         To use PICOP, enable stress balance computation.\n");
    m_shelf_base_temperature->copy_from(m_pico->shelf_base_temperature());
    m_shelf_base_mass_flux->copy_from(m_pico->shelf_base_mass_flux());
    return;
  }

  if (inputs.hydrology == nullptr) {
    // Use outputs from PICOP if the stress balance is not available
    m_log->message(3,
                   "WARNING: PICOP requires hydrology routing for plume transport calculations.\n"
                   );
  }

  m_log->message(3, "  PICOP: Computing plume-based melt rates...\n");
  
  PicopPhysics picop_physics(*m_config);

  compute_shelf_base_elevation(inputs, m_shelf_base_elevation);

  profiling().begin("ocean.compute_grounding_line_elevation");
  compute_grounding_line_elevation(inputs, m_grounding_line_elevation);
  profiling().end("ocean.compute_grounding_line_elevation");
  
  profiling().begin("ocean.compute_local_slope");
  compute_local_slope(inputs, m_local_slope);
  profiling().end("ocean.compute_local_slope");
  
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

namespace {
//! A grounding-line subglacial-discharge outflow.
struct Outflow {
  double x, y;         // location (m)
  double q_sg0;        // discharge flux at the outflow (m^2 s^-1)
  double radius;       // governing length scale 5L' (m)
  double dir_x, dir_y; // normalized ice-flow direction at the outflow (downstream)
};

//! Gather every rank's outflow list onto all ranks (MPI_Allgatherv).
/*! The number of outflows is O(grounding-line length / dx) -- a few hundred -- so this
    is cheap and lets the q_sg painting below avoid wide-stencil communication. */
std::vector<Outflow> gather_outflows(MPI_Comm com, const std::vector<Outflow> &local) {
  int size = 1;
  MPI_Comm_size(com, &size);

  const int n_local = static_cast<int>(local.size());
  std::vector<int> counts(size);
  MPI_Allgather(&n_local, 1, MPI_INT, counts.data(), 1, MPI_INT, com);

  // 6 doubles per outflow
  std::vector<int> dcounts(size), ddispls(size);
  int total = 0;
  for (int r = 0; r < size; ++r) {
    dcounts[r] = counts[r] * 6;
    ddispls[r] = total;
    total += dcounts[r];
  }

  std::vector<double> sendbuf(static_cast<size_t>(n_local) * 6);
  for (int k = 0; k < n_local; ++k) {
    sendbuf[6 * k + 0] = local[k].x;
    sendbuf[6 * k + 1] = local[k].y;
    sendbuf[6 * k + 2] = local[k].q_sg0;
    sendbuf[6 * k + 3] = local[k].radius;
    sendbuf[6 * k + 4] = local[k].dir_x;
    sendbuf[6 * k + 5] = local[k].dir_y;
  }

  std::vector<double> recvbuf(static_cast<size_t>(total));
  MPI_Allgatherv(sendbuf.data(), n_local * 6, MPI_DOUBLE,
                 recvbuf.data(), dcounts.data(), ddispls.data(), MPI_DOUBLE, com);

  std::vector<Outflow> all(recvbuf.size() / 6);
  for (size_t k = 0; k < all.size(); ++k) {
    all[k] = { recvbuf[6 * k + 0], recvbuf[6 * k + 1],
               recvbuf[6 * k + 2], recvbuf[6 * k + 3],
               recvbuf[6 * k + 4], recvbuf[6 * k + 5] };
  }
  return all;
}
} // end of anonymous namespace

// Semi-Lagrangian transport helpers (defined below, used by build_discharge_field).
void transport_step(const array::Scalar1 &U_old, const array::CellType &cell_type,
                    const array::Vector &flow_direction, array::Scalar &U_new);
void transport_step_distance(const array::Scalar1 &U_old, const array::CellType &cell_type,
                             const array::Vector &flow_direction, array::Scalar &U_new);

//! Build the discharge flux field q_sg(x,y) on floating cells.
/*!
 * PISM routes subglacial water under grounded ice only, so the discharge enters the
 * ocean cavity at the grounding line. We (1) find grounding-line outflows and their flux
 * q_sg0 (PISM's hydrology flux() is already m^2 s^-1, so no channel-width conversion is
 * needed), (2) compute each outflow's governing length scale 5L' = 5 q_sg0 / m_fw, and
 * (3) paint q_sg onto floating cells within 5L' with quadratic decay to zero at 5L',
 * taking the max where outflows overlap (Pelle et al. 2023, Defining q_sg(x,y)).
 */
void Picop::build_discharge_field(const Inputs &inputs,
                                  const PicopPhysics &physics,
                                  const array::Scalar &T_a,
                                  const array::Scalar &S_a) {

  // Averaging radius for outflow fields: Pelle et al. (2023) average alpha, z_gl, T_a,
  // S_a within 5 km of the discharge outflow before computing the length scale 5L'.
  const double avg_radius = 5000.0;  // m

  const auto &cell_type = inputs.geometry->cell_type;
  const auto &surf      = inputs.geometry->ice_surface_elevation;
  const auto &thk       = inputs.geometry->ice_thickness;

  array::Scalar &water_flux = m_work;
  compute_magnitude(inputs.hydrology->flux(), water_flux);

  // Pass 1: grounding-line discharge outflows (location + flux). The length scale
  // (radius) is filled in below from 5-km-averaged fields.
  std::vector<Outflow> local;
  {
    array::AccessScope scope{&cell_type, &water_flux, &m_flow_direction};
    for (auto p : m_grid->points()) {
      const int i = p.i(), j = p.j();

      // Outflow = grounded-ice cell at the grounding line with nonzero flux.
      if (not (cell_type.grounded_ice(i, j) and cell_type.next_to_floating_ice(i, j))) {
        continue;
      }
      const double q0 = water_flux(i, j);  // already m^2 s^-1
      if (q0 <= 0.0) {
        continue;
      }
      // Ice-flow direction at the outflow (downstream, into the cavity). Populated by
      // compute_grounding_line_elevation(), which runs before compute_melt_rate().
      const auto D = m_flow_direction(i, j);
      local.push_back({ m_grid->x(i), m_grid->y(j), q0, 0.0, D.u, D.v });
    }
  }

  std::vector<Outflow> all = gather_outflows(m_grid->com, local);
  const int N = static_cast<int>(all.size());

  if (N == 0) {
    m_log->message(2,
                   "PICOP discharge: no outflows found"
                   " (no grounded cell next to floating ice with positive subglacial water flux)\n");
    m_discharge_flux.set(0.0);
    return;
  }

  // Pass 1.5: accumulate floating-cell fields within avg_radius of each outflow, summed
  // across ranks so the average spans rank boundaries. Buffer layout, per outflow k:
  // [0*N+k] alpha, [1*N+k] z_gl, [2*N+k] t_a, [3*N+k] s_a, [4*N+k] count.
  std::vector<double> local_acc(5 * N, 0.0), acc(5 * N, 0.0);
  {
    array::AccessScope scope{&cell_type, &surf, &thk,
                             &m_grounding_line_elevation, &m_local_slope, &T_a, &S_a};
    for (auto p : m_grid->points()) {
      const int i = p.i(), j = p.j();
      if (not cell_type.floating_ice(i, j)) {
        continue;
      }
      const double x = m_grid->x(i), y = m_grid->y(j);
      const double z_b   = surf(i, j) - thk(i, j);
      const double z_gl  = std::min(m_grounding_line_elevation(i, j), z_b);
      const double alpha = m_local_slope(i, j);
      const double s_a   = S_a(i, j);
      double t_a = T_a(i, j);
      const double t_min = physics.characteristic_freezing_point(s_a, 0);
      if (t_a < t_min) {
        t_a = t_min;
      }
      for (int k = 0; k < N; ++k) {
        if (std::hypot(x - all[k].x, y - all[k].y) < avg_radius) {
          local_acc[0 * N + k] += alpha;
          local_acc[1 * N + k] += z_gl;
          local_acc[2 * N + k] += t_a;
          local_acc[3 * N + k] += s_a;
          local_acc[4 * N + k] += 1.0;
        }
      }
    }
  }
  GlobalSum(m_grid->com, local_acc.data(), acc.data(), 5 * N);

  // Compute each outflow's governing length scale 5L' from the averaged fields.
  for (int k = 0; k < N; ++k) {
    const double n = acc[4 * N + k];
    if (n <= 0.0) {
      all[k].radius = 0.0;  // no floating cell within 5 km (shouldn't happen at the GL)
      continue;
    }
    const double alpha = acc[0 * N + k] / n;
    const double z_gl  = acc[1 * N + k] / n;
    const double t_a   = acc[2 * N + k] / n;
    const double s_a   = acc[3 * N + k] / n;

    const double t_f_gl   = physics.characteristic_freezing_point(s_a, z_gl);
    const double Gamma_TS = physics.effective_heat_exchange_coefficient(t_a, t_f_gl, alpha);
    const double m_fw0    = physics.fresh_water_melt_rate(all[k].q_sg0, s_a, t_a,
                                                          Gamma_TS, t_f_gl, alpha);
    all[k].radius = physics.governing_length_scale(all[k].q_sg0, m_fw0);  // 5L'
  }

  {
    double max_q0 = 0.0, max_radius = 0.0;
    int n_active = 0;
    for (const auto &o : all) {
      max_q0     = std::max(max_q0, o.q_sg0);
      max_radius = std::max(max_radius, o.radius);
      if (o.radius > 0.0) {
        n_active += 1;
      }
    }
    m_log->message(3,
                   "PICOP discharge: %d outflow(s), %d with 5L' > 0;"
                   " max q_sg0 = %.3e m2/s, max 5L' = %.1f m\n",
                   N, n_active, max_q0, max_radius);
  }

  // Pass 2: distribute q_sg onto floating cells (method selected by ocean.picop.discharge_method).
  m_discharge_flux.set(0.0);

  // Along-flow method: transport the discharge downstream along ice flow instead of painting
  // an isotropic disk. Three tracers are relaxed over the floating region from the outflows:
  //   q0 = source discharge q_sg0, L5 = source 5L', s = along-flow path distance from source.
  // Then q_sg = q0 * (1 - s/L5)^2 (clamped), i.e. the same quadratic decay as the isotropic
  // paint but with distance measured along the flowline. Each floating cell is fed by its
  // single upstream outflow (the semi-Lagrangian characteristic), not a max over all disks.
  if (m_discharge_method == DISCHARGE_ALONG_FLOW) {
    // Dirichlet BCs: zero everywhere except the outflow cells (held constant by transport_step
    // since they are grounded). Locate each outflow's owned cell by inverting its coordinates.
    m_disch_q0.set(0.0);
    m_disch_L5.set(0.0);
    m_disch_s.set(0.0);

    const double x0 = m_grid->x(0), y0 = m_grid->y(0);
    const double dx = m_grid->dx(), dy = m_grid->dy();
    const int xs = m_grid->xs(), xm = m_grid->xm();
    const int ys = m_grid->ys(), ym = m_grid->ym();
    {
      array::AccessScope scope{&m_disch_q0, &m_disch_L5, &m_disch_s};
      for (const auto &o : all) {
        if (o.radius <= 0.0) {
          continue;
        }
        const int i = static_cast<int>(std::lround((o.x - x0) / dx));
        const int j = static_cast<int>(std::lround((o.y - y0) / dy));
        if (i >= xs and i < xs + xm and j >= ys and j < ys + ym) {
          m_disch_q0(i, j) = o.q_sg0;
          m_disch_L5(i, j) = o.radius;
          m_disch_s(i, j)  = 0.0;
        }
      }
    }
    m_disch_q0.update_ghosts();
    m_disch_L5.update_ghosts();
    m_disch_s.update_ghosts();

    // Relaxation (same structure as compute_grounding_line_elevation): converge on q0, the
    // bounded tracer; L5 and s propagate along the same characteristics.
    const int max_iter = 500;
    const double rtol  = 1e-3;
    array::Scalar &scratch = m_work;
    double residual = 0.0;
    for (int iter = 0; iter < max_iter; ++iter) {
      transport_step(m_disch_q0, cell_type, m_flow_direction, scratch);

      residual = 0.0;
      {
        array::AccessScope scope{&m_disch_q0, &scratch};
        for (auto p : m_grid->points()) {
          const int i = p.i(), j = p.j();
          const double denom = std::max(std::abs(m_disch_q0(i, j)), 1e-12);
          residual = std::max(residual, std::abs(scratch(i, j) - m_disch_q0(i, j)) / denom);
        }
      }
      residual = GlobalMax(m_grid->com, residual);
      m_disch_q0.copy_from(scratch);

      transport_step(m_disch_L5, cell_type, m_flow_direction, scratch);
      m_disch_L5.copy_from(scratch);

      transport_step_distance(m_disch_s, cell_type, m_flow_direction, scratch);
      m_disch_s.copy_from(scratch);

      if (residual < rtol) {
        m_log->message(3,
                       "PICOP discharge along-flow transport converged at iteration %03d"
                       " (max rel. change %f)\n", iter, residual);
        break;
      }
    }
    if (residual >= rtol) {
      m_log->message(2,
                     "PICOP discharge along-flow transport reached max iterations %03d"
                     " (max rel. change %f)\n", max_iter, residual);
    }

    // Assemble q_sg from the transported tracers, with quadratic decay in along-flow distance.
    // Sub-grid floor (mirrors the isotropic method's r_nn floor): when 5L' is smaller than the
    // grid spacing the decay would zero the plume at the very first cell, so floating cells
    // within r_nn along-flow distance of an outflow receive the full discharge instead.
    const double r_nn = 1.5 * std::max(dx, dy);
    array::AccessScope scope{&cell_type, &m_disch_q0, &m_disch_L5, &m_disch_s, &m_discharge_flux};
    for (auto p : m_grid->points()) {
      const int i = p.i(), j = p.j();
      if (not cell_type.floating_ice(i, j)) {
        continue;
      }
      const double L5 = m_disch_L5(i, j);
      if (L5 > 0.0) {
        double w;
        if (m_disch_s(i, j) < r_nn) {
          w = 1.0;  // sub-grid floor: first cell(s) along flow get the full discharge
        } else {
          w = std::max(0.0, std::min(1.0, 1.0 - m_disch_s(i, j) / L5));
        }
        m_discharge_flux(i, j) = m_disch_q0(i, j) * w * w;
      }
    }
    return;
  }

  // Isotropic method (optionally restricted to downstream cells): paint q_sg on floating cells,
  // taking the max over overlapping outflows. A cell within 5L' of an outflow gets the
  // discharge flux with quadratic decay (0 at 5L').
  //
  // When the plume is sub-grid (5L' < grid spacing) the decay reaches no floating cell center,
  // so we additionally deposit the full q_sg0 into the outflow's immediate floating neighbor(s)
  // (within r_nn): at this resolution the plume is confined to the first cell.
  const bool gate = (m_discharge_method == DISCHARGE_DOWNSTREAM_GATE);
  const double r_nn = 1.5 * std::max(m_grid->dx(), m_grid->dy());
  {
    array::AccessScope scope{&cell_type, &m_discharge_flux};
    for (auto p : m_grid->points()) {
      const int i = p.i(), j = p.j();
      if (not cell_type.floating_ice(i, j)) {
        continue;
      }
      const double x = m_grid->x(i), y = m_grid->y(j);
      double q = 0.0;
      for (const auto &o : all) {
        if (o.radius <= 0.0) {
          continue;
        }
        // Directional gate: paint only floating cells downstream of the outflow, i.e. where
        // (cell - outflow) . flow_direction > 0. This confines the plume influence to the
        // along-flow direction instead of an isotropic disk (which lit up the shelf sides).
        // If the outflow has no flow direction, fall back to the isotropic behavior.
        if (gate) {
          const double dir_mag2 = o.dir_x * o.dir_x + o.dir_y * o.dir_y;
          if (dir_mag2 > 0.0) {
            const double dot = (x - o.x) * o.dir_x + (y - o.y) * o.dir_y;
            if (dot <= 0.0) {
              continue;
            }
          }
        }
        const double d = std::hypot(x - o.x, y - o.y);
        double w2 = 0.0;
        if (d < o.radius) {
          const double w = 1.0 - d / o.radius;  // 1 at the outflow, 0 at 5L'
          w2 = w * w;                            // quadratic decay
        }
        if (d < r_nn) {
          w2 = 1.0;  // sub-grid floor: nearest floating cell(s) receive the full flux
        }
        if (w2 > 0.0) {
          q = std::max(q, o.q_sg0 * w2);
        }
      }
      m_discharge_flux(i, j) = q;
    }
  }
}

void Picop::compute_melt_rate(const Inputs &inputs,
                              const PicopPhysics &physics,
                              const array::Scalar &T_a,
                              const array::Scalar &S_a,
                              array::Scalar1 &result)  {

  const auto &ice_surface_elevation = inputs.geometry->ice_surface_elevation;
  const auto &ice_thickness = inputs.geometry->ice_thickness;
  const auto &cell_type = inputs.geometry->cell_type;

  // Build q_sg(x,y) from grounding-line discharge outflows (fills m_discharge_flux).
  build_discharge_field(inputs, physics, T_a, S_a);

  array::AccessScope scope{&T_a, &S_a, &cell_type, &m_discharge_flux,
                           &ice_surface_elevation, &ice_thickness,
                           &m_local_slope, &m_grounding_line_elevation,
                           &m_fresh_water_melt_rate, &result};

  for (auto p : m_grid->points()) {
    int i = p.i(), j = p.j();

    if (cell_type.floating_ice(i, j)) {
      const double z_b = ice_surface_elevation(i, j) - ice_thickness(i, j);
      // Clamp z_gl: grounding line origin cannot be above the local shelf base
      // (same as ISSM: if(zgl > z_base) zgl = z_base)
      const double z_gl = std::min(m_grounding_line_elevation(i, j), z_b);
      const double alpha = m_local_slope(i, j);

      const double s_a = S_a(i, j);
      double t_a = T_a(i, j);
      const double t_min = physics.characteristic_freezing_point(s_a, 0);
      /* Low bound for Toc to ensure X_hat is between 0 and 1 */
      if (t_a < t_min) {
        t_a = t_min;
      }
      const double t_f_gl = physics.characteristic_freezing_point(s_a, z_gl);
      const double Gamma_TS = physics.effective_heat_exchange_coefficient(t_a, t_f_gl, alpha);
      const double l = physics.length_scaling(t_a, t_f_gl, Gamma_TS, alpha);
      const double g_alpha = physics.geometric_scaling(Gamma_TS, alpha);
      const double X_hat = std::min(std::max(physics.dimensionless_coordinate(z_b, z_gl, l), 0.0), 1.0);
      const double M_X = physics.dimensionless_melt_curve(X_hat);
      const double M = physics.melt_function(t_a, t_f_gl, g_alpha);
      const double q_sg = m_discharge_flux(i, j);
      // Fresh-water (subglacial discharge plume) melt contribution, optionally disabled
      // via ocean.picop.add_fresh_water_melt. When disabled, the melt rate reduces to the
      // ambient (PICO-style) contribution M_X * M.
      const double m_fw =
          m_add_fresh_water_melt
              ? physics.fresh_water_melt_rate(q_sg, s_a, t_a, Gamma_TS, t_f_gl, alpha)
              : 0.0;
      m_fresh_water_melt_rate(i, j) = m_fw;
      // Equation 16 in Pelle et al (2023). Guard the 0/0 that occurs when both the
      // ambient melt M and the discharge melt m_fw vanish (e.g. a flat shelf base with
      // no subglacial discharge); the limit of the ratio in that case is 0.
      const double denom = M + m_fw;
      result(i, j) = denom > 0.0 ? (M_X * M * M) / denom + m_fw : m_fw;
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

  for (auto p : grid->points()) {
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

/*!
 * Like transport_step(), but accumulates along-flow path distance: on floating ice the
 * transported value is incremented by the physical length of one semi-Lagrangian step
 * (the departure point is one flow-vector away from the cell center). Used to advect the
 * along-flow distance `s` from each discharge outflow. Non-floating cells are held constant.
 */
void transport_step_distance(const array::Scalar1 &U_old, const array::CellType &cell_type,
                             const array::Vector &flow_direction, array::Scalar &U_new) {

  auto grid = U_new.grid();
  const double dx = grid->dx(), dy = grid->dy();

  array::AccessScope scope{ &U_old, &cell_type, &U_new, &flow_direction };

  for (auto p : grid->points()) {
    const int i = p.i(), j = p.j();

    if (cell_type.floating_ice(i, j)) {
      auto U = U_old.box(i, j);
      auto D = flow_direction(i, j);

      // |D| == 1 (normalized), so the step covers one flow-vector; its physical length is
      // hypot(D.u*dx, D.v*dy) (== dx for square cells).
      const double ds = std::hypot(D.u * dx, D.v * dy);

      U_new(i, j) = interpolate(U, -D.u, -D.v) + ds;
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
  for (auto p : m_grid->points()) {
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

  double residual = 0.0;
  for (int iter = 0; iter < max_iter; ++iter) {

    transport_step(result_old, cell_type, m_flow_direction, result_new);

    // elevation of a plume origin that reached (x,y) cannot be above the elevation of the
    // bottom of the ice at (x,y)
    for (auto p : m_grid->points()) {
      int i = p.i(), j = p.j();

      if (cell_type.floating_ice(i, j)) {
        double ice_bottom_elevation = ice_surface_elevation(i, j) - ice_thickness(i, j);

        if (result_new(i, j) > ice_bottom_elevation) {
          result_new(i, j) = ice_bottom_elevation;
        }
      }
    }

    residual = 0.0;
    for (auto p : m_grid->points()) {
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
  }

  // Warn if maximum iterations reached without convergence
  if (residual >= rtol) {
    m_log->message(2, "grounding line elevation maximum number of iterations reached %03d, max rel. change = %f\n", max_iter, residual);
  }
}

void Picop::compute_shelf_base_elevation(const Inputs &inputs,
                                         array::Scalar1 &result) {

  const auto &cell_type = inputs.geometry->cell_type;
  const auto &ice_surface_elevation = inputs.geometry->ice_surface_elevation;
  const auto &ice_thickness = inputs.geometry->ice_thickness;

  array::AccessScope scope{&cell_type,
                           &ice_surface_elevation, &ice_thickness, &result };

  for (auto p : m_grid->points()) {
    int i = p.i(), j = p.j();      
      result(i, j) = ice_surface_elevation(i, j) - ice_thickness(i, j);
  }
  result.update_ghosts();
}


void Picop::compute_local_slope(const Inputs &inputs,
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
  
  for (auto p : m_grid->points()) {
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
    
    double slope = atan(sqrt(h_x*h_x + h_y*h_y));
    result(i, j) = slope;
    }
  }
}
// Write diagnostic variables to extra files if requested
DiagnosticList Picop::spatial_diagnostics_impl() const {

  DiagnosticList result = {
    { "picop_basal_melt_rate", Diagnostic::wrap(m_basal_melt_rate) },
    { "picop_fresh_water_melt_rate", Diagnostic::wrap(m_fresh_water_melt_rate) },
    { "picop_discharge_flux", Diagnostic::wrap(m_discharge_flux) },
    { "picop_disch_q0", Diagnostic::wrap(m_disch_q0) },
    { "picop_disch_L5", Diagnostic::wrap(m_disch_L5) },
    { "picop_disch_s", Diagnostic::wrap(m_disch_s) },
    { "picop_grounding_line_elevation", Diagnostic::wrap(m_grounding_line_elevation) },
    { "picop_local_slope", Diagnostic::wrap(m_local_slope) },
    { "picop_temperature", Diagnostic::wrap(m_theta_ocean) },
    { "picop_salinity", Diagnostic::wrap(m_salinity_ocean) },
    { "picop_shelf_base_elevation", Diagnostic::wrap(m_shelf_base_elevation) },
  };

  return combine(result, OceanModel::spatial_diagnostics_impl());
}


} // end of namespace ocean
} // end of namespace pism
