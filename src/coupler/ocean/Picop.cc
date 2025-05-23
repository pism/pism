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

#include <gsl/gsl_math.h> // GSL_NAN

#include "pism/coupler/util/options.hh"
#include "pism/geometry/Geometry.hh"
#include "pism/geometry/UNO.hh"
#include "pism/util/Config.hh"
#include "pism/util/Grid.hh"
#include "pism/util/Mask.hh"
#include "pism/util/Time.hh"
#include "pism/stressbalance/StressBalance.hh"

#include "pism/coupler/ocean/PicoPhysics.hh"
#include "pism/coupler/ocean/Picop.hh"
#include "pism/coupler/ocean/PicopPhysics.hh"
#include "pism/util/array/Forcing.hh"
#include "pism/util/Logger.hh"

namespace pism {

namespace ocean {

Picop::Picop(std::shared_ptr<const Grid> grid)
  : CompleteOceanModel(grid, std::shared_ptr<OceanModel>()),
    m_pico(std::make_shared<Pico>(grid)),
    m_uno(new pism::UNO(grid, pism::PISM_UNO_UPWIND1)),
    m_basal_melt_rate(m_grid, "picop_basal_melt_rate"),
    m_grounding_line_elevation(grid, "picop_grounding_line_elevation"),
    m_shelf_base_elevation(grid, "picop_shelf_base_elevation"),
    m_slope(grid, "picop_basal_slope"),
    m_theta_ocean(m_pico->get_temperature()),
    m_salinity_ocean(m_pico->get_salinity()),
    m_flow_direction(grid, "ice_flow_direction")
{

  ForcingOptions opt(*m_grid->ctx(), "ocean.picop");

  m_grounding_line_elevation.metadata(0).long_name("grounding line elevation").units("m");
  m_grounding_line_elevation.metadata()["_FillValue"] = { 0.0 };
  m_grounding_line_elevation.set(0.0);

  m_shelf_base_elevation.metadata(0).long_name("shelf base elevation").units("m");
  m_shelf_base_elevation.metadata()["_FillValue"] = { 0.0 };
  m_shelf_base_elevation.set(0.0);

}

void Picop::init_impl(const Geometry &geometry) {
  (void) geometry;

  m_pico->init(geometry);
  m_log->message(2, "* Initializing the Plume extension of PICO for the ocean ...\n");

  PicoPhysics pico_physics(*m_config);
  PicopPhysics picop_physics(*m_config);
    
  // compute_shelf_base_elevation(geometry, m_shelf_base_elevation);
  // compute_grounding_line_elevation(geometry, m_grounding_line_elevation);
  // compute_slope(geometry, m_shelf_base_elevation, m_slope);
  // compute_melt_rate(picop_physics, m_theta_ocean, m_salinity_ocean,
  //                   m_grounding_line_elevation, m_shelf_base_elevation,
  //                   m_slope, m_basal_melt_rate);

  // m_shelf_base_mass_flux->copy_from(m_basal_melt_rate);
  // m_shelf_base_mass_flux->scale(pico_physics.ice_density());

  double
    ice_density   = m_config->get_number("constants.ice.density"),
    water_density = m_config->get_number("constants.sea_water.density"),
    g             = m_config->get_number("constants.standard_gravity");

  compute_average_water_column_pressure(geometry, ice_density, water_density, g,
                                        *m_water_column_pressure);
  
}

void Picop::define_model_state_impl(const File &output) const {


  OceanModel::define_model_state_impl(output);
}

void Picop::write_model_state_impl(const File &output) const {


  OceanModel::write_model_state_impl(output);
}

void Picop::update_impl(const Inputs &inputs, double t, double dt) {

  m_pico->update(inputs, t, dt);

  const auto &geometry = *inputs.geometry;
  
  PicoPhysics pico_physics(*m_config);
  PicopPhysics picop_physics(*m_config);

  compute_shelf_base_elevation(geometry, m_shelf_base_elevation);
  
  compute_grounding_line_elevation(inputs, m_grounding_line_elevation);
  
  compute_slope(geometry, m_shelf_base_elevation, m_slope);
  
  compute_melt_rate(picop_physics, m_theta_ocean, m_salinity_ocean,
                    m_grounding_line_elevation, m_shelf_base_elevation,
                    m_slope, m_basal_melt_rate);

  m_shelf_base_mass_flux->copy_from(m_basal_melt_rate);
  m_shelf_base_mass_flux->scale(pico_physics.ice_density());

  double
    ice_density   = m_config->get_number("constants.ice.density"),
    water_density = m_config->get_number("constants.sea_water.density"),
    g             = m_config->get_number("constants.standard_gravity");

  compute_average_water_column_pressure(geometry, ice_density, water_density, g,
                                        *m_water_column_pressure);

}


MaxTimestep Picop::max_timestep_impl(double t) const {
  (void) t;

  return MaxTimestep("ocean picop");
}

void Picop::compute_melt_rate(const PicopPhysics &physics,
                              const array::Scalar &T_a,
                              const array::Scalar &S_a,
                              array::Scalar1 &grounding_line_elevation,
                              array::Scalar1 &shelf_base_elevation,
                              array::Scalar1 &slope,
                              array::Scalar1 &melt_rate) const {

  array::AccessScope scope{&grounding_line_elevation,
                           &T_a, &S_a,
                           &shelf_base_elevation,
                           &slope, &melt_rate};

  for (auto p = m_grid->points(); p; p.next()) {
    int i = p.i(), j = p.j();

    const double z_b = shelf_base_elevation(i, j);
    const double z_gl = grounding_line_elevation(i, j);
    const double alpha = slope(i, j);
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

void Picop::compute_shelf_base_elevation(const Geometry &geometry,
                                         array::Scalar1 &shelf_base_elevation) const {
  
  const array::Scalar &surface    = geometry.ice_surface_elevation;
  const array::Scalar &H          = geometry.ice_thickness;
  const array::Scalar &cell_type  = geometry.cell_type;

  array::AccessScope scope{&surface, &H, &cell_type, &shelf_base_elevation};

  // Step 1: Initialize zgl0 at grounding line: bed elevation
  for (auto p = m_grid->points(); p; p.next()) {
    int i = p.i(), j = p.j();
    shelf_base_elevation(i, j) = (cell_type(i, j) == MASK_FLOATING) ? surface(i, j) - H(i, j) : 0.0;
  }
}

void Picop::compute_grounding_line_elevation(const Inputs &inputs,
                                             array::Scalar1 &grounding_line_elevation) {

  const array::Scalar &bed           = inputs.geometry->bed_elevation;
  const array::Scalar &H             = inputs.geometry->ice_thickness;
  const array::CellType1 &cell_type  = inputs.geometry->cell_type;
  const array::Scalar &z_s           = inputs.geometry->sea_level_elevation;
  const array::Vector &adv_vel       = inputs.stress_balance->advective_velocity();;

  array::AccessScope scope{&bed, &H, &z_s, &cell_type, &grounding_line_elevation};

  // compute_magnitude(adv_vel, m_work);
  
  // Step 1: Initialize zgl0 at grounding line: bed elevation
  //         Normalize velocities
  for (auto p = m_grid->points(); p; p.next()) {
    int i = p.i(), j = p.j();
    grounding_line_elevation(i, j) = (cell_type(i, j) == MASK_GROUNDED) ? bed(i, j) : 0.0;
    double flow_speed = adv_vel(i,j).magnitude();
    if (flow_speed > 0.0) {
      m_flow_direction(i, j) = adv_vel(i, j) / flow_speed;
    } else {
      m_flow_direction(i, j) = 0.0;
    }
  }


  // const double tol = 1e-4;
  // const double alpha = 0.1;
  // const double max_iter = 500;
  
  // using std::sqrt;
  
  // for (int iter = 0; iter < max_iter; ++iter) {

  //   double residual = 0.0;

  //   const array::Scalar &result_old = m_uno->x();
  //   m_uno->update(alpha, cell_type, grounding_line_elevation, adv_vel, true);
  //   const array::Scalar &result_new = m_uno->x();

  //   for (auto p = m_grid->points(); p; p.next()) {
  //     int i = p.i(), j = p.j();
  //     residual += result_new(i, j) * result_new(i, j) - result_old(i, j) * result_old(i, j);
  //   }

  //   residual = std::sqrt(GlobalSum(m_grid->com, residual));
  //   if (residual < tol) break;
  // }
    
}

void Picop::compute_slope(const Geometry &geometry,
                          array::Scalar1 &shelf_base_elevation,
                          array::Scalar &slope) const {

  const array::CellType1 &cell_type = geometry.cell_type;
  
  array::AccessScope scope{&cell_type, &shelf_base_elevation};

  double
    dx = m_grid->dx(),
    dy = m_grid->dy();

  
  for (auto p = m_grid->points(); p; p.next()) {
    int i = p.i(), j = p.j();
    auto h = shelf_base_elevation.star(i, j);
    double s_n = (h.c - h.n) / (dx * dx + dy * dy);

    slope(i, j) = s_n;
  }
  

  
}

// Write diagnostic variables to extra files if requested
DiagnosticList Picop::diagnostics_impl() const {

  DiagnosticList result = {
    { "picop_grounding_line_elevation", Diagnostic::wrap(m_grounding_line_elevation) },
  };

  return combine(result, OceanModel::diagnostics_impl());
}


} // end of namespace ocean
} // end of namespace pism
