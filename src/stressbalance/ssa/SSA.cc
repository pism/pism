// Copyright (C) 2004--2019, 2021, 2022, 2023, 2024, 2025 Constantine Khroulev, Ed Bueler, Jed Brown, Torsten Albrecht
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

#include "pism/stressbalance/ssa/SSA.hh"
#include "pism/util/EnthalpyConverter.hh"
#include "pism/rheology/FlowLawFactory.hh"
#include "pism/util/Mask.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/io/File.hh"
#include "pism/util/array/CellType.hh"
#include "pism/stressbalance/StressBalance.hh"
#include "pism/geometry/Geometry.hh"
#include "pism/util/Logger.hh"

namespace pism {
namespace stressbalance {

SSAStrengthExtension::SSAStrengthExtension(const Config &config) {
  m_min_thickness = config.get_number("stress_balance.ssa.strength_extension.min_thickness");
  m_constant_nu = config.get_number("stress_balance.ssa.strength_extension.constant_nu");
}

//! Set strength = (viscosity times thickness).
/*! Determines nu by input strength and current min_thickness. */
void SSAStrengthExtension::set_notional_strength(double my_nuH) {
  if (my_nuH <= 0.0) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "nuH must be positive, got %f", my_nuH);
  }
  m_constant_nu = my_nuH / m_min_thickness;
}

//! Set minimum thickness to trigger use of extension.
/*! Preserves strength (nuH) by also updating using current nu.  */
void SSAStrengthExtension::set_min_thickness(double my_min_thickness) {
  if (my_min_thickness <= 0.0) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "min_thickness must be positive, got %f",
                                  my_min_thickness);
  }
  double nuH = m_constant_nu * m_min_thickness;
  m_min_thickness = my_min_thickness;
  m_constant_nu = nuH / m_min_thickness;
}

//! Returns strength = (viscosity times thickness).
double SSAStrengthExtension::get_notional_strength() const {
  return m_constant_nu * m_min_thickness;
}

//! Returns minimum thickness to trigger use of extension.
double SSAStrengthExtension::get_min_thickness() const {
  return m_min_thickness;
}


SSA::SSA(std::shared_ptr<const Grid> g)
  : ShallowStressBalance(g),
    m_velocity_global(m_grid, "bar")
{

  m_e_factor = m_config->get_number("stress_balance.ssa.enhancement_factor");

  strength_extension = new SSAStrengthExtension(*m_config);

  // override velocity metadata
  m_velocity.metadata(0).set_name("u_ssa");
  m_velocity.metadata(0).long_name("SSA model ice velocity in the X direction");

  m_velocity.metadata(1).set_name("v_ssa");
  m_velocity.metadata(1).long_name("SSA model ice velocity in the Y direction");

  {
    rheology::FlowLawFactory ice_factory("stress_balance.ssa.", m_config, m_EC);
    ice_factory.remove(ICE_GOLDSBY_KOHLSTEDT);
    m_flow_law = ice_factory.create();
  }
}

SSA::~SSA() {
  if (strength_extension != NULL) {
    delete strength_extension;
    strength_extension = NULL;
  }
}


//! \brief Initialize a generic regular-grid SSA solver.
void SSA::init_impl() {

  ShallowStressBalance::init_impl();

  m_log->message(2, "* Initializing the SSA stress balance...\n");
  m_log->message(2,
             "  [using the %s flow law]\n", m_flow_law->name().c_str());

  InputOptions opts = process_input_options(m_grid->com, m_config);

  // Check if PISM is being initialized from an output file from a previous run
  // and read the initial guess (unless asked not to).
  if (opts.type == INIT_RESTART) {
    if (m_config->get_flag("stress_balance.ssa.read_initial_guess")) {
      File input_file(m_grid->com, opts.filename, io::PISM_GUESS, io::PISM_READONLY);
      bool u_ssa_found = input_file.variable_exists("u_ssa");
      bool v_ssa_found = input_file.variable_exists("v_ssa");
      unsigned int start = input_file.nrecords() - 1;

      if (u_ssa_found and v_ssa_found) {
        m_log->message(3, "Reading u_ssa and v_ssa...\n");

        m_velocity.read(input_file, start);
      }
    }
  } else {
    m_velocity.set(0.0); // default initial guess
  }
}

//! \brief Update the SSA solution.
void SSA::update(const Inputs &inputs, bool full_update) {

  if (full_update) {
    solve(inputs);
    compute_basal_frictional_heating(m_velocity,
                                     *inputs.basal_yield_stress,
                                     inputs.geometry->cell_type,
                                     m_basal_frictional_heating);
  }
}


/*!
 * Estimate velocity at ice-free cells near the ice margin using interpolation from
 * immediate neighbors that are icy.
 *
 * This is used to improve the initial guess of ice viscosity at marginal locations when
 * ice advances: otherwise we would use the *zero* velocity (if CFBC is "on"), and that is
 * a poor estimate at best.
 *
 * Note: icy cells of `velocity` are treated as read-only, and ice-free marginal cells are
 * write-only. This means that it's okay for `velocity` to be a input-output argument: we
 * don't use of the values modified by this method.
 */
void SSA::extrapolate_velocity(const array::CellType1 &cell_type,
                               array::Vector1 &velocity) const {
  array::AccessScope list{&cell_type, &velocity};

  for (auto p = m_grid->points(); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (cell_type.ice_free(i, j) and cell_type.next_to_ice(i, j)) {

      auto M = cell_type.star_int(i, j);
      auto vel = velocity.star(i, j);

      Vector2d sum{0.0, 0.0};
      int N = 0;
      for (auto d : {North, East, South, West}) {
        if (mask::icy(M[d])) {
          sum += vel[d];
          ++N;
        }
      }
      velocity(i, j) = sum / std::max(N, 1);
    }
  }
  velocity.update_ghosts();
}


std::string SSA::stdout_report() const {
  return m_stdout_ssa;
}

void SSA::define_model_state_impl(const File &output) const {
  m_velocity.define(output, io::PISM_DOUBLE);
}

void SSA::write_model_state_impl(const File &output) const {
  m_velocity.write(output);
}

} // end of namespace stressbalance
} // end of namespace pism
