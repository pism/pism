// Copyright (C) 2004--2019, 2021, 2022, 2023 Constantine Khroulev, Ed Bueler, Jed Brown, Torsten Albrecht
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

#include "SSA.hh"
#include "pism/basalstrength/basal_resistance.hh"
#include "pism/util/EnthalpyConverter.hh"
#include "pism/rheology/FlowLawFactory.hh"
#include "pism/util/Mask.hh"
#include "pism/util/Vars.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/io/File.hh"
#include "pism/util/pism_options.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/util/array/CellType.hh"
#include "pism/stressbalance/StressBalance.hh"
#include "pism/geometry/Geometry.hh"

#include "SSA_diagnostics.hh"

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
    m_mask(m_grid, "ssa_mask"),
    m_taud(m_grid, "taud"),
    m_velocity_global(m_grid, "bar")
{

  m_e_factor = m_config->get_number("stress_balance.ssa.enhancement_factor");

  strength_extension = new SSAStrengthExtension(*m_config);

  // grounded_dragging_floating integer mask
  m_mask.set_attrs("diagnostic", "ice-type (ice-free/grounded/floating/ocean) integer mask",
                   "", "", "", 0);
  m_mask.metadata()["flag_values"] =
    {MASK_ICE_FREE_BEDROCK, MASK_GROUNDED, MASK_FLOATING, MASK_ICE_FREE_OCEAN};
  m_mask.metadata()["flag_meanings"] =
    "ice_free_bedrock grounded_ice floating_ice ice_free_ocean";

  m_taud.set_attrs("diagnostic",
                   "X-component of the driving shear stress at the base of ice",
                   "Pa", "Pa", "", 0);
  m_taud.set_attrs("diagnostic",
                   "Y-component of the driving shear stress at the base of ice",
                   "Pa", "Pa", "", 1);

  // override velocity metadata
  m_velocity.metadata(0).set_name("u_ssa");
  m_velocity.metadata(0)["long_name"] = "SSA model ice velocity in the X direction";

  m_velocity.metadata(1).set_name("v_ssa");
  m_velocity.metadata(1)["long_name"] = "SSA model ice velocity in the Y direction";

  m_da = m_velocity_global.dm();

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
      File input_file(m_grid->com, opts.filename, PISM_GUESS, PISM_READONLY);
      bool u_ssa_found = input_file.find_variable("u_ssa");
      bool v_ssa_found = input_file.find_variable("v_ssa");
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

  // update the cell type mask using the ice-free thickness threshold for stress balance
  // computations
  {
    const double H_threshold = m_config->get_number("stress_balance.ice_free_thickness_standard");
    GeometryCalculator gc(*m_config);
    gc.set_icefree_thickness(H_threshold);

    gc.compute_mask(inputs.geometry->sea_level_elevation,
                    inputs.geometry->bed_elevation,
                    inputs.geometry->ice_thickness,
                    m_mask);
  }

  if (full_update) {
    solve(inputs);
    compute_basal_frictional_heating(m_velocity,
                                     *inputs.basal_yield_stress,
                                     m_mask,
                                     m_basal_frictional_heating);
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
static int weight(bool margin_bc,
                  int M_ij, int M_n,
                  double h_ij, double h_n,
                  int N_ij, int N_n) {
  using mask::grounded;
  using mask::icy;
  using mask::floating_ice;
  using mask::ice_free;
  using mask::ice_free_ocean;

  // grounding lines and calving fronts
  if ((grounded(M_ij) and floating_ice(M_n)) or
      (floating_ice(M_ij) and grounded(M_n)) or
      (floating_ice(M_ij) and ice_free_ocean(M_n))) {
    return 0;
  }

  // fjord walls, nunataks, headwalls
  if ((icy(M_ij) and ice_free(M_n) and h_n > h_ij) or
      (ice_free(M_ij) and icy(M_n) and h_ij > h_n)) {
    return 0;
  }

  // This condition has to match the one used to implement the calving front stress
  // boundary condition in SSAFD::assemble_rhs().
  if (margin_bc and
      ((icy(M_ij) and ice_free(M_n)) or
       (ice_free(M_ij) and icy(M_n)))) {
    return 0;
  }

  // boundaries of the "no model" area
  if ((N_ij == 0 and N_n == 1) or (N_ij == 1 and N_n == 0)) {
    return 0;
  }

  return 1;
}

//! \brief Compute the gravitational driving stress.
/*!
Computes the gravitational driving stress at the base of the ice:
\f[ \tau_d = - \rho g H \nabla h \f]
 */
void SSA::compute_driving_stress(const array::Scalar &ice_thickness,
                                 const array::Scalar1 &surface_elevation,
                                 const array::CellType1 &cell_type,
                                 const array::Scalar1 *no_model_mask,
                                 array::Vector &result) const {

  using mask::ice_free_ocean;
  using mask::floating_ice;

  bool cfbc = m_config->get_flag("stress_balance.calving_front_stress_bc");
  bool surface_gradient_inward = m_config->get_flag("stress_balance.ssa.compute_surface_gradient_inward");

  double
    dx = m_grid->dx(),
    dy = m_grid->dy();

  array::AccessScope list{&surface_elevation, &cell_type, &ice_thickness, &result};

  if (no_model_mask) {
    list.add(*no_model_mask);
  }

  for (auto p = m_grid->points(); p; p.next()) {
    const int i = p.i(), j = p.j();

    const double pressure = m_EC->pressure(ice_thickness(i, j)); // FIXME issue #15
    if (pressure <= 0.0) {
      result(i, j) = 0.0;
      continue;
    }

    // Special case for verification tests.
    if (surface_gradient_inward) {
      double
        h_x = diff_x_p(surface_elevation, i, j),
        h_y = diff_y_p(surface_elevation, i, j);
      result(i, j) = - pressure * Vector2d(h_x, h_y);
      continue;
    }

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

    auto M = cell_type.star(i, j);
    auto h = surface_elevation.star(i, j);
    stencils::Star<int> N(0);

    if (no_model_mask) {
      N = no_model_mask->star_int(i, j);
    }

    // x-derivative
    double h_x = 0.0;
    {
      double
        west = weight(cfbc, M.c, M.w, h.c, h.w, N.c, N.w),
        east = weight(cfbc, M.c, M.e, h.c, h.e, N.c, N.e);

      if (east + west > 0) {
        h_x = 1.0 / ((west + east) * dx) * (west * (h.c - h.w) + east * (h.e - h.c));
        if (floating_ice(M.c) and (ice_free_ocean(M.e) or ice_free_ocean(M.w)))  {
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
      double
        south = weight(cfbc, M.c, M.s, h.c, h.s, N.c, N.s),
        north = weight(cfbc, M.c, M.n, h.c, h.n, N.c, N.n);

      if (north + south > 0) {
        h_y = 1.0 / ((south + north) * dy) * (south * (h.c - h.s) + north * (h.n - h.c));
        if (floating_ice(M.c) and (ice_free_ocean(M.s) or ice_free_ocean(M.n)))  {
          // at the ice front: use constant extrapolation to approximate the value outside
          // the ice extent
          h_y /= 2.0;
        }
      } else {
        h_y = 0.0;
      }
    }

    result(i, j) = - pressure * Vector2d(h_x, h_y);
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

      auto M = cell_type.star(i, j);
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


//! \brief Set the initial guess of the SSA velocity.
void SSA::set_initial_guess(const array::Vector &guess) {
  m_velocity.copy_from(guess);
}

const array::Vector& SSA::driving_stress() const {
  return m_taud;
}


void SSA::define_model_state_impl(const File &output) const {
  m_velocity.define(output);
}

void SSA::write_model_state_impl(const File &output) const {
  m_velocity.write(output);
}

DiagnosticList SSA::diagnostics_impl() const {
  DiagnosticList result = ShallowStressBalance::diagnostics_impl();

  // replace these diagnostics
  result["taud"] = Diagnostic::Ptr(new SSA_taud(this));
  result["taud_mag"] = Diagnostic::Ptr(new SSA_taud_mag(this));

  return result;
}

SSA_taud::SSA_taud(const SSA *m)
  : Diag<SSA>(m) {

  // set metadata:
  m_vars = {SpatialVariableMetadata(m_sys, "taud_x"),
            SpatialVariableMetadata(m_sys, "taud_y")};

  set_attrs("X-component of the driving shear stress at the base of ice", "",
            "Pa", "Pa", 0);
  set_attrs("Y-component of the driving shear stress at the base of ice", "",
            "Pa", "Pa", 1);

  for (auto &v : m_vars) {
    v["comment"] = "this is the driving stress used by the SSA solver";
  }
}

array::Array::Ptr SSA_taud::compute_impl() const {

  array::Vector::Ptr result(new array::Vector(m_grid, "result"));
  result->metadata(0) = m_vars[0];
  result->metadata(1) = m_vars[1];

  result->copy_from(model->driving_stress());

  return result;
}

SSA_taud_mag::SSA_taud_mag(const SSA *m)
  : Diag<SSA>(m) {

  // set metadata:
  m_vars = {SpatialVariableMetadata(m_sys, "taud_mag")};

  set_attrs("magnitude of the driving shear stress at the base of ice", "",
            "Pa", "Pa", 0);
  m_vars[0]["comment"] =
    "this is the magnitude of the driving stress used by the SSA solver";
}

array::Array::Ptr SSA_taud_mag::compute_impl() const {

  // Allocate memory:
  array::Scalar::Ptr result(new array::Scalar(m_grid, "taud_mag"));
  result->metadata() = m_vars[0];

  compute_magnitude(model->driving_stress(), *result);

  return result;
}

} // end of namespace stressbalance
} // end of namespace pism
