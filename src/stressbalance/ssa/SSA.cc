// Copyright (C) 2004--2019, 2021 Constantine Khroulev, Ed Bueler, Jed Brown, Torsten Albrecht
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
#include "pism/util/IceModelVec2CellType.hh"
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


SSA::SSA(IceGrid::ConstPtr g)
  : ShallowStressBalance(g),
    m_mask(m_grid, "ssa_mask", WITH_GHOSTS, m_config->get_number("grid.max_stencil_width")),
    m_taud(m_grid, "taud", WITHOUT_GHOSTS),
    m_velocity_global(m_grid, "bar", WITHOUT_GHOSTS)
{

  m_e_factor = m_config->get_number("stress_balance.ssa.enhancement_factor");

  strength_extension = new SSAStrengthExtension(*m_config);

  // grounded_dragging_floating integer mask
  m_mask.set_attrs("diagnostic", "ice-type (ice-free/grounded/floating/ocean) integer mask",
                   "", "", "", 0);
  std::vector<double> mask_values = {MASK_ICE_FREE_BEDROCK, MASK_GROUNDED, MASK_FLOATING, MASK_ICE_FREE_OCEAN};
  m_mask.metadata().set_numbers("flag_values", mask_values);
  m_mask.metadata().set_string("flag_meanings",
                              "ice_free_bedrock grounded_ice floating_ice ice_free_ocean");

  m_taud.set_attrs("diagnostic",
                   "X-component of the driving shear stress at the base of ice",
                   "Pa", "Pa", "", 0);
  m_taud.set_attrs("diagnostic",
                   "Y-component of the driving shear stress at the base of ice",
                   "Pa", "Pa", "", 1);

  // override velocity metadata
  m_velocity.metadata(0).set_name("u_ssa");
  m_velocity.metadata(0).set_string("long_name", "SSA model ice velocity in the X direction");

  m_velocity.metadata(1).set_name("v_ssa");
  m_velocity.metadata(1).set_string("long_name", "SSA model ice velocity in the Y direction");

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

If configuration parameter `sia.surface_gradient_method` = `eta` then the surface
gradient \f$\nabla h\f$ is computed by the gradient of the transformed variable
\f$\eta= H^{(2n+2)/n}\f$ (frequently, \f$\eta= H^{8/3}\f$). The idea is that
this quantity is more regular at ice sheet margins, and so we get a better
surface gradient. When the thickness at a grid point is very small (below \c
minThickEtaTransform in the procedure), the formula is slightly modified to
give a lower driving stress. The transformation is not used in floating ice.
 */
void SSA::compute_driving_stress(const IceModelVec2S &ice_thickness,
                                 const IceModelVec2S &surface_elevation,
                                 const IceModelVec2CellType &cell_type,
                                 const IceModelVec2Int *no_model_mask,
                                 IceModelVec2V &result) const {

  bool cfbc = m_config->get_flag("stress_balance.calving_front_stress_bc");
  bool surface_gradient_inward = m_config->get_flag("stress_balance.ssa.compute_surface_gradient_inward");

  double
    dx = m_grid->dx(),
    dy = m_grid->dy();

  IceModelVec::AccessList list{&surface_elevation, &cell_type, &ice_thickness, &result};

  if (no_model_mask) {
    list.add(*no_model_mask);
  }

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    const double pressure = m_EC->pressure(ice_thickness(i, j)); // FIXME issue #15
    if (pressure <= 0.0) {
      result(i, j) = 0.0;
      continue;
    }

    // Special case for verification tests.
    if (surface_gradient_inward) {
      double
        h_x = surface_elevation.diff_x_p(i, j),
        h_y = surface_elevation.diff_y_p(i, j);
      result(i, j) = - pressure * Vector2(h_x, h_y);
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

    auto M = cell_type.int_star(i, j);
    auto h = surface_elevation.star(i, j);
    StarStencil<int> N(0);

    if (no_model_mask) {
      N = no_model_mask->int_star(i, j);
    }

    // x-derivative
    double h_x = 0.0;
    {
      double
        west = weight(cfbc, M.ij, M.w, h.ij, h.w, N.ij, N.w),
        east = weight(cfbc, M.ij, M.e, h.ij, h.e, N.ij, N.e);

      if (east + west > 0) {
        h_x = 1.0 / ((west + east) * dx) * (west * (h.ij - h.w) + east * (h.e - h.ij));
      } else {
        h_x = 0.0;
      }
    }

    // y-derivative
    double h_y = 0.0;
    {
      double
        south = weight(cfbc, M.ij, M.s, h.ij, h.s, N.ij, N.s),
        north = weight(cfbc, M.ij, M.n, h.ij, h.n, N.ij, N.n);

      if (north + south > 0) {
        h_y = 1.0 / ((south + north) * dy) * (south * (h.ij - h.s) + north * (h.n - h.ij));
      } else {
        h_y = 0.0;
      }
    }

    result(i, j) = - pressure * Vector2(h_x, h_y);
  }
}

std::string SSA::stdout_report() const {
  return m_stdout_ssa;
}


//! \brief Set the initial guess of the SSA velocity.
void SSA::set_initial_guess(const IceModelVec2V &guess) {
  m_velocity.copy_from(guess);
}

const IceModelVec2V& SSA::driving_stress() const {
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
    v.set_string("comment",
                 "this is the driving stress used by the SSA solver");
  }
}

IceModelVec::Ptr SSA_taud::compute_impl() const {

  IceModelVec2V::Ptr result(new IceModelVec2V(m_grid, "result", WITHOUT_GHOSTS));
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
  m_vars[0].set_string("comment",
                     "this is the magnitude of the driving stress used by the SSA solver");
}

IceModelVec::Ptr SSA_taud_mag::compute_impl() const {

  // Allocate memory:
  IceModelVec2S::Ptr result(new IceModelVec2S(m_grid, "taud_mag", WITHOUT_GHOSTS));
  result->metadata() = m_vars[0];

  result->set_to_magnitude(model->driving_stress());

  return result;
}

} // end of namespace stressbalance
} // end of namespace pism
