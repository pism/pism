// Copyright (C) 2004--2016 Constantine Khroulev, Ed Bueler, Jed Brown, Torsten Albrecht
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
#include "base/basalstrength/basal_resistance.hh"
#include "base/enthalpyConverter.hh"
#include "base/rheology/FlowLawFactory.hh"
#include "base/util/Mask.hh"
#include "base/util/PISMVars.hh"
#include "base/util/error_handling.hh"
#include "base/util/io/PIO.hh"
#include "base/util/pism_options.hh"
#include "base/util/pism_utilities.hh"
#include "base/util/IceModelVec2CellType.hh"

#include "SSA_diagnostics.hh"

namespace pism {
namespace stressbalance {

SSAStrengthExtension::SSAStrengthExtension(const Config &config) {
  m_min_thickness = config.get_double("min_thickness_strength_extension_ssa");
  m_constant_nu = config.get_double("constant_nu_strength_extension_ssa");
}

//! Set strength = (viscosity times thickness).
/*! Determines nu by input strength and current min_thickness. */
void SSAStrengthExtension::set_notional_strength(double my_nuH) {
  if (my_nuH <= 0.0) {
    throw RuntimeError::formatted("nuH must be positive, got %f", my_nuH);
  }
  m_constant_nu = my_nuH / m_min_thickness;
}

//! Set minimum thickness to trigger use of extension.
/*! Preserves strength (nuH) by also updating using current nu.  */
void SSAStrengthExtension::set_min_thickness(double my_min_thickness) {
  if (my_min_thickness <= 0.0) {
    throw RuntimeError::formatted("min_thickness must be positive, got %f",
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


SSA::SSA(IceGrid::ConstPtr g, EnthalpyConverter::Ptr e)
  : ShallowStressBalance(g, e)
{
  m_thickness = NULL;
  m_tauc = NULL;
  m_surface = NULL;
  m_bed = NULL;
  m_ice_enthalpy = NULL;
  m_gl_mask = NULL;

  strength_extension = new SSAStrengthExtension(*m_config);

  const unsigned int WIDE_STENCIL = m_config->get_double("grid_max_stencil_width");

  // grounded_dragging_floating integer mask
  m_mask.create(m_grid, "ssa_mask", WITH_GHOSTS, WIDE_STENCIL);
  m_mask.set_attrs("diagnostic", "ice-type (ice-free/grounded/floating/ocean) integer mask",
                  "", "");
  std::vector<double> mask_values(4);
  mask_values[0] = MASK_ICE_FREE_BEDROCK;
  mask_values[1] = MASK_GROUNDED;
  mask_values[2] = MASK_FLOATING;
  mask_values[3] = MASK_ICE_FREE_OCEAN;
  m_mask.metadata().set_doubles("flag_values", mask_values);
  m_mask.metadata().set_string("flag_meanings",
                              "ice_free_bedrock grounded_ice floating_ice ice_free_ocean");

  m_taud.create(m_grid, "taud", WITHOUT_GHOSTS);
  m_taud.set_attrs("diagnostic",
                 "X-component of the driving shear stress at the base of ice",
                 "Pa", "", 0);
  m_taud.set_attrs("diagnostic",
                 "Y-component of the driving shear stress at the base of ice",
                 "Pa", "", 1);


  // override velocity metadata
  m_velocity.metadata(0).set_name("u_ssa");
  m_velocity.metadata(0).set_string("long_name", "SSA model ice velocity in the X direction");

  m_velocity.metadata(1).set_name("v_ssa");
  m_velocity.metadata(1).set_string("long_name", "SSA model ice velocity in the Y direction");

  m_velocity_global.create(m_grid, "bar", WITHOUT_GHOSTS);

  m_da = m_velocity_global.get_dm();

  {
    rheology::FlowLawFactory ice_factory("ssa_", m_config, m_EC);
    ice_factory.remove(ICE_GOLDSBY_KOHLSTEDT);
    m_flow_law = ice_factory.create();
  }
}

SSA::~SSA() {
  if (m_flow_law != NULL) {
    delete m_flow_law;
    m_flow_law = NULL;
  }
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

  if (m_config->get_boolean("sub_groundingline")) {
    m_gl_mask = m_grid->variables().get_2d_scalar("gl_mask");
  }

  m_thickness = m_grid->variables().get_2d_scalar("land_ice_thickness");
  m_tauc      = m_grid->variables().get_2d_scalar("tauc");

  try {
    m_surface = m_grid->variables().get_2d_scalar("surface_altitude");
  } catch (RuntimeError) {
    m_surface = NULL;
  }

  m_bed      = m_grid->variables().get_2d_scalar("bedrock_altitude");
  m_ice_enthalpy = m_grid->variables().get_3d_scalar("enthalpy");

  // Check if PISM is being initialized from an output file from a previous run
  // and read the initial guess (unless asked not to).
  options::String input_file("-i", "PISM input file");
  bool bootstrap = options::Bool("-bootstrap", "enable bootstrapping heuristics");

  if (input_file.is_set() and not bootstrap) {
    bool u_ssa_found = false, v_ssa_found = false;
    unsigned int start = 0;
    PIO nc(m_grid->com, "guess_mode");

    bool dont_read_initial_guess = options::Bool("-dontreadSSAvels",
                                                 "don't read the initial guess");

    nc.open(input_file, PISM_READONLY);
    u_ssa_found = nc.inq_var("u_ssa");
    v_ssa_found = nc.inq_var("v_ssa");
    start = nc.inq_nrecords() - 1;
    nc.close();

    if (u_ssa_found && v_ssa_found && (not dont_read_initial_guess)) {
      m_log->message(3, "Reading u_ssa and v_ssa...\n");

      m_velocity.read(input_file, start);
    }

  } else {
    m_velocity.set(0.0); // default initial guess
  }

  if (m_config->get_boolean("ssa_dirichlet_bc")) {
    m_bc_mask = m_grid->variables().get_2d_mask("bc_mask");
    m_bc_values = m_grid->variables().get_2d_vector("vel_ssa_bc");
  }
}

//! \brief Update the SSA solution.
void SSA::update(bool fast, double sea_level, const IceModelVec2S &melange_back_pressure) {

  m_sea_level = sea_level;
  (void) melange_back_pressure;

  // update the cell type mask using the ice-free thickness threshold for stress balance
  // computations
  {
    const double H_threshold = m_config->get_double("mask_icefree_thickness_stress_balance_standard");
    GeometryCalculator gc(*m_config);
    gc.set_icefree_thickness(H_threshold);

    gc.compute_mask(sea_level, *m_bed, *m_thickness, m_mask);
  }

  if (not fast) {
    solve();
    compute_basal_frictional_heating(m_velocity, *m_tauc, m_mask,
                                     m_basal_frictional_heating);
  }
}

//! \brief Compute the gravitational driving stress.
/*!
Computes the gravitational driving stress at the base of the ice:
\f[ \tau_d = - \rho g H \nabla h \f]

If configuration parameter `surface_gradient_method` = `eta` then the surface
gradient \f$\nabla h\f$ is computed by the gradient of the transformed variable
\f$\eta= H^{(2n+2)/n}\f$ (frequently, \f$\eta= H^{8/3}\f$). The idea is that
this quantity is more regular at ice sheet margins, and so we get a better
surface gradient. When the thickness at a grid point is very small (below \c
minThickEtaTransform in the procedure), the formula is slightly modified to
give a lower driving stress. The transformation is not used in floating ice.
 */
void SSA::compute_driving_stress(IceModelVec2V &result) {
  const IceModelVec2S &thk = *m_thickness; // to improve readability (below)

  const double n = m_flow_law->exponent(), // frequently n = 3
    etapow  = (2.0 * n + 2.0)/n,  // = 8/3 if n = 3
    invpow  = 1.0 / etapow,  // = 3/8
    dinvpow = (- n - 2.0) / (2.0 * n + 2.0); // = -5/8
  const double minThickEtaTransform = 5.0; // m
  const double dx=m_grid->dx(), dy=m_grid->dy();

  bool cfbc = m_config->get_boolean("calving_front_stress_boundary_condition");
  bool compute_surf_grad_inward_ssa = m_config->get_boolean("compute_surf_grad_inward_ssa");
  bool use_eta = (m_config->get_string("surface_gradient_method") == "eta");

  IceModelVec::AccessList list;
  list.add(*m_surface);
  list.add(*m_bed);
  list.add(m_mask);
  list.add(thk);
  list.add(result);

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    const double pressure = m_EC->pressure(thk(i,j)); // FIXME issue #15
    if (pressure <= 0.0) {
      result(i,j).u = 0.0;
      result(i,j).v = 0.0;
    } else {
      double h_x = 0.0, h_y = 0.0;
      // FIXME: we need to handle grid periodicity correctly.
      if (m_mask.grounded(i,j) && (use_eta == true)) {
        // in grounded case, differentiate eta = H^{8/3} by chain rule
        if (thk(i,j) > 0.0) {
          const double myH = (thk(i,j) < minThickEtaTransform ?
                              minThickEtaTransform : thk(i,j));
          const double eta = pow(myH, etapow), factor = invpow * pow(eta, dinvpow);
          h_x = factor * (pow(thk(i+1,j),etapow) - pow(thk(i-1,j),etapow)) / (2*dx);
          h_y = factor * (pow(thk(i,j+1),etapow) - pow(thk(i,j-1),etapow)) / (2*dy);
        }
        // now add bed slope to get actual h_x,h_y
        // FIXME: there is no reason to assume user's bed is periodized
        h_x += m_bed->diff_x(i,j);
        h_y += m_bed->diff_y(i,j);
      } else {  // floating or eta transformation is not used
        if (compute_surf_grad_inward_ssa) {
          // Special case for verification tests.
          h_x = m_surface->diff_x_p(i,j);
          h_y = m_surface->diff_y_p(i,j);
        } else {              // general case

          // To compute the x-derivative we use
          // * away from the grounding line -- 2nd order centered difference
          //
          // * at the grounded cell near the grounding line -- 1st order
          //   one-sided difference using the grounded neighbor
          //
          // * at the floating cell near the grounding line -- 1st order
          //   one-sided difference using the floating neighbor
          //
          // All three cases can be combined by writing h_x as the weighted
          // average of one-sided differences, with weights of 0 if a finite
          // difference is not used and 1 if it is.
          //
          // The y derivative is handled the same way.

          // x-derivative
          {
            double west = 1, east = 1;
            if ((m_mask.grounded(i,j) && m_mask.floating_ice(i+1,j)) ||
                (m_mask.floating_ice(i,j) && m_mask.grounded(i+1,j)) ||
                (m_mask.floating_ice(i,j) && m_mask.ice_free_ocean(i+1,j))) {
              east = 0;
            }
            if ((m_mask.grounded(i,j) && m_mask.floating_ice(i-1,j)) ||
                (m_mask.floating_ice(i,j) && m_mask.grounded(i-1,j)) ||
                (m_mask.floating_ice(i,j) && m_mask.ice_free_ocean(i-1,j))) {
              west = 0;
            }

            // This driving stress computation has to match the calving front
            // stress boundary condition in SSAFD::assemble_rhs().
            if (cfbc) {
              if (m_mask.icy(i,j) && m_mask.ice_free(i+1,j)) {
                east = 0;
              }
              if (m_mask.icy(i,j) && m_mask.ice_free(i-1,j)) {
                west = 0;
              }
            }

            if (east + west > 0) {
              h_x = 1.0 / (west + east) * (west * m_surface->diff_x_stagE(i-1,j) +
                                           east * m_surface->diff_x_stagE(i,j));
            } else {
              h_x = 0.0;
            }
          }

          // y-derivative
          {
            double south = 1, north = 1;
            if ((m_mask.grounded(i,j) && m_mask.floating_ice(i,j+1)) ||
                (m_mask.floating_ice(i,j) && m_mask.grounded(i,j+1)) ||
                (m_mask.floating_ice(i,j) && m_mask.ice_free_ocean(i,j+1))) {
              north = 0;
            }
            if ((m_mask.grounded(i,j) && m_mask.floating_ice(i,j-1)) ||
                (m_mask.floating_ice(i,j) && m_mask.grounded(i,j-1)) ||
                (m_mask.floating_ice(i,j) && m_mask.ice_free_ocean(i,j-1))) {
              south = 0;
            }

            // This driving stress computation has to match the calving front
            // stress boundary condition in SSAFD::assemble_rhs().
            if (cfbc) {
              if (m_mask.icy(i,j) && m_mask.ice_free(i,j+1)) {
                north = 0;
              }
              if (m_mask.icy(i,j) && m_mask.ice_free(i,j-1)) {
                south = 0;
              }
            }

            if (north + south > 0) {
              h_y = 1.0 / (south + north) * (south * m_surface->diff_y_stagN(i,j-1) +
                                             north * m_surface->diff_y_stagN(i,j));
            } else {
              h_y = 0.0;
            }
          }

        } // end of "general case"

      } // end of "floating or eta transformation is not used"

      result(i,j).u = - pressure * h_x;
      result(i,j).v = - pressure * h_y;
    } // end of "(pressure > 0)"
  }
}

std::string SSA::stdout_report() {
  return m_stdout_ssa;
}


//! \brief Set the initial guess of the SSA velocity.
void SSA::set_initial_guess(const IceModelVec2V &guess) {
  m_velocity.copy_from(guess);
}


void SSA::add_vars_to_output_impl(const std::string &/*keyword*/, std::set<std::string> &result) {
  result.insert("vel_ssa");
}


void SSA::define_variables_impl(const std::set<std::string> &vars, const PIO &nc, IO_Type nctype) {

  if (set_contains(vars, "vel_ssa")) {
    m_velocity.define(nc, nctype);
  }
}


void SSA::write_variables_impl(const std::set<std::string> &vars, const PIO &nc) {

  if (set_contains(vars, "vel_ssa")) {
    m_velocity.write(nc);
  }
}

void SSA::get_diagnostics_impl(std::map<std::string, Diagnostic::Ptr> &dict,
                          std::map<std::string, TSDiagnostic::Ptr> &ts_dict) {

  ShallowStressBalance::get_diagnostics_impl(dict, ts_dict);

  // replace these diagnostics
  dict["taud"] = Diagnostic::Ptr(new SSA_taud(this));
  dict["taud_mag"] = Diagnostic::Ptr(new SSA_taud_mag(this));
}

SSA_taud::SSA_taud(SSA *m)
  : Diag<SSA>(m) {

  m_dof = 2;

  // set metadata:
  m_vars.push_back(SpatialVariableMetadata(m_sys, "taud_x"));
  m_vars.push_back(SpatialVariableMetadata(m_sys, "taud_y"));

  set_attrs("X-component of the driving shear stress at the base of ice", "",
            "Pa", "Pa", 0);
  set_attrs("Y-component of the driving shear stress at the base of ice", "",
            "Pa", "Pa", 1);

  for (int k = 0; k < m_dof; ++k) {
    m_vars[k].set_string("comment",
                       "this is the driving stress used by the SSA solver");
  }
}

IceModelVec::Ptr SSA_taud::compute_impl() {

  IceModelVec2V::Ptr result(new IceModelVec2V);
  result->create(m_grid, "result", WITHOUT_GHOSTS);
  result->metadata(0) = m_vars[0];
  result->metadata(1) = m_vars[1];

  model->compute_driving_stress(*result);

  return result;
}

SSA_taud_mag::SSA_taud_mag(SSA *m)
  : Diag<SSA>(m) {

  // set metadata:
  m_vars.push_back(SpatialVariableMetadata(m_sys, "taud_mag"));

  set_attrs("magnitude of the driving shear stress at the base of ice", "",
            "Pa", "Pa", 0);
  m_vars[0].set_string("comment",
                     "this is the magnitude of the driving stress used by the SSA solver");
}

IceModelVec::Ptr SSA_taud_mag::compute_impl() {

  // Allocate memory:
  IceModelVec2S::Ptr result(new IceModelVec2S);
  result->create(m_grid, "taud_mag", WITHOUT_GHOSTS);
  result->metadata() = m_vars[0];
  result->write_in_glaciological_units = true;

  IceModelVec2V::Ptr taud = IceModelVec2V::ToVector(SSA_taud(model).compute());

  result->set_to_magnitude(*taud);

  return result;
}

//! Evaluate the ocean pressure difference term in the calving-front BC.
double SSA::ocean_pressure_difference(bool shelf, bool dry_mode, double H, double bed,
                                      double sea_level, double rho_ice, double rho_ocean,
                                      double g) {
  if (shelf) {
    // floating shelf
    return 0.5 * rho_ice * g * (1.0 - (rho_ice / rho_ocean)) * H * H;
  } else {
    // grounded terminus
    if (bed >= sea_level or dry_mode) {
      return 0.5 * rho_ice * g * H * H;
    } else {
      return 0.5 * rho_ice * g * (H * H - (rho_ocean / rho_ice) * pow(sea_level - bed, 2.0));
    }
  }
}

} // end of namespace stressbalance
} // end of namespace pism
