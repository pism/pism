// Copyright (C) 2004--2014 Constantine Khroulev, Ed Bueler, Jed Brown, Torsten Albrecht
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
#include "Mask.hh"
#include "basal_resistance.hh"
#include "PISMVars.hh"
#include "pism_options.hh"
#include "flowlaw_factory.hh"
#include "PIO.hh"
#include "enthalpyConverter.hh"
#include "error_handling.hh"

namespace pism {

SSA::SSA(IceGrid &g, EnthalpyConverter &e, const Config &c)
  : ShallowStressBalance(g, e, c)
{
  mask = NULL;
  thickness = NULL;
  tauc = NULL;
  surface = NULL;
  bed = NULL;
  enthalpy = NULL;
  driving_stress_x = NULL;
  driving_stress_y = NULL;
  gl_mask = NULL;

  strength_extension = new SSAStrengthExtension(config);
  allocate();
}

SSA::~SSA() { 
  if (deallocate() != 0) {
    PetscPrintf(grid.com, "FATAL ERROR: SSA de-allocation failed.\n");
    abort();
  }
  delete strength_extension;
}


//! \brief Initialize a generic regular-grid SSA solver.
PetscErrorCode SSA::init(Vars &vars) {
  PetscErrorCode ierr;

  ierr = ShallowStressBalance::init(vars); CHKERRQ(ierr);

  ierr = verbPrintf(2,grid.com,"* Initializing the SSA stress balance...\n"); CHKERRQ(ierr);
  ierr = verbPrintf(2, grid.com,
                    "  [using the %s flow law]\n", flow_law->name().c_str()); CHKERRQ(ierr);
  
  if (config.get_flag("sub_groundingline")) {
    gl_mask = vars.get_2d_scalar("gl_mask");
  }

  mask      = vars.get_2d_mask("mask");
  thickness = vars.get_2d_scalar("land_ice_thickness");
  tauc      = vars.get_2d_scalar("tauc");

  try {
    surface = vars.get_2d_scalar("surface_altitude");
  } catch (RuntimeError) {
    driving_stress_x = vars.get_2d_scalar("ssa_driving_stress_x");
    driving_stress_y = vars.get_2d_scalar("ssa_driving_stress_y");
  }

  bed      = vars.get_2d_scalar("bedrock_altitude");
  enthalpy = vars.get_3d_scalar("enthalpy");
  
  // Check if PISM is being initialized from an output file from a previous run
  // and read the initial guess (unless asked not to).
  bool i_set;
  std::string filename;
  ierr = OptionsString("-i", "PISM input file",
                           filename, i_set); CHKERRQ(ierr);

  if (i_set) {
    bool dont_read_initial_guess, u_ssa_found, v_ssa_found;
    unsigned int start;
    PIO nc(grid, "guess_mode");

    ierr = OptionsIsSet("-dontreadSSAvels", dont_read_initial_guess); CHKERRQ(ierr);

    nc.open(filename, PISM_READONLY);
    u_ssa_found = nc.inq_var("u_ssa");
    v_ssa_found = nc.inq_var("v_ssa");
    start = nc.inq_nrecords() - 1;
    nc.close();

    if (u_ssa_found && v_ssa_found &&
        (! dont_read_initial_guess)) {
      ierr = verbPrintf(3,grid.com,"Reading u_ssa and v_ssa...\n"); CHKERRQ(ierr);

      ierr = m_velocity.read(filename, start); CHKERRQ(ierr);
    }

  } else {
    ierr = m_velocity.set(0.0); CHKERRQ(ierr); // default initial guess
  }

  if (config.get_flag("ssa_dirichlet_bc")) {
    bc_locations = vars.get_2d_mask("bcflag");
    m_vel_bc = vars.get_2d_vector("vel_ssa_bc");
  }

  return 0;
}

//! \brief Allocate objects which any SSA solver would use.
PetscErrorCode SSA::allocate() {
  PetscErrorCode ierr;

  ierr = taud.create(grid, "taud", WITHOUT_GHOSTS); CHKERRQ(ierr);
  ierr = taud.set_attrs("diagnostic",
                        "X-component of the driving shear stress at the base of ice",
                        "Pa", "", 0); CHKERRQ(ierr);
  ierr = taud.set_attrs("diagnostic",
                        "Y-component of the driving shear stress at the base of ice",
                        "Pa", "", 1); CHKERRQ(ierr);


  // override velocity metadata
  std::vector<std::string> long_names;
  long_names.push_back("SSA model ice velocity in the X direction");
  long_names.push_back("SSA model ice velocity in the Y direction");
  ierr = m_velocity.rename("_ssa",long_names,""); CHKERRQ(ierr);

  ierr = m_velocity_global.create(grid, "bar", WITHOUT_GHOSTS); CHKERRQ(ierr);

  m_da = m_velocity_global.get_dm();

  {
    IceFlowLawFactory ice_factory(grid.com, "ssa_", config, &EC);
    ice_factory.removeType(ICE_GOLDSBY_KOHLSTEDT);

    ierr = ice_factory.setType(config.get_string("ssa_flow_law")); CHKERRQ(ierr);

    ierr = ice_factory.setFromOptions(); CHKERRQ(ierr);
    ierr = ice_factory.create(&flow_law); CHKERRQ(ierr);
  }

  return 0;
}


PetscErrorCode SSA::deallocate() {

  if (flow_law != NULL) {
    delete flow_law;
    flow_law = NULL;
  }

  return 0;
}


//! \brief Update the SSA solution.
PetscErrorCode SSA::update(bool fast, IceModelVec2S &melange_back_pressure) {
  PetscErrorCode ierr;

  if (fast)
    return 0;

  (void) melange_back_pressure;

  ierr = solve();
  if (ierr != 0) {
    PetscPrintf(grid.com, "PISM ERROR: SSA solver failed.\n");
    return ierr;
  }

  ierr = compute_basal_frictional_heating(m_velocity, *tauc, *mask,
                                          basal_frictional_heating); CHKERRQ(ierr);
  
  return 0;
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
PetscErrorCode SSA::compute_driving_stress(IceModelVec2V &result) {
  IceModelVec2S &thk = *thickness; // to improve readability (below)

  const double n = flow_law->exponent(), // frequently n = 3
    etapow  = (2.0 * n + 2.0)/n,  // = 8/3 if n = 3
    invpow  = 1.0 / etapow,  // = 3/8
    dinvpow = (- n - 2.0) / (2.0 * n + 2.0); // = -5/8
  const double minThickEtaTransform = 5.0; // m
  const double dx=grid.dx, dy=grid.dy;

  bool cfbc = config.get_flag("calving_front_stress_boundary_condition");
  bool compute_surf_grad_inward_ssa = config.get_flag("compute_surf_grad_inward_ssa");
  bool use_eta = (config.get_string("surface_gradient_method") == "eta");

  MaskQuery m(*mask);

  IceModelVec::AccessList list;
  list.add(*surface);
  list.add(*bed);
  list.add(*mask);
  list.add(thk);
  list.add(result);

  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    const double pressure = EC.getPressureFromDepth(thk(i,j)); // FIXME issue #15
    if (pressure <= 0.0) {
      result(i,j).u = 0.0;
      result(i,j).v = 0.0;
    } else {
      double h_x = 0.0, h_y = 0.0;
      // FIXME: we need to handle grid periodicity correctly.
      if (m.grounded(i,j) && (use_eta == true)) {
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
        h_x += bed->diff_x(i,j);
        h_y += bed->diff_y(i,j);
      } else {  // floating or eta transformation is not used
        if (compute_surf_grad_inward_ssa) {
          // Special case for verification tests.
          h_x = surface->diff_x_p(i,j);
          h_y = surface->diff_y_p(i,j);
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
            if ((m.grounded(i,j) && m.floating_ice(i+1,j)) || (m.floating_ice(i,j) && m.grounded(i+1,j)) ||
                (m.floating_ice(i,j) && m.ice_free_ocean(i+1,j)))
              east = 0;
            if ((m.grounded(i,j) && m.floating_ice(i-1,j)) || (m.floating_ice(i,j) && m.grounded(i-1,j)) ||
                (m.floating_ice(i,j) && m.ice_free_ocean(i-1,j)))
              west = 0;

            // This driving stress computation has to match the calving front
            // stress boundary condition in SSAFD::assemble_rhs().
            if (cfbc) {
              if (m.icy(i,j) && m.ice_free(i+1,j))
                east = 0;
              if (m.icy(i,j) && m.ice_free(i-1,j))
                west = 0;
            }

            if (east + west > 0)
              h_x = 1.0 / (west + east) * (west * surface->diff_x_stagE(i-1,j) +
                                           east * surface->diff_x_stagE(i,j));
            else
              h_x = 0.0;
          }

          // y-derivative
          {
            double south = 1, north = 1;
            if ((m.grounded(i,j) && m.floating_ice(i,j+1)) || (m.floating_ice(i,j) && m.grounded(i,j+1)) ||
                (m.floating_ice(i,j) && m.ice_free_ocean(i,j+1)))
              north = 0;
            if ((m.grounded(i,j) && m.floating_ice(i,j-1)) || (m.floating_ice(i,j) && m.grounded(i,j-1)) ||
                (m.floating_ice(i,j) && m.ice_free_ocean(i,j-1)))
              south = 0;

            // This driving stress computation has to match the calving front
            // stress boundary condition in SSAFD::assemble_rhs().
            if (cfbc) {
              if (m.icy(i,j) && m.ice_free(i,j+1))
                north = 0;
              if (m.icy(i,j) && m.ice_free(i,j-1))
                south = 0;
            }

            if (north + south > 0)
              h_y = 1.0 / (south + north) * (south * surface->diff_y_stagN(i,j-1) +
                                             north * surface->diff_y_stagN(i,j));
            else
              h_y = 0.0;
          }

        } // end of "general case"

      } // end of "floating or eta transformation is not used"

      result(i,j).u = - pressure * h_x;
      result(i,j).v = - pressure * h_y;
    } // end of "(pressure > 0)"
  }

  return 0;
}

PetscErrorCode SSA::stdout_report(std::string &result) {
  result = stdout_ssa;
  return 0;
}


//! \brief Set the initial guess of the SSA velocity.
PetscErrorCode SSA::set_initial_guess(IceModelVec2V &guess) {
  PetscErrorCode ierr;
  ierr = m_velocity.copy_from(guess); CHKERRQ(ierr);
  return 0;
}


void SSA::add_vars_to_output(const std::string &/*keyword*/, std::set<std::string> &result) {
  result.insert("vel_ssa");
}


PetscErrorCode SSA::define_variables(const std::set<std::string> &vars, const PIO &nc, IO_Type nctype) {
  PetscErrorCode ierr;

  if (set_contains(vars, "vel_ssa")) {
    ierr = m_velocity.define(nc, nctype); CHKERRQ(ierr);
  }

  return 0;
}


PetscErrorCode SSA::write_variables(const std::set<std::string> &vars, const PIO &nc) {
  PetscErrorCode ierr;

  if (set_contains(vars, "vel_ssa")) {
    ierr = m_velocity.write(nc); CHKERRQ(ierr);
  }

  return 0;
}

void SSA::get_diagnostics(std::map<std::string, Diagnostic*> &dict,
                          std::map<std::string, TSDiagnostic*> &ts_dict) {

  ShallowStressBalance::get_diagnostics(dict, ts_dict);

  if (dict["taud"] != NULL) {
    delete dict["taud"];
  }
  dict["taud"] = new SSA_taud(this, grid, *variables);

  if (dict["taud_mag"] != NULL) {
    delete dict["taud_mag"];
  }
  dict["taud_mag"] = new SSA_taud_mag(this, grid, *variables);
}

SSA_taud::SSA_taud(SSA *m, IceGrid &g, Vars &my_vars)
  : Diag<SSA>(m, g, my_vars) {

  dof = 2;
  vars.resize(dof,  NCSpatialVariable(g.get_unit_system()));
  // set metadata:
  vars[0].init_2d("taud_x", grid);
  vars[1].init_2d("taud_y", grid);

  set_attrs("X-component of the driving shear stress at the base of ice", "",
            "Pa", "Pa", 0);
  set_attrs("Y-component of the driving shear stress at the base of ice", "",
            "Pa", "Pa", 1);

  for (int k = 0; k < dof; ++k)
    vars[k].set_string("comment",
                       "this is the driving stress used by the SSA solver");
}

PetscErrorCode SSA_taud::compute(IceModelVec* &output) {
  PetscErrorCode ierr;

  IceModelVec2V *result = new IceModelVec2V;
  ierr = result->create(grid, "result", WITHOUT_GHOSTS); CHKERRQ(ierr);
  result->metadata() = vars[0];
  result->metadata(1) = vars[1];

  ierr = model->compute_driving_stress(*result); CHKERRQ(ierr);

  output = result;
  return 0;
}

SSA_taud_mag::SSA_taud_mag(SSA *m, IceGrid &g, Vars &my_vars)
  : Diag<SSA>(m, g, my_vars) {

  // set metadata:
  vars[0].init_2d("taud_mag", grid);

  set_attrs("magnitude of the driving shear stress at the base of ice", "",
            "Pa", "Pa", 0);
  vars[0].set_string("comment",
                     "this is the magnitude of the driving stress used by the SSA solver");
}

PetscErrorCode SSA_taud_mag::compute(IceModelVec* &output) {
  PetscErrorCode ierr;

  // Allocate memory:
  IceModelVec2S *result = new IceModelVec2S;
  ierr = result->create(grid, "taud_mag", WITHOUT_GHOSTS); CHKERRQ(ierr);
  result->metadata() = vars[0];
  result->write_in_glaciological_units = true;

  IceModelVec* tmp;
  SSA_taud diag(model, grid, variables);

  ierr = diag.compute(tmp);

  IceModelVec2V *taud = dynamic_cast<IceModelVec2V*>(tmp);
  if (taud == NULL)
    SETERRQ(grid.com, 1, "expected an IceModelVec2V, but dynamic_cast failed");

  ierr = taud->magnitude(*result); CHKERRQ(ierr);

  delete tmp;

  output = result;
  return 0;
}


} // end of namespace pism
