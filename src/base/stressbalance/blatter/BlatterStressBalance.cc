// Copyright (C) 2010-2016 Ed Bueler and Constantine Khroulev
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

#include "BlatterStressBalance.hh"
#include "coupler/PISMOcean.hh"
#include "base/util/IceGrid.hh"
#include "base/util/PISMVars.hh"
#include "base/basalstrength/basal_resistance.hh"
#include "FE3DTools.h"
#include "base/enthalpyConverter.hh"
#include "base/rheology/FlowLaw.hh"
#include "base/rheology/FlowLawFactory.hh"
#include "base/util/PISMVars.hh"
#include "base/util/error_handling.hh"
#include "base/util/pism_const.hh"
#include "base/util/pism_utilities.hh"

namespace pism {
namespace stressbalance {

/*!
 * FIXMEs:
 *
 * We need to allow spatially-variable (and depth-dependent) ice hardness.
 *
 * We need to compute the volumetric strain heating.
 *
 */

//! C-wrapper for PISM's IceFlowLaw::viscosity().
extern "C"
PetscErrorCode viscosity(void *ctx, double hardness, double gamma,
                         double *eta, double *deta) {
  BlatterQ1Ctx *blatter_ctx = (BlatterQ1Ctx*)ctx;
  BlatterStressBalance *blatter_stress_balance = (BlatterStressBalance*)blatter_ctx->extra;

  try {
    blatter_stress_balance->flow_law()->effective_viscosity(hardness, gamma, eta, deta);
  } catch (...) {
    MPI_Comm com = blatter_stress_balance->grid()->ctx()->com();
    handle_fatal_errors(com);
    SETERRQ(com, 1, "A PISM callback failed");
  }
  return 0;
}

//! C-wrapper for PISM's IceBasalResistancePlasticLaw::dragWithDerivative().
extern "C"
PetscErrorCode drag(void *ctx, double tauc, double u, double v,
                    double *taud, double *dtaub) {
  BlatterQ1Ctx *blatter_ctx = (BlatterQ1Ctx*)ctx;
  BlatterStressBalance *blatter_stress_balance = (BlatterStressBalance*)blatter_ctx->extra;

  try {
    blatter_stress_balance->sliding_law()->drag_with_derivative(tauc, u, v, taud, dtaub);
  } catch (...) {
    MPI_Comm com = blatter_stress_balance->grid()->ctx()->com();
    handle_fatal_errors(com);
    SETERRQ(com, 1, "A PISM callback failed");
  }
  return 0;
}

/*! @brief A no-cost wrapper around 2D arrays. Hides the indexing order. */
template <typename T>
class PointerWrapper2D {
  T **m_data;
public:
  PointerWrapper2D() {
    m_data = NULL;
  }

  T*** address() {
    return &m_data;
  }

  T& operator()(int i, int j) {
    // NOTE: this indexing order is important!
    return m_data[i][j];
  }
};

/*! @brief A no-cost wrapper around 3D arrays. Hides the indexing order. */
template <typename T>
class PointerWrapper3D {
  T ***m_data;
public:
  PointerWrapper3D() {
    m_data = NULL;
  }

  T**** address() {
    return &m_data;
  }

  T& operator()(int i, int j, int k) {
    // NOTE: this indexing order is important!
    return m_data[i][j][k];
  }
};

BlatterStressBalance::BlatterStressBalance(IceGrid::ConstPtr g,
                                           EnthalpyConverter::Ptr e)
  : ShallowStressBalance(g, e), m_min_thickness(10.0) {

  Config::ConstPtr config = g->ctx()->config();

  int blatter_Mz = (int)config->get_double("blatter_Mz");
  m_da2 = g->get_dm(1, (int)config->get_double("grid_max_stencil_width"));

  m_ctx.Lx = 2.0 * g->Lx();
  m_ctx.Ly = 2.0 * g->Ly();
  m_ctx.dirichlet_scale = 1.0;
  m_ctx.rhog = config->get_double("ice_density") * config->get_double("standard_gravity");
  m_ctx.no_slip = PETSC_TRUE;	// FIXME (at least make configurable)
  m_ctx.nonlinear.viscosity = viscosity;
  m_ctx.nonlinear.drag = drag;
  m_ctx.extra = this;
  initialize_Q12D(m_ctx.Q12D.chi, m_ctx.Q12D.dchi);
  initialize_Q13D(m_ctx.Q13D.chi, m_ctx.Q13D.dchi);

  PetscErrorCode ierr = BlatterQ1_create(g->com, *m_da2, blatter_Mz,
                                         &this->m_ctx, this->m_snes.rawptr());
  PISM_CHK(ierr, "BlatterQ1_create");

  // now allocate u, v, and strain_heating (strain heating)
  m_u.create(grid(), "uvel", WITH_GHOSTS);
  m_u.set_attrs("diagnostic", "horizontal velocity of ice in the X direction",
              "m s-1", "land_ice_x_velocity");
  m_u.metadata().set_string("glaciological_units", "m year-1");
  m_u.write_in_glaciological_units = true;

  m_v.create(grid(), "vvel", WITH_GHOSTS);
  m_v.set_attrs("diagnostic", "horizontal velocity of ice in the Y direction",
              "m s-1", "land_ice_y_velocity");
  m_v.metadata().set_string("glaciological_units", "m year-1");
  m_v.write_in_glaciological_units = true;

  // strain_heating
  m_strain_heating.create(grid(), "strainheat", WITHOUT_GHOSTS); // never diff'ed in hor dirs
  m_strain_heating.set_attrs("internal",
                             "rate of strain heating in ice (dissipation heating)",
                             "W m-3", "");
  m_strain_heating.metadata().set_string("glaciological_units", "mW m-3");

  std::vector<double> sigma(blatter_Mz);
  double dz = 1.0 / (blatter_Mz - 1);
  for (int i = 0; i < blatter_Mz; ++i) {
    sigma[i] = i * dz;
  }
  sigma.back() = 1.0;

  std::map<std::string,std::string> z_attrs;
  z_attrs["axis"]          = "Z";
  z_attrs["long_name"]     = "scaled Z-coordinate in the ice (z_base=0, z_surface=1)";
  z_attrs["units"]         = "1";
  z_attrs["positive"]      = "up";

  // storage for u and v on the sigma vertical grid (for restarting)
  m_u_sigma.create(grid(), "uvel_sigma", "z_sigma", sigma, z_attrs);
  m_u_sigma.set_attrs("diagnostic",
                    "horizontal velocity of ice in the X direction on the sigma vertical grid",
                    "m s-1", "");
  m_u_sigma.metadata().set_string("glaciological_units", "m year-1");
  m_u_sigma.write_in_glaciological_units = false;

  m_v_sigma.create(grid(), "vvel_sigma", "z_sigma", sigma, z_attrs);
  m_v_sigma.set_attrs("diagnostic",
                    "horizontal velocity of ice in the Y direction on the sigma vertical grid",
                    "m s-1", "");
  m_v_sigma.metadata().set_string("glaciological_units", "m year-1");
  m_v_sigma.write_in_glaciological_units = false;

  {
    rheology::FlowLawFactory ice_factory("blatter_", config, e);
    ice_factory.remove(ICE_GOLDSBY_KOHLSTEDT);

    ice_factory.set_default(config->get_string("blatter_flow_law"));

    m_flow_law = ice_factory.create();
  }
}

BlatterStressBalance::~BlatterStressBalance() {
}

void BlatterStressBalance::init_impl() {

  const Vars &vars = m_grid->variables();

  m_bed_elevation = vars.get_2d_scalar("bedrock_altitude");

  m_ice_thickness = vars.get_2d_scalar("land_ice_thickness");

  m_tauc = vars.get_2d_scalar("tauc");

  m_enthalpy = vars.get_3d_scalar("enthalpy");
}

void BlatterStressBalance::update(bool fast, const IceModelVec2S &melange_back_pressure) {

  (void) melange_back_pressure;

  if (fast) {
    throw RuntimeError("'fast' mode not meaningful for BlatterStressBalance");
  }

  // setup
  setup();

  // solve
  PetscErrorCode ierr = SNESSolve(this->m_snes, PETSC_NULL, PETSC_NULL);
  PISM_CHK(ierr, "SNESSolve");

  // Transfer solution from the FEM mesh to the regular grid used in the rest
  // of PISM and compute the vertically-averaged velocity.
  transfer_velocity();

  // Copy solution from the SNES to m_u_sigma and m_v_sigma.
  copy_velocity(FROM_SNES_STORAGE);

  compute_volumetric_strain_heating();
}

/*! \brief Set up model parameters on the fine grid. */
/*!
 * This method expects bed_elevation, ice_thickness, and tauc to have width=1 ghosts.
 *
 * We should also compute ice hardness on the "sigma" grid here.
 */
void BlatterStressBalance::setup() {
  PetscErrorCode ierr;
  Vec param_vec;
  PointerWrapper2D<PrmNode> parameters;
  DM da;
  double
    ice_density = m_config->get_double("ice_density"),
    sea_water_density = m_config->get_double("sea_water_density"),
    alpha = ice_density / sea_water_density;

  ierr = SNESGetDM(this->m_snes, &da); PISM_CHK(ierr, "SNESGetDM");

  ierr = BlatterQ1_begin_2D_parameter_access(da, PETSC_FALSE, &param_vec, parameters.address());
  PISM_CHK(ierr, "BlatterQ1_begin_2D_parameter_access");

  const IceModelVec2S
    &bed       = *m_bed_elevation,
    &thickness = *m_ice_thickness,
    &tauc      = *m_tauc;

  IceModelVec::AccessList list;
  list.add(bed);
  list.add(thickness);
  list.add(tauc);

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    // compute the elevation of the bottom surface of the ice
    if (bed(i,j) > -alpha * thickness(i,j)) {
      // grounded
      parameters(i, j).ice_bottom = bed(i,j);
    } else {
      // floating
      parameters(i, j).ice_bottom = -alpha * thickness(i,j);
    }

    parameters(i, j).thickness = thickness(i,j);

    // fudge ice thickness (FIXME!!!)
    if (thickness(i,j) < m_min_thickness)
      parameters(i, j).thickness += m_min_thickness;

    parameters(i, j).tauc = tauc(i,j);
  }

  ierr = BlatterQ1_end_2D_parameter_access(da, PETSC_FALSE, &param_vec, parameters.address());
  PISM_CHK(ierr, "BlatterQ1_end_2D_parameter_access");
}

//! Initialize ice hardness on the "sigma" grid.
void BlatterStressBalance::initialize_ice_hardness() {
  PetscErrorCode ierr;

  PointerWrapper3D<PetscScalar> hardness;
  unsigned int Mz_fem = static_cast<unsigned int>(m_config->get_double("blatter_Mz"));
  DM da;
  Vec hardness_vec;

  ierr = SNESGetDM(this->m_snes, &da); PISM_CHK(ierr, "SNESGetDM");

  ierr = BlatterQ1_begin_hardness_access(da, PETSC_FALSE, &hardness_vec, hardness.address());
  PISM_CHK(ierr, "BlatterQ1_begin_hardness_access");

  IceModelVec::AccessList list;
  list.add(*m_enthalpy);
  list.add(*m_ice_thickness);

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    double thk = (*m_ice_thickness)(i,j);

    // fudge ice thickness (FIXME!!!)
    if (thk < m_min_thickness)
      thk += m_min_thickness;

    double dz_fem = thk / (Mz_fem - 1);
    const double *E = m_enthalpy->get_column(i, j);

    // compute ice hardness on the sigma grid
    for (unsigned int k = 0; k < Mz_fem; ++k) {
      double z_fem = k * dz_fem,
        depth = thk - z_fem,
        pressure = m_EC->pressure(depth),
        E_local;
      unsigned int k0 = m_grid->kBelowHeight(z_fem);

      const std::vector<double> &zlevels = m_enthalpy->get_levels();
      const unsigned int Mz = m_grid->Mz();

      if (k0 + 1 < Mz) {
        double lambda = (z_fem - zlevels[k0]) / (zlevels[k0+1] - zlevels[k0]);

        E_local = (1.0 - lambda) * E[k0] + lambda * E[k0 + 1];
      } else {
        // should never happen
        E_local = E[Mz-1];
      }

      hardness(i, j, k) = m_flow_law->hardness(E_local, pressure);
    }
  }

  ierr = BlatterQ1_end_hardness_access(da, PETSC_FALSE, &hardness_vec, hardness.address());
  PISM_CHK(ierr, "BlatterQ1_end_hardness_access");
}

//! Transfer the velocity field from the FEM "sigma" grid onto PISM's grid.
/*!
 * We also compute vertically-averaged ice velocity here.
 */
void BlatterStressBalance::transfer_velocity() {
  PetscErrorCode ierr;

  PointerWrapper3D<Vector2> U;
  double *u_ij, *v_ij;
  DM da;
  Vec X;
  unsigned int Mz_fem = static_cast<unsigned int>(m_config->get_double("blatter_Mz"));

  ierr = SNESGetDM(this->m_snes, &da); PISM_CHK(ierr, "SNESGetDM");
  ierr = SNESGetSolution(this->m_snes, &X); PISM_CHK(ierr, "SNESGetSolution");

  ierr = DMDAVecGetArray(da, X, U.address()); PISM_CHK(ierr, "DMDAVecGetArray");

  IceModelVec::AccessList list;
  list.add(m_u);
  list.add(m_v);
  list.add(*m_ice_thickness);
  list.add(m_velocity);

  const unsigned int Mz = m_grid->Mz();
  const std::vector<double> &zlevels = m_u.get_levels();

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    u_ij = m_u.get_column(i, j);
    v_ij = m_v.get_column(i, j);

    double thk = (*m_ice_thickness)(i,j);

    // fudge ice thickness (FIXME!!!)
    if (thk < m_min_thickness)
      thk += m_min_thickness;

    double dz_fem = thk / (Mz_fem - 1);

    // compute vertically-averaged velocity using trapezoid rule
    double ubar = 0, vbar = 0;
    for (unsigned int k = 0; k < Mz_fem - 1; ++k) {
      ubar += U(i, j, k).u + U(i, j, k+1).u;
      vbar += U(i, j, k).v + U(i, j, k+1).v;
    }
    // finish the traperoidal rule (1/2 * dz) and compute the average:
    m_velocity(i,j).u = ubar * (0.5*dz_fem) / thk;
    m_velocity(i,j).v = vbar * (0.5*dz_fem) / thk;

    // compute 3D horizontal velocity
    unsigned int current_level = 0;
    for (unsigned int k = 0; k < Mz; ++k) {

      // find the FEM grid level just below the current PISM grid level
      while ((current_level + 1) * dz_fem < zlevels[k])
        current_level++;

      if (current_level + 1 < Mz_fem) {
        // use linear interpolation
        double z0 = current_level * dz_fem,
          lambda = (zlevels[k] - z0) / dz_fem;

        u_ij[k] = (U(i, j, current_level).u * (1 - lambda) +
                   U(i, j, current_level+1).u * lambda);

        v_ij[k] = (U(i, j, current_level).v * (1 - lambda) +
                   U(i, j, current_level+1).v * lambda);

      } else {
        // extrapolate above the surface
        u_ij[k] = U(i, j, Mz_fem-1).u;
        v_ij[k] = U(i, j, Mz_fem-1).v;
      }

    }	// k-loop
  } // loop over map-plane grid points

  ierr = DMDAVecRestoreArray(da, X, U.address()); PISM_CHK(ierr, "DMDAVecRestoreArray");

  m_u.update_ghosts();
  m_v.update_ghosts();
}

//! Copy velocity from a dof=2 vector to special storage (to save it for re-starting).
void BlatterStressBalance::copy_velocity(Direction direction) {
  PetscErrorCode ierr;

  PointerWrapper3D<Vector2> U;
  DM da;
  Vec X;
  unsigned int Mz_fem = static_cast<unsigned int>(m_config->get_double("blatter_Mz"));

  ierr = SNESGetDM(this->m_snes, &da); PISM_CHK(ierr, "SNESGetDM");
  ierr = SNESGetSolution(this->m_snes, &X); PISM_CHK(ierr, "SNESGetSolution");

  ierr = DMDAVecGetArray(da, X, U.address()); PISM_CHK(ierr, "DMDAVecGetArray");

  IceModelVec::AccessList list;
  list.add(m_u_sigma);
  list.add(m_v_sigma);

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    double
      *u_ij = m_u_sigma.get_column(i, j),
      *v_ij = m_v_sigma.get_column(i, j);

    if (direction == FROM_SNES_STORAGE) {
      for (unsigned int k = 0; k < Mz_fem; ++k) {
        u_ij[k] = U(i, j, k).u;
        v_ij[k] = U(i, j, k).v;
      }
    } else {
      for (unsigned int k = 0; k < Mz_fem; ++k) {
        U(i, j, k).u = u_ij[k];
        U(i, j, k).v = v_ij[k];
      }
    }
  } // loop over map-plane grid points

  ierr = DMDAVecRestoreArray(da, X, U.address()); PISM_CHK(ierr, "DMDAVecRestoreArray");
}

const IceModelVec3& BlatterStressBalance::velocity_u() const {
  return m_u;
}

const IceModelVec3& BlatterStressBalance::velocity_v() const {
  return m_v;
}

const IceModelVec3Custom& BlatterStressBalance::velocity_u_sigma() const {
  return m_u_sigma;
}

const IceModelVec3Custom& BlatterStressBalance::velocity_v_sigma() const {
  return m_v_sigma;
}

void BlatterStressBalance::compute_volumetric_strain_heating() {
  // FIXME: implement this
  m_strain_heating.set(0.0);
}

void BlatterStressBalance::add_vars_to_output_impl(const std::string &keyword,
                                                   std::set<std::string> &result) {
  (void) keyword;

  result.insert("u_sigma");
  result.insert("v_sigma");
}

//! Defines requested couplings fields.
void BlatterStressBalance::define_variables_impl(const std::set<std::string> &vars,
                                                 const PIO &nc,
                                                 IO_Type nctype) {
  if (set_contains(vars, "u_sigma")) {
    m_u_sigma.define(nc, nctype);
  }

  if (set_contains(vars, "v_sigma")) {
    m_v_sigma.define(nc, nctype);
  }
}

//! Writes requested couplings fields to file.
void BlatterStressBalance::write_variables_impl(const std::set<std::string> &vars, const PIO &nc) {

  if (set_contains(vars, "u_sigma") or set_contains(vars, "v_sigma")) {
    copy_velocity(FROM_SNES_STORAGE);
  }

  if (set_contains(vars, "u_sigma")) {
    m_u_sigma.write(nc);
  }

  if (set_contains(vars, "v_sigma")) {
    m_v_sigma.write(nc);
  }
}

//! \brief Produce a report string for the standard output.
std::string BlatterStressBalance::stdout_report() {
  return "";
}


} // end of namespace stressbalance
} // end of namespace pism
