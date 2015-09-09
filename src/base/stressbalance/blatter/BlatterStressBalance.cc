// Copyright (C) 2010-2015 Ed Bueler and Constantine Khroulev
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
#include "base/rheology/flowlaws.hh"
#include "base/rheology/flowlaw_factory.hh"
#include "base/util/PISMVars.hh"
#include "base/util/error_handling.hh"

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
void viscosity(void *ctx, double hardness, double gamma,
	       double *eta, double *deta) {
  BlatterQ1Ctx *blatter_ctx = (BlatterQ1Ctx*)ctx;
  BlatterStressBalance *blatter_stress_balance = (BlatterStressBalance*)blatter_ctx->extra;

  blatter_stress_balance->m_flow_law->effective_viscosity(hardness, gamma, eta, deta);
}

//! C-wrapper for PISM's IceBasalResistancePlasticLaw::dragWithDerivative().
void drag(void *ctx, double tauc, double u, double v,
	  double *taud, double *dtaub) {
  BlatterQ1Ctx *blatter_ctx = (BlatterQ1Ctx*)ctx;
  BlatterStressBalance *blatter_stress_balance = (BlatterStressBalance*)blatter_ctx->extra;

  blatter_stress_balance->basal_sliding_law->drag_with_derivative(tauc, u, v, taud, dtaub);
}

BlatterStressBalance::BlatterStressBalance(IceGrid::ConstPtr g,
					   EnthalpyConverter::Ptr e)
  : ShallowStressBalance(g, e), min_thickness(10.0)
{

  Config::ConstPtr config = g->ctx()->config();

  int blatter_Mz = (int)config->get_double("blatter_Mz");
  m_da2 = g->get_dm(1, (int)config->get_double("grid_max_stencil_width"));

  ctx.Lx = 2.0 * g->Lx();
  ctx.Ly = 2.0 * g->Ly();
  ctx.dirichlet_scale = 1.0;
  ctx.rhog = config->get_double("ice_density") * config->get_double("standard_gravity");
  ctx.no_slip = PETSC_TRUE;	// FIXME (at least make configurable)
  ctx.nonlinear.viscosity = viscosity;
  ctx.nonlinear.drag = drag;
  ctx.extra = this;
  initialize_Q12D(ctx.Q12D.chi, ctx.Q12D.dchi);
  initialize_Q13D(ctx.Q13D.chi, ctx.Q13D.dchi);

  PetscErrorCode ierr = BlatterQ1_create(g->com, m_da2->get(), blatter_Mz,
                                         &this->ctx, &this->snes);
  PISM_CHK(ierr, "BlatterQ1_create");

  // now allocate u, v, and strain_heating (strain heating)
  u.create(grid(), "uvel", WITH_GHOSTS);
  u.set_attrs("diagnostic", "horizontal velocity of ice in the X direction",
              "m s-1", "land_ice_x_velocity");
  u.metadata().set_string("glaciological_units", "m year-1");
  u.write_in_glaciological_units = true;

  v.create(grid(), "vvel", WITH_GHOSTS);
  v.set_attrs("diagnostic", "horizontal velocity of ice in the Y direction",
              "m s-1", "land_ice_y_velocity");
  v.metadata().set_string("glaciological_units", "m year-1");
  v.write_in_glaciological_units = true;

  // strain_heating
  strain_heating.create(grid(), "strainheat", WITHOUT_GHOSTS); // never diff'ed in hor dirs
  strain_heating.set_attrs("internal",
                           "rate of strain heating in ice (dissipation heating)",
                           "W m-3", "");
  strain_heating.metadata().set_string("glaciological_units", "mW m-3");

  std::vector<double> sigma(blatter_Mz);
  double dz = 1.0 / (blatter_Mz - 1);
  for (int i = 0; i < blatter_Mz; ++i)
    sigma[i] = i * dz;
  sigma.back() = 1.0;

  std::map<std::string,std::string> z_attrs;
  z_attrs["axis"]          = "Z";
  z_attrs["long_name"]     = "scaled Z-coordinate in the ice (z_base=0, z_surface=1)";
  z_attrs["units"]         = "1";
  z_attrs["positive"]      = "up";

  // storage for u and v on the sigma vertical grid (for restarting)
  u_sigma.create(grid(), "uvel_sigma", "z_sigma", sigma, z_attrs);
  u_sigma.set_attrs("diagnostic",
                    "horizontal velocity of ice in the X direction on the sigma vertical grid",
                    "m s-1", "");
  u_sigma.metadata().set_string("glaciological_units", "m year-1");
  u_sigma.write_in_glaciological_units = false;

  v_sigma.create(grid(), "vvel_sigma", "z_sigma", sigma, z_attrs);
  v_sigma.set_attrs("diagnostic",
                    "horizontal velocity of ice in the Y direction on the sigma vertical grid",
                    "m s-1", "");
  v_sigma.metadata().set_string("glaciological_units", "m year-1");
  v_sigma.write_in_glaciological_units = false;

  {
    rheology::FlowLawFactory ice_factory("blatter_", config, e);
    ice_factory.remove_type(ICE_GOLDSBY_KOHLSTEDT);

    ice_factory.set_default_type(config->get_string("blatter_flow_law"));

    m_flow_law = ice_factory.create();
  }
}

BlatterStressBalance::~BlatterStressBalance() {
  PetscErrorCode ierr = SNESDestroy(&this->snes); PISM_CHK(ierr, "SNESDestroy");
}

PetscErrorCode BlatterStressBalance::init() {

  const Vars &vars = m_grid->variables();

  bed_elevation = vars.get_2d_scalar("bedrock_altitude");

  ice_thickness = vars.get_2d_scalar("land_ice_thickness");

  tauc = vars.get_2d_scalar("tauc");

  enthalpy = vars.get_3d_scalar("enthalpy");

  return 0;
}

void BlatterStressBalance::update(bool fast, const IceModelVec2S &melange_back_pressure) {

  (void) melange_back_pressure;

  if (fast) {
    throw RuntimeError("'fast' mode not meaningful for BlatterStressBalance");
  }

  // setup
  setup();

  // solve
  PetscErrorCode ierr = SNESSolve(this->snes, PETSC_NULL, PETSC_NULL);
  PISM_CHK(ierr, "SNESSolve");

  // Transfer solution from the FEM mesh to the regular grid used in the rest
  // of PISM and compute the vertically-averaged velocity.
  transfer_velocity();

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
  PrmNode **parameters;
  DM da;
  double
    ice_density = m_config.get_double("ice_density"),
    sea_water_density = m_config.get("sea_water_density"),
    alpha = ice_density / sea_water_density;

  ierr = SNESGetDM(this->snes, &da); CHKERRQ(ierr);

  ierr = BlatterQ1_begin_2D_parameter_access(da, &parameters); CHKERRQ(ierr);

  IceModelVec::AccessList list;
  list.add(*bed_elevation);
  list.add(*ice_thickness);
  list.add(*tauc);

  int GHOSTS = 1;
  for (PointsWithGhosts p(*m_grid, GHOSTS); p; p.next()) {
    const int i = p.i(), j = p.j();

    // compute the elevation of the bottom surface of the ice
    if ((*bed_elevation)(i,j) > -alpha * (*ice_thickness)(i,j)) {
      // grounded
      parameters[i][j].ice_bottom = (*bed_elevation)(i,j);
    } else {
      // floating
      parameters[i][j].ice_bottom = -alpha * (*ice_thickness)(i,j);
    }

    parameters[i][j].thickness = (*ice_thickness)(i,j);

    // fudge ice thickness (FIXME!!!)
    if ((*ice_thickness)(i,j) < min_thickness)
      parameters[i][j].thickness += min_thickness;

    parameters[i][j].tauc = (*tauc)(i,j);
  }

  ierr = BlatterQ1_end_2D_parameter_access(da, &parameters); CHKERRQ(ierr);

  return 0;
}

//! Initialize ice hardness on the "sigma" grid.
PetscErrorCode BlatterStressBalance::initialize_ice_hardness() {
  PetscErrorCode ierr;
  PetscScalar ***hardness, *E;
  unsigned int Mz_fem = static_cast<unsigned int>(config.get("blatter_Mz"));
  DM da;


  ierr = SNESGetDM(this->snes, &da); CHKERRQ(ierr);

  ierr = BlatterQ1_begin_hardness_access(da, &hardness); CHKERRQ(ierr);
  ierr = enthalpy->begin_access(); CHKERRQ(ierr);
  ierr = ice_thickness->begin_access(); CHKERRQ(ierr);

  int GHOSTS = 1;
  for (int   i = grid.xs - GHOSTS; i < grid.xs+grid.xm + GHOSTS; ++i) {
    for (int j = grid.ys - GHOSTS; j < grid.ys+grid.ym + GHOSTS; ++j) {
      double thk = (*ice_thickness)(i,j);

      // fudge ice thickness (FIXME!!!)
      if (thk < min_thickness)
	thk += min_thickness;

      double dz_fem = thk / (Mz_fem - 1);
      ierr = enthalpy->getInternalColumn(i, j, &E); CHKERRQ(ierr);

      // compute ice hardness on the sigma grid
      for (unsigned int k = 0; k < Mz_fem; ++k) {
	double z_fem = k * dz_fem,
	  depth = thk - z_fem,
	  pressure = EC.getPressureFromDepth(depth),
	  E_local;
	unsigned int k0 = grid.kBelowHeight(z_fem);
	
	if (k0 + 1 < grid.Mz) {
	  double lambda = (z_fem - grid.zlevels[k0]) / (grid.zlevels[k0+1] - grid.zlevels[k0]);
	  
	  E_local = (1.0 - lambda) * E[k0] + lambda * E[k0 + 1];
	} else {
	  // should never happen
	  E_local = E[grid.Mz-1];
	}

	hardness[i][j][k] = flow_law->hardness_parameter(E_local, pressure);
      }
    }
  }

  ierr = ice_thickness->end_access(); CHKERRQ(ierr);
  ierr = enthalpy->end_access(); CHKERRQ(ierr);
  ierr = BlatterQ1_end_hardness_access(da, &hardness); CHKERRQ(ierr);

  return 0;
}

//! Transfer the velocity field from the FEM "sigma" grid onto PISM's grid.
/*!
 * We also compute vertically-averaged ice velocity here.
 */
void BlatterStressBalance::transfer_velocity() {
  PetscErrorCode ierr;

  Vector2 ***U;
  double *u_ij, *v_ij;
  DM da;
  Vec X;
  unsigned int Mz_fem = static_cast<unsigned int>(config.get("blatter_Mz"));

  ierr = SNESGetDM(this->snes, &da); PISM_CHK(ierr, "SNESGetDM");
  ierr = SNESGetSolution(this->snes, &X); PISM_CHK(ierr, "SNESGetSolution");

  ierr = DMDAVecGetArray(da, X, &U); PISM_CHK(ierr, "DMDAVecGetArray");

  IceModelVec::AccessList list;
  list.add(u);
  list.add(v);
  list.add(*ice_thickness);
  list.add(m_velocity);
 
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    u_ij = u.get_column(i, j, &u_ij);
    v_ij = v.get_column(i, j);

    double thk = (*ice_thickness)(i,j);

    // fudge ice thickness (FIXME!!!)
    if (thk < min_thickness)
      thk += min_thickness;

    double dz_fem = thk / (Mz_fem - 1);

    // compute vertically-averaged velocity using trapezoid rule
    double ubar = 0, vbar = 0;
    for (unsigned int k = 0; k < Mz_fem - 1; ++k) {
      ubar += U[i][j][k].u + U[i][j][k+1].u;
      vbar += U[i][j][k].v + U[i][j][k+1].v;
    }
    // finish the traperoidal rule (1/2 * dz) and compute the average:
    m_velocity(i,j).u = ubar * (0.5*dz_fem) / thk;
    m_velocity(i,j).v = vbar * (0.5*dz_fem) / thk;
      
    // compute 3D horizontal velocity
    unsigned int current_level = 0;
    for (unsigned int k = 0; k < grid.Mz; ++k) {

      // find the FEM grid level just below the current PISM grid level
      while ((current_level + 1) * dz_fem < grid.zlevels[k])
        current_level++;

      if (current_level + 1 < Mz_fem) {
        // use linear interpolation
        double z0 = current_level * dz_fem,
          lambda = (grid.zlevels[k] - z0) / dz_fem;

        u_ij[k] = (U[i][j][current_level].u * (1 - lambda) +
                   U[i][j][current_level+1].u * lambda);

        v_ij[k] = (U[i][j][current_level].v * (1 - lambda) +
                   U[i][j][current_level+1].v * lambda);

      } else {
        // extrapolate above the surface
        u_ij[k] = U[i][j][Mz_fem-1].u;
        v_ij[k] = U[i][j][Mz_fem-1].v;
      }

    }	// k-loop
  } // loop over map-plane grid points

  ierr = DMDAVecRestoreArray(da, X, &U); PISM_CHK(ierr, "DMDAVecRestoreArray");

  u.update_ghosts();
  v.update_ghosts();
}

//! Copy velocity from a dof=2 vector to special storage (to save it for re-starting).
void BlatterStressBalance::save_velocity() {
  PetscErrorCode ierr;

  Vector2 ***U;
  PetscScalar *u_ij, *v_ij;
  DM da;
  Vec X;
  unsigned int Mz_fem = static_cast<unsigned int>(config.get("blatter_Mz"));

  ierr = SNESGetDM(this->snes, &da); CHKERRQ(ierr);
  ierr = SNESGetSolution(this->snes, &X); CHKERRQ(ierr);

  ierr = DMDAVecGetArray(da, X, &U); CHKERRQ(ierr);

  ierr = u_sigma.begin_access(); CHKERRQ(ierr);
  ierr = v_sigma.begin_access(); CHKERRQ(ierr);
 
  for (int   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (int j = grid.ys; j < grid.ys+grid.ym; ++j) {
      ierr = u_sigma.getInternalColumn(i, j, &u_ij); CHKERRQ(ierr);
      ierr = v_sigma.getInternalColumn(i, j, &v_ij); CHKERRQ(ierr);

      for (unsigned int k = 0; k < Mz_fem; ++k) {
	u_ij[k] = U[i][j][k].u;
	v_ij[k] = U[i][j][k].v;
      }
    } // j-loop
  } // i-loop

  ierr = v_sigma.end_access(); CHKERRQ(ierr);
  ierr = u_sigma.end_access(); CHKERRQ(ierr);

  ierr = DMDAVecRestoreArray(da, X, &U); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode BlatterStressBalance::compute_volumetric_strain_heating() {
  PetscErrorCode ierr;
  // FIXME: implement this
  ierr = strain_heating.set(0.0); CHKERRQ(ierr);
  
  return 0;
}

void BlatterStressBalance::add_vars_to_output(const std::string &/*keyword*/, std::set<std::string> &result) {
  result.insert("u_sigma");
  result.insert("v_sigma");
}

//! Defines requested couplings fields.
PetscErrorCode BlatterStressBalance::define_variables(const std::set<std::string> &vars, const PIO &nc,
                                                      IO_Type nctype) {
  PetscErrorCode ierr;

  if (set_contains(vars, "u_sigma")) {
    ierr = u_sigma.define(nc, nctype); CHKERRQ(ierr);
  }

  if (set_contains(vars, "v_sigma")) {
    ierr = v_sigma.define(nc, nctype); CHKERRQ(ierr);
  }

  return 0;
}

//! Writes requested couplings fields to file.
PetscErrorCode BlatterStressBalance::write_variables(const std::set<std::string> &vars, const PIO &nc) {
  PetscErrorCode ierr;

  if (set_contains(vars, "u_sigma") || set_contains(vars, "v_sigma")) {
    ierr = save_velocity(); CHKERRQ(ierr);
  }

  if (set_contains(vars, "u_sigma")) {
    ierr = u_sigma.write(nc); CHKERRQ(ierr);
  }

  if (set_contains(vars, "v_sigma")) {
    ierr = v_sigma.write(nc); CHKERRQ(ierr);
  }

  return 0;
}

//! \brief Produce a report string for the standard output.
std::string BlatterStressBalance::stdout_report() {
  return "";
}


} // end of namespace stressbalance
} // end of namespace pism
