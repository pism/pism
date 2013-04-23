// Copyright (C) 2010-2013 Ed Bueler and Constantine Khroulev
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
#include "PISMOcean.hh"
#include "IceGrid.hh"
#include "PISMVars.hh"
#include "basal_resistance.hh"
#include "FE3DTools.h"
#include "enthalpyConverter.hh"

/*!
 * FIXMEs:
 *
 * We need to allow spatially-variable (and depth-dependent) ice hardness.
 *
 * We need to compute the volumetric strain heating.
 *
 */

//! C-wrapper for PISM's IceFlowLaw::viscosity().
void viscosity(void *ctx, PetscReal hardness, PetscReal gamma,
	       PetscReal *eta, PetscReal *deta) {
  BlatterQ1Ctx *blatter_ctx = (BlatterQ1Ctx*)ctx;
  BlatterStressBalance *blatter_stress_balance = (BlatterStressBalance*)blatter_ctx->extra;

  blatter_stress_balance->flow_law->effective_viscosity(hardness, gamma, eta, deta);
}

//! C-wrapper for PISM's IceBasalResistancePlasticLaw::dragWithDerivative().
void drag(void *ctx, PetscReal tauc, PetscReal u, PetscReal v,
	  PetscReal *taud, PetscReal *dtaub) {
  BlatterQ1Ctx *blatter_ctx = (BlatterQ1Ctx*)ctx;
  BlatterStressBalance *blatter_stress_balance = (BlatterStressBalance*)blatter_ctx->extra;

  blatter_stress_balance->basal.dragWithDerivative(tauc, u, v, taud, dtaub);
}

BlatterStressBalance::BlatterStressBalance(IceGrid &g,
					   IceBasalResistancePlasticLaw &b,
					   EnthalpyConverter &e,
					   const NCConfigVariable &conf)
  : ShallowStressBalance(g, b, e, conf), min_thickness(10.0)
{
  if (allocate_blatter() != 0) {
    PetscPrintf(grid.com, "FATAL ERROR: BlatterStressBalance allocation failed.\n");
    PISMEnd();
  }
}

BlatterStressBalance::~BlatterStressBalance()
{
  if (deallocate_blatter() != 0) {
    PetscPrintf(grid.com, "FATAL ERROR: BlatterStressBalance deallocation failed.\n");
    PISMEnd();
  }
}

PetscErrorCode BlatterStressBalance::allocate_blatter() {
  PetscErrorCode ierr;
  DM da2;
  int blatter_Mz = (int)config.get("blatter_Mz");
  ierr = grid.get_dm(1, grid.max_stencil_width, da2); CHKERRQ(ierr);

  ctx.Lx = 2.0 * grid.Lx;
  ctx.Ly = 2.0 * grid.Ly;
  ctx.dirichlet_scale = 1.0;
  ctx.rhog = config.get("ice_density") * config.get("standard_gravity");
  ctx.no_slip = PETSC_TRUE;	// FIXME (at least make configurable)
  ctx.nonlinear.viscosity = viscosity;
  ctx.nonlinear.drag = drag;
  ctx.extra = this;
  initialize_Q12D(ctx.Q12D.chi, ctx.Q12D.dchi);
  initialize_Q13D(ctx.Q13D.chi, ctx.Q13D.dchi);

  ierr = BlatterQ1_create(grid.com, da2, blatter_Mz, &this->ctx, &this->snes); CHKERRQ(ierr);

  // now allocate u, v, and strain_heating (strain heating)
  ierr =     u.create(grid, "uvel", true); CHKERRQ(ierr);
  ierr =     u.set_attrs("diagnostic", "horizontal velocity of ice in the X direction",
			  "m s-1", "land_ice_x_velocity"); CHKERRQ(ierr);
  ierr =     u.set_glaciological_units("m year-1"); CHKERRQ(ierr);
  u.write_in_glaciological_units = true;

  ierr =     v.create(grid, "vvel", true); CHKERRQ(ierr);
  ierr =     v.set_attrs("diagnostic", "horizontal velocity of ice in the Y direction",
			  "m s-1", "land_ice_y_velocity"); CHKERRQ(ierr);
  ierr =     v.set_glaciological_units("m year-1"); CHKERRQ(ierr);
  v.write_in_glaciological_units = true;

  // strain_heating
  ierr = strain_heating.create(grid, "strainheat", false); CHKERRQ(ierr); // never diff'ed in hor dirs
  ierr = strain_heating.set_attrs("internal",
                          "rate of strain heating in ice (dissipation heating)",
	        	  "W m-3", ""); CHKERRQ(ierr);
  ierr = strain_heating.set_glaciological_units("mW m-3"); CHKERRQ(ierr);

  // storage for u and v on the sigma vertical grid (for restarting)
  ierr =     u_sigma.create(grid, "uvel_sigma", false, blatter_Mz); CHKERRQ(ierr);
  ierr =     u_sigma.set_attrs("diagnostic",
			       "horizontal velocity of ice in the X direction on the sigma vertical grid",
			       "m s-1", ""); CHKERRQ(ierr);
  ierr =     u_sigma.set_glaciological_units("m year-1"); CHKERRQ(ierr);
  u_sigma.write_in_glaciological_units = false;

  ierr =     v_sigma.create(grid, "vvel_sigma", false, blatter_Mz); CHKERRQ(ierr);
  ierr =     v_sigma.set_attrs("diagnostic",
			       "horizontal velocity of ice in the Y direction on the sigma vertical grid",
			       "m s-1", ""); CHKERRQ(ierr);
  ierr =     v_sigma.set_glaciological_units("m year-1"); CHKERRQ(ierr);
  v_sigma.write_in_glaciological_units = false;

  {
    IceFlowLawFactory ice_factory(grid.com, "blatter_",
                                  grid.get_unit_system(),
                                  config, &EC);
    ice_factory.removeType(ICE_GOLDSBY_KOHLSTEDT);

    ierr = ice_factory.setType(config.get_string("blatter_flow_law")); CHKERRQ(ierr);

    ierr = ice_factory.setFromOptions(); CHKERRQ(ierr);
    ierr = ice_factory.create(&flow_law); CHKERRQ(ierr);
  }

  return 0;
}

PetscErrorCode BlatterStressBalance::deallocate_blatter() {
  PetscErrorCode ierr;

  ierr = SNESDestroy(&this->snes);CHKERRQ(ierr);

  delete flow_law;

  return 0;
}


PetscErrorCode BlatterStressBalance::init(PISMVars &vars) {

  variables = &vars;

  bed_elevation = dynamic_cast<IceModelVec2S*>(vars.get("bedrock_altitude"));
  if (bed_elevation == NULL) SETERRQ(grid.com, 1, "bedrock_altitude is not available");

  ice_thickness = dynamic_cast<IceModelVec2S*>(vars.get("land_ice_thickness"));
  if (ice_thickness == NULL) SETERRQ(grid.com, 1, "land_ice_thickness is not available");

  tauc = dynamic_cast<IceModelVec2S*>(vars.get("tauc"));
  if (tauc == NULL) SETERRQ(grid.com, 1, "tauc is not available");

  enthalpy = dynamic_cast<IceModelVec3*>(vars.get("enthalpy"));
  if (enthalpy == NULL) SETERRQ(grid.com, 1, "enthalpy is not available");

  return 0;
}

PetscErrorCode BlatterStressBalance::update(bool fast) {
  PetscErrorCode ierr;

  if (fast) {
    ierr = verbPrintf(1,grid.com,
       "PISM ERROR:  'fast' mode not meaningful for BlatterStressBalance\n"
       "  ENDING ...\n\n"); CHKERRQ(ierr);
    PISMEnd();
  }

  // setup
  ierr =  setup(); CHKERRQ(ierr);

  // solve
  ierr = SNESSolve(this->snes, PETSC_NULL, PETSC_NULL);CHKERRQ(ierr);

  // Transfer solution from the FEM mesh to the regular grid used in the rest
  // of PISM and compute the vertically-averaged velocity.
  ierr = transfer_velocity(); CHKERRQ(ierr);

  ierr =  compute_volumetric_strain_heating(); CHKERRQ(ierr);

  return 0;
}

/*! \brief Set up model parameters on the fine grid. */
/*!
 * This method expects bed_elevation, ice_thickness, and tauc to have width=1 ghosts.
 *
 * We should also compute ice hardness on the "sigma" grid here.
 */
PetscErrorCode BlatterStressBalance::setup() {
  PetscErrorCode ierr;
  PrmNode **parameters;
  DM da;
  PetscReal
    ice_density = config.get("ice_density"),
    sea_water_density = config.get("sea_water_density"),
    alpha = ice_density / sea_water_density;

  ierr = SNESGetDM(this->snes, &da); CHKERRQ(ierr);

  ierr = BlatterQ1_begin_2D_parameter_access(da, &parameters); CHKERRQ(ierr);

  ierr = bed_elevation->begin_access(); CHKERRQ(ierr);
  ierr = ice_thickness->begin_access(); CHKERRQ(ierr);
  ierr = tauc->begin_access(); CHKERRQ(ierr);

  PetscInt GHOSTS = 1;
  for (PetscInt   i = grid.xs - GHOSTS; i < grid.xs+grid.xm + GHOSTS; ++i) {
    for (PetscInt j = grid.ys - GHOSTS; j < grid.ys+grid.ym + GHOSTS; ++j) {

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
  }

  ierr = bed_elevation->end_access(); CHKERRQ(ierr);
  ierr = ice_thickness->end_access(); CHKERRQ(ierr);
  ierr = tauc->end_access(); CHKERRQ(ierr);

  ierr = BlatterQ1_end_2D_parameter_access(da, &parameters); CHKERRQ(ierr);

  return 0;
}

//! Initialize ice hardness on the "sigma" grid.
PetscErrorCode BlatterStressBalance::initialize_ice_hardness() {
  PetscErrorCode ierr;
  PetscScalar ***hardness, *E;
  int Mz_fem = static_cast<int>(config.get("blatter_Mz"));
  DM da;


  ierr = SNESGetDM(this->snes, &da); CHKERRQ(ierr);

  ierr = BlatterQ1_begin_hardness_access(da, &hardness); CHKERRQ(ierr);
  ierr = enthalpy->begin_access(); CHKERRQ(ierr);
  ierr = ice_thickness->begin_access(); CHKERRQ(ierr);

  PetscInt GHOSTS = 1;
  for (PetscInt   i = grid.xs - GHOSTS; i < grid.xs+grid.xm + GHOSTS; ++i) {
    for (PetscInt j = grid.ys - GHOSTS; j < grid.ys+grid.ym + GHOSTS; ++j) {
      PetscScalar thk = (*ice_thickness)(i,j);

      // fudge ice thickness (FIXME!!!)
      if (thk < min_thickness)
	thk += min_thickness;

      PetscScalar dz_fem = thk / (Mz_fem - 1);
      ierr = enthalpy->getInternalColumn(i, j, &E); CHKERRQ(ierr);

      // compute ice hardness on the sigma grid
      for (int k = 0; k < Mz_fem; ++k) {
	PetscReal z_fem = k * dz_fem,
	  depth = thk - z_fem,
	  pressure = EC.getPressureFromDepth(depth),
	  E_local;
	int k0 = grid.kBelowHeight(z_fem);
	
	if (k0 + 1 < grid.Mz) {
	  PetscReal lambda = (z_fem - grid.zlevels[k0]) / (grid.zlevels[k0+1] - grid.zlevels[k0]);
	  
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
PetscErrorCode BlatterStressBalance::transfer_velocity() {
  PetscErrorCode ierr;

  PISMVector2 ***U;
  PetscScalar *u_ij, *v_ij;
  DM da;
  Vec X;
  int Mz_fem = static_cast<int>(config.get("blatter_Mz"));

  ierr = SNESGetDM(this->snes, &da); CHKERRQ(ierr);
  ierr = SNESGetSolution(this->snes, &X); CHKERRQ(ierr);

  ierr = DMDAVecGetArray(da, X, &U); CHKERRQ(ierr);

  ierr = u.begin_access(); CHKERRQ(ierr);
  ierr = v.begin_access(); CHKERRQ(ierr);
  ierr = ice_thickness->begin_access(); CHKERRQ(ierr);
  ierr = m_velocity.begin_access(); CHKERRQ(ierr);
 
  for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
      ierr = u.getInternalColumn(i, j, &u_ij); CHKERRQ(ierr);
      ierr = v.getInternalColumn(i, j, &v_ij); CHKERRQ(ierr);

      PetscReal thk = (*ice_thickness)(i,j);

      // fudge ice thickness (FIXME!!!)
      if (thk < min_thickness)
	thk += min_thickness;

      PetscInt current_level = 0;
      PetscScalar dz_fem = thk / (Mz_fem - 1);

      // compute vertically-averaged velocity using trapezoid rule
      PetscReal ubar = 0, vbar = 0;
      for (int k = 0; k < Mz_fem - 1; ++k) {
	ubar += U[i][j][k].u + U[i][j][k+1].u;
	vbar += U[i][j][k].v + U[i][j][k+1].v;
      }
      // finish the traperoidal rule (1/2 * dz) and compute the average:
      m_velocity(i,j).u = ubar * (0.5*dz_fem) / thk;
      m_velocity(i,j).v = vbar * (0.5*dz_fem) / thk;
      
      // compute 3D horizontal velocity
      for (int k = 0; k < grid.Mz; ++k) {

	// find the FEM grid level just below the current PISM grid level
	while ((current_level + 1) * dz_fem < grid.zlevels[k])
	  current_level++;

	if (current_level + 1 < Mz_fem) {
	  // use linear interpolation
	  PetscReal z0 = current_level * dz_fem,
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
    } // j-loop
  } // i-loop

  ierr = m_velocity.end_access(); CHKERRQ(ierr);
  ierr = ice_thickness->end_access(); CHKERRQ(ierr);
  ierr = v.end_access(); CHKERRQ(ierr);
  ierr = u.end_access(); CHKERRQ(ierr);

  ierr = DMDAVecRestoreArray(da, X, &U); CHKERRQ(ierr);

  ierr = u.update_ghosts(); CHKERRQ(ierr);
  ierr = v.update_ghosts(); CHKERRQ(ierr);

  return 0;
}

//! Copy velocity from a dof=2 vector to special storage (to save it for re-starting).
PetscErrorCode BlatterStressBalance::save_velocity() {
  PetscErrorCode ierr;

  PISMVector2 ***U;
  PetscScalar *u_ij, *v_ij;
  DM da;
  Vec X;
  int Mz_fem = static_cast<int>(config.get("blatter_Mz"));

  ierr = SNESGetDM(this->snes, &da); CHKERRQ(ierr);
  ierr = SNESGetSolution(this->snes, &X); CHKERRQ(ierr);

  ierr = DMDAVecGetArray(da, X, &U); CHKERRQ(ierr);

  ierr = u_sigma.begin_access(); CHKERRQ(ierr);
  ierr = v_sigma.begin_access(); CHKERRQ(ierr);
 
  for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
      ierr = u_sigma.getInternalColumn(i, j, &u_ij); CHKERRQ(ierr);
      ierr = v_sigma.getInternalColumn(i, j, &v_ij); CHKERRQ(ierr);

      for (int k = 0; k < Mz_fem; ++k) {
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

PetscErrorCode BlatterStressBalance::extend_the_grid(PetscInt old_Mz) {
  PetscErrorCode ierr;
  ierr = u.extend_vertically(old_Mz, 0.0); CHKERRQ(ierr);
  ierr = v.extend_vertically(old_Mz, 0.0); CHKERRQ(ierr);
  ierr = strain_heating.extend_vertically(old_Mz, 0.0); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode BlatterStressBalance::compute_volumetric_strain_heating() {
  PetscErrorCode ierr;
  // FIXME: implement this
  ierr = strain_heating.set(0.0); CHKERRQ(ierr);
  
  return 0;
}

void BlatterStressBalance::add_vars_to_output(string /*keyword*/,
					      map<string,NCSpatialVariable> &result) {
  result["u_sigma"] = u_sigma.get_metadata();
  result["v_sigma"] = v_sigma.get_metadata();
}

//! Defines requested couplings fields.
PetscErrorCode BlatterStressBalance::define_variables(set<string> vars, const PIO &nc,
						      PISM_IO_Type nctype) {
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
PetscErrorCode BlatterStressBalance::write_variables(set<string> vars, const PIO &nc) {
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
