// Copyright (C) 2010, 2011, 2012 Constantine Khroulev and Ed Bueler
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

#include "SSB_Modifier.hh"
#include "flowlaws.hh"
#include "PISMVars.hh"
#include "IceGrid.hh"
#include "flowlaw_factory.hh"

PetscErrorCode SSB_Modifier::allocate() {
  PetscErrorCode ierr;

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

  ierr = Sigma.create(grid, "strainheat", false); CHKERRQ(ierr); // never diff'ed in hor dirs
  ierr = Sigma.set_attrs("internal",
                          "rate of strain heating in ice (dissipation heating)",
	        	  "W m-3", ""); CHKERRQ(ierr);
  ierr = Sigma.set_glaciological_units("mW m-3"); CHKERRQ(ierr);

  ierr = diffusive_flux.create(grid, "diffusive_flux", true, 1); CHKERRQ(ierr);
  ierr = diffusive_flux.set_attrs("internal", 
                                  "diffusive (SIA) flux components on the staggered grid",
                                  "", ""); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode SSB_Modifier::extend_the_grid(PetscInt old_Mz) {
  PetscErrorCode ierr;

  ierr =     u.extend_vertically(old_Mz, 0.0); CHKERRQ(ierr);
  ierr =     v.extend_vertically(old_Mz, 0.0); CHKERRQ(ierr);
  ierr = Sigma.extend_vertically(old_Mz, 0.0); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode SSBM_Trivial::init(PISMVars &vars) {
  PetscErrorCode ierr;

  ierr = SSB_Modifier::init(vars); CHKERRQ(ierr);

  enthalpy = dynamic_cast<IceModelVec3*>(vars.get("enthalpy"));
  if (enthalpy == NULL) SETERRQ(grid.com, 1, "enthalpy is not available");

  thickness = dynamic_cast<IceModelVec2S*>(vars.get("land_ice_thickness"));
  if (thickness == NULL) SETERRQ(grid.com, 1, "land_ice_thickness is not available");

  return 0;
}

SSBM_Trivial::SSBM_Trivial(IceGrid &g, EnthalpyConverter &e, const NCConfigVariable &c)
  : SSB_Modifier(g, e, c)
{
  IceFlowLawFactory ice_factory(grid.com, "", config, &EC);

  ice_factory.setType(config.get_string("sia_flow_law").c_str());

  ice_factory.setFromOptions();
  ice_factory.create(&flow_law);
}

SSBM_Trivial::~SSBM_Trivial()
{
  if (flow_law != NULL) {
    delete flow_law;
    flow_law = NULL;
  }
}


//! \brief Distribute the input velocity throughout the column.
/*!
 * Things to update:
 * - 3D-distributed horizontal velocity
 * - maximum horizontal velocity
 * - diffusive ice flux
 * - maximum diffusivity
 * - strain heating (Sigma)
 */
PetscErrorCode SSBM_Trivial::update(IceModelVec2V *vel_input,
                                    IceModelVec2S *D2_input,
                                    bool fast) {
  PetscErrorCode ierr;

  if (fast)
    return 0;

  // horizontal velocity and its maximum:
  ierr = u.begin_access(); CHKERRQ(ierr);
  ierr = v.begin_access(); CHKERRQ(ierr);
  ierr = vel_input->begin_access(); CHKERRQ(ierr);
  PetscReal my_u_max = 0, my_v_max = 0;
  for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
      ierr = u.setColumn(i,j, (*vel_input)(i,j).u); CHKERRQ(ierr);
      ierr = v.setColumn(i,j, (*vel_input)(i,j).v); CHKERRQ(ierr);

      my_u_max = PetscMax(my_u_max, PetscAbs((*vel_input)(i,j).u));
      my_v_max = PetscMax(my_v_max, PetscAbs((*vel_input)(i,j).u));
    }
  }
  ierr = vel_input->end_access(); CHKERRQ(ierr);
  ierr = v.end_access(); CHKERRQ(ierr);
  ierr = u.end_access(); CHKERRQ(ierr);  

  // Communicate to get ghosts (needed to compute w):
  ierr = u.beginGhostComm(); CHKERRQ(ierr);
  ierr = v.beginGhostComm(); CHKERRQ(ierr);
  ierr = u.endGhostComm(); CHKERRQ(ierr);
  ierr = v.endGhostComm(); CHKERRQ(ierr);

  ierr = PISMGlobalMax(&my_u_max, &u_max, grid.com); CHKERRQ(ierr);
  ierr = PISMGlobalMax(&my_v_max, &v_max, grid.com); CHKERRQ(ierr);

  // diffusive flux and maximum diffusivity
  ierr = diffusive_flux.set(0.0); CHKERRQ(ierr);
  D_max = 0.0;

  // strain heating
  ierr = compute_sigma(D2_input, Sigma); CHKERRQ(ierr);

  return 0;
}


//! \brief Compute the volumetric strain heating.
/*!
 * Uses:
 * - delta on the staggered grid, which should be initialized by the update(true) call.
 * - enthalpy
 * - surface gradient on the staggered grid
 * - ice thickness relative to the smoothed bed
 */
PetscErrorCode SSBM_Trivial::compute_sigma(IceModelVec2S *D2_input, IceModelVec3 &result) {
  PetscErrorCode ierr;
  PetscScalar *E, *sigma;
  const PetscReal
    n_glen  = flow_law->exponent(),
    Sig_pow = (1.0 + n_glen) / (2.0 * n_glen),
    enhancement_factor = config.get("ssa_enhancement_factor"),
    standard_gravity = config.get("standard_gravity"),
    ice_rho = config.get("ice_density");

  ierr = enthalpy->begin_access(); CHKERRQ(ierr);
  ierr = result.begin_access(); CHKERRQ(ierr);
  ierr = thickness->begin_access(); CHKERRQ(ierr);
  ierr = D2_input->begin_access(); CHKERRQ(ierr);

  for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
      ierr = enthalpy->getInternalColumn(i,j,&E); CHKERRQ(ierr);
      ierr = result.getInternalColumn(i,j,&sigma); CHKERRQ(ierr);

      PetscReal thk = (*thickness)(i,j);
      // in the ice:
      PetscInt ks = grid.kBelowHeight(thk);
        for (PetscInt k=0; k<ks; ++k) {
          // Use hydrostatic pressure; presumably this is not quite right in context
          //   of shelves and streams; here we hard-wire the Glen law
          PetscReal depth = thk - grid.zlevels[k],
            pressure = ice_rho * standard_gravity * depth, // FIXME issue #15
          // Account for the enhancement factor.
          //   Note, enhancement factor is not used in SSA anyway.
          //   Should we get rid of it completely?  If not, what is most consistent here?
            BofT    = flow_law->hardness_parameter(E[k], pressure) * pow(enhancement_factor,-1/n_glen);
          sigma[k] = 2.0 * BofT * pow((*D2_input)(i,j), Sig_pow);
        }

        // above the ice:
        for (PetscInt k=ks+1; k<grid.Mz; ++k) {
          sigma[k] = 0.0;
        }
    }
  }

  ierr = D2_input->end_access(); CHKERRQ(ierr);
  ierr = thickness->end_access(); CHKERRQ(ierr);
  ierr = result.end_access(); CHKERRQ(ierr);  
  ierr = enthalpy->end_access(); CHKERRQ(ierr);

  return 0;
}
