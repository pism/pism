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

#include "PISMStressBalance.hh"
#include "ShallowStressBalance.hh"
#include "SSB_Modifier.hh"
#include "PISMOcean.hh"
#include "IceGrid.hh"

PISMStressBalance::PISMStressBalance(IceGrid &g,
                                     ShallowStressBalance *sb,
                                     SSB_Modifier *ssb_mod,
                                     PISMOceanModel *ocean_model,
                                     const NCConfigVariable &conf)
  : PISMComponent_Diag(g, conf), stress_balance(sb), modifier(ssb_mod), ocean(ocean_model) {

  basal_melt_rate = NULL;
  variables = NULL;

  allocate();
}

PISMStressBalance::~PISMStressBalance() {
  delete stress_balance;
  delete modifier;
}

PetscErrorCode PISMStressBalance::allocate() {
  PetscErrorCode ierr;

  // allocate the vertical velocity field:
  ierr = w.create(grid, "wvel_rel", false); CHKERRQ(ierr);
  ierr = w.set_attrs("diagnostic",
                     "vertical velocity of ice, relative to base of ice directly below",
                     "m s-1", ""); CHKERRQ(ierr);
  w.time_independent = false;
  ierr = w.set_glaciological_units("m year-1"); CHKERRQ(ierr);
  w.write_in_glaciological_units = true;

  return 0;
}

//! \brief Initialize the PISMStressBalance object.
PetscErrorCode PISMStressBalance::init(PISMVars &vars) {
  PetscErrorCode ierr;

  variables = &vars;

  ierr = stress_balance->init(vars); CHKERRQ(ierr);   
  ierr = modifier->init(vars); CHKERRQ(ierr); 

  return 0;
}

PetscErrorCode PISMStressBalance::set_boundary_conditions(IceModelVec2Int &locations,
                                                          IceModelVec2V &velocities) {
  PetscErrorCode ierr;
  ierr = stress_balance->set_boundary_conditions(locations, velocities); CHKERRQ(ierr);
  return 0;
}

//! \brief Set the basal melt rate. (If not NULL, it will be included in the
//! computation of the vertical valocity).
PetscErrorCode PISMStressBalance::set_basal_melt_rate(IceModelVec2S *bmr_input) {
  basal_melt_rate = bmr_input;
  return 0;
}

//! \brief Performs the shallow stress balance computation.
PetscErrorCode PISMStressBalance::update(bool fast) {
  PetscErrorCode ierr;
  IceModelVec2V *velocity_2d;
  IceModelVec2S *D2;
  IceModelVec3  *u, *v;

  // Tell the ShallowStressBalance object about the current sea level:
  if (ocean) {
    PetscReal sea_level;
    ierr = ocean->sea_level_elevation(sea_level); CHKERRQ(ierr);
    stress_balance->set_sea_level_elevation(sea_level);
  }

  ierr = stress_balance->update(fast); CHKERRQ(ierr);

  ierr = stress_balance->get_advective_2d_velocity(velocity_2d); CHKERRQ(ierr); 
  ierr = stress_balance->get_D2(D2); CHKERRQ(ierr);

  ierr = modifier->update(velocity_2d, D2, fast); CHKERRQ(ierr);

  if (!fast) {
    ierr = modifier->get_horizontal_3d_velocity(u, v); CHKERRQ(ierr);

    ierr = compute_vertical_velocity(u, v, basal_melt_rate, w); CHKERRQ(ierr); 
  }

  return 0;
}

PetscErrorCode PISMStressBalance::get_advective_2d_velocity(IceModelVec2V* &result) {
  PetscErrorCode ierr;
  ierr = stress_balance->get_advective_2d_velocity(result); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode PISMStressBalance::get_diffusive_flux(IceModelVec2Stag* &result) {
  PetscErrorCode ierr;
  ierr = modifier->get_diffusive_flux(result); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode PISMStressBalance::get_max_diffusivity(PetscReal &D) {
  PetscErrorCode ierr;
  ierr = modifier->get_max_diffusivity(D); CHKERRQ(ierr);  
  return 0;
}

PetscErrorCode PISMStressBalance::get_max_2d_velocity(PetscReal &u_max, PetscReal &v_max) {
  PetscErrorCode ierr;
  ierr = stress_balance->get_max_2d_velocity(u_max, v_max); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode PISMStressBalance::get_3d_velocity(IceModelVec3* &u, IceModelVec3* &v, IceModelVec3* &w_out) {
  PetscErrorCode ierr;
  ierr = modifier->get_horizontal_3d_velocity(u, v); CHKERRQ(ierr);
  w_out = &w;
  return 0;
}

PetscErrorCode PISMStressBalance::get_max_3d_velocity(PetscReal &u, PetscReal &v, PetscReal &w_out) {
  PetscErrorCode ierr;
  ierr = modifier->get_max_horizontal_velocity(u, v); CHKERRQ(ierr);
  w_out = w_max;
  return 0;
}

PetscErrorCode PISMStressBalance::get_basal_frictional_heating(IceModelVec2S* &result) {
  PetscErrorCode ierr;
  ierr = stress_balance->get_basal_frictional_heating(result); CHKERRQ(ierr); 
  return 0;
}

PetscErrorCode PISMStressBalance::get_volumetric_strain_heating(IceModelVec3* &result) {
  PetscErrorCode ierr;
  ierr = modifier->get_volumetric_strain_heating(result); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode PISMStressBalance::get_principal_strain_rates(
                IceModelVec2S &result_e1, IceModelVec2S &result_e2) {
  PetscErrorCode ierr;
  ierr = stress_balance->compute_principal_strain_rates(result_e1, result_e2); CHKERRQ(ierr);
  return 0;
}


//! \brief Extend the grid vertically.
PetscErrorCode PISMStressBalance::extend_the_grid(PetscInt old_Mz) {
  PetscErrorCode ierr;

  ierr = w.extend_vertically(old_Mz, 0.0); CHKERRQ(ierr);

  ierr = stress_balance->extend_the_grid(old_Mz); CHKERRQ(ierr);

  ierr = modifier->extend_the_grid(old_Mz); CHKERRQ(ierr);

  return 0;
}

//! Compute vertical velocity using incompressibility of the ice.
/*!
The vertical velocity \f$w(x,y,z,t)\f$ is the velocity <i>relative to the
location of the base of the ice column</i>.  That is, the vertical velocity
computed here is identified as \f$\tilde w(x,y,s,t)\f$ in the page
\ref vertchange.

Thus \f$w<0\f$ here means that that
that part of the ice is getting closer to the base of the ice, and so on.
The slope of the bed (i.e. relative to the geoid) and/or the motion of the
bed (i.e. from bed deformation) do not affect the vertical velocity.

In fact the following statement is exactly true if the basal melt rate is zero:
the vertical velocity at a point in the ice is positive (negative) if and only
if the average horizontal divergence of the horizontal velocity, in the portion
of the ice column below that point, is negative (positive).
In particular, because \f$z=0\f$ is the location of the base of the ice
always, the only way to have \f$w(x,y,0,t) \ne 0\f$ is to have a basal melt
rate.

Incompressibility itself says
   \f[ \nabla\cdot\mathbf{U} + \frac{\partial w}{\partial z} = 0. \f]
This is immediately equivalent to the integral
   \f[ w(x,y,z,t) = - \int_{b(x,y,t)}^{z} \nabla\cdot\mathbf{U}\,d\zeta
                           + w_b(x,y,t). \f]
Here the value \f$w_b(x,y,t)\f$ is either zero or the negative of the basal melt rate
according to the value of the flag \c include_bmr_in_continuity.

The vertical integral is computed by the trapezoid rule.
 */
PetscErrorCode PISMStressBalance::compute_vertical_velocity(IceModelVec3 *u, IceModelVec3 *v,
                                                            IceModelVec2S *bmr,
                                                            IceModelVec3 &result) {
  PetscErrorCode ierr;
  const PetscScalar dx = grid.dx, dy = grid.dy;

  ierr = u->begin_access(); CHKERRQ(ierr);
  ierr = v->begin_access(); CHKERRQ(ierr);
  ierr = result.begin_access(); CHKERRQ(ierr);

  if (bmr) {
    ierr = bmr->begin_access(); CHKERRQ(ierr);
  }

  PetscScalar *w_ij, *u_im1, *u_ip1, *v_jm1, *v_jp1;

  PetscReal my_w_max = 0.0;
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      ierr = result.getInternalColumn(i,j,&w_ij); CHKERRQ(ierr);

      ierr = u->getInternalColumn(i-1,j,&u_im1); CHKERRQ(ierr);
      ierr = u->getInternalColumn(i+1,j,&u_ip1); CHKERRQ(ierr);

      ierr = v->getInternalColumn(i,j-1,&v_jm1); CHKERRQ(ierr);
      ierr = v->getInternalColumn(i,j+1,&v_jp1); CHKERRQ(ierr);

      // at the base:
      if (bmr) {
        w_ij[0] = - (*bmr)(i,j);
      } else {
        w_ij[0] = 0.0;
      }
      my_w_max = PetscMax(my_w_max, PetscAbs(w_ij[0]));
      
      // within the ice and above:
      PetscScalar OLDintegrand
             = (u_ip1[0] - u_im1[0]) / (2.0*dx) + (v_jp1[0] - v_jm1[0]) / (2.0*dy);
      for (PetscInt k = 1; k < grid.Mz; ++k) {
        const PetscScalar NEWintegrand
             = (u_ip1[k] - u_im1[k]) / (2.0*dx) + (v_jp1[k] - v_jm1[k]) / (2.0*dy);
        const PetscScalar dz = grid.zlevels[k] - grid.zlevels[k-1];
        w_ij[k] = w_ij[k-1] - 0.5 * (NEWintegrand + OLDintegrand) * dz;
        OLDintegrand = NEWintegrand;

        my_w_max = PetscMax(my_w_max, PetscAbs(w_ij[k]));
      }
    }
  }

  if (bmr) {
    ierr = bmr->end_access(); CHKERRQ(ierr);
  }

  ierr = u->end_access(); CHKERRQ(ierr);
  ierr = v->end_access(); CHKERRQ(ierr);
  ierr = result.end_access(); CHKERRQ(ierr);

  ierr = PISMGlobalMax(&my_w_max, &w_max, grid.com); CHKERRQ(ierr);
  
  return 0;
}

PetscErrorCode PISMStressBalance::stdout_report(string &result) {
  PetscErrorCode ierr;
  string tmp1, tmp2;
  
  ierr = stress_balance->stdout_report(tmp1); CHKERRQ(ierr);

  ierr = modifier->stdout_report(tmp2); CHKERRQ(ierr);

  result = tmp1 + tmp2;

  return 0;
}

PetscErrorCode PISMStressBalance::define_variables(set<string> vars, const PIO &nc,
                                                   PISM_IO_Type nctype) {
  PetscErrorCode ierr;

  ierr = stress_balance->define_variables(vars, nc, nctype); CHKERRQ(ierr);
  ierr = modifier->define_variables(vars, nc, nctype); CHKERRQ(ierr);

  return 0;
}


PetscErrorCode PISMStressBalance::write_variables(set<string> vars, string filename) {
  PetscErrorCode ierr;

  ierr = stress_balance->write_variables(vars, filename); CHKERRQ(ierr);
  ierr = modifier->write_variables(vars, filename); CHKERRQ(ierr);

  return 0;
}

void PISMStressBalance::add_vars_to_output(string keyword, map<string,NCSpatialVariable> &result) {

  stress_balance->add_vars_to_output(keyword, result);
  modifier->add_vars_to_output(keyword, result);

}

