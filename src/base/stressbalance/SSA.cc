// Copyright (C) 2004--2011 Constantine Khroulev, Ed Bueler and Jed Brown
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

#include "SSA.hh"

SSA::SSA(IceGrid &g, IceBasalResistancePlasticLaw &b,
         IceFlowLaw &i, EnthalpyConverter &e,
         const NCConfigVariable &c)
  : ShallowStressBalance(g, b, i, e, c)
{
  mask = NULL;
  thickness = NULL;
  tauc = NULL;
  surface = NULL;
  bed = NULL;
  enthalpy = NULL;

  allocate();
}


//! \brief Initialize a generic regular-grid SSA solver.
PetscErrorCode SSA::init(PISMVars &vars) {
  PetscErrorCode ierr;

  ierr = ShallowStressBalance::init(vars); CHKERRQ(ierr);

  ierr = verbPrintf(2,grid.com,"* Initializing the SSA stress balance...\n"); CHKERRQ(ierr);

  mask = dynamic_cast<IceModelVec2Mask*>(vars.get("mask"));
  if (mask == NULL) SETERRQ(1, "mask is not available");

  thickness = dynamic_cast<IceModelVec2S*>(vars.get("land_ice_thickness"));
  if (thickness == NULL) SETERRQ(1, "land_ice_thickness is not available");

  tauc = dynamic_cast<IceModelVec2S*>(vars.get("tauc"));
  if (tauc == NULL) SETERRQ(1, "tauc is not available");

  surface = dynamic_cast<IceModelVec2S*>(vars.get("surface_altitude"));
  if (surface == NULL) SETERRQ(1, "surface_altitude is not available");

  bed = dynamic_cast<IceModelVec2S*>(vars.get("bedrock_altitude"));
  if (bed == NULL) SETERRQ(1, "bedrock_altitude is not available");

  enthalpy = dynamic_cast<IceModelVec3*>(vars.get("enthalpy"));
  if (enthalpy == NULL) SETERRQ(1, "enthalpy is not available");


  // Check if PISM is being initialized from an output file from a previous run
  // and read the initial guess (unless asked not to).
  bool i_set;
  string filename;
  ierr = PISMOptionsString("-i", "PISM input file",
                           filename, i_set); CHKERRQ(ierr);

  if (i_set) {
    bool dont_read_initial_guess, ubar_ssa_found, vbar_ssa_found;
    int start;
    NCTool nc(grid.com, grid.rank);

    ierr = PISMOptionsIsSet("-dontreadSSAvels", dont_read_initial_guess); CHKERRQ(ierr);

    ierr = nc.open_for_reading(filename.c_str()); CHKERRQ(ierr);
    ierr = nc.find_variable("ubar_ssa", NULL, ubar_ssa_found); CHKERRQ(ierr); 
    ierr = nc.find_variable("vbar_ssa", NULL, vbar_ssa_found); CHKERRQ(ierr); 
    ierr = nc.get_dim_length("t", &start); CHKERRQ(ierr);
    ierr = nc.close(); CHKERRQ(ierr); 
    start -= 1;

    if (ubar_ssa_found && vbar_ssa_found &&
        (! dont_read_initial_guess)) {
      ierr = verbPrintf(3,grid.com,"Reading ubar_ssa and vbar_ssa...\n"); CHKERRQ(ierr);

      ierr = velocity.read(filename.c_str(), start); CHKERRQ(ierr); 
    }
    
  } else {
    ierr = velocity.set(0.0); CHKERRQ(ierr); // default initial guess
  }

  event_ssa = grid.profiler->create("ssa_update", "time spent solving the SSA");

  return 0;
}

//! \brief Allocate objects which any SSA solver would use.
PetscErrorCode SSA::allocate() {
  PetscErrorCode ierr;

  const PetscScalar power = 1.0 / ice.exponent();
  char unitstr[TEMPORARY_STRING_LENGTH];
  snprintf(unitstr, sizeof(unitstr), "Pa s%f", power);
  ierr = hardness.create(grid, "hardness", false); CHKERRQ(ierr);
  ierr = hardness.set_attrs("diagnostic",
                            "vertically-averaged ice hardness",
                            unitstr, ""); CHKERRQ(ierr);

  ierr = nuH.create(grid, "nuH", true); CHKERRQ(ierr);
  ierr = nuH.set_attrs("internal",
                       "ice thickness times effective viscosity",
                       "Pa s m", ""); CHKERRQ(ierr);

  ierr = nuH_old.create(grid, "nuH_old", true); CHKERRQ(ierr);
  ierr = nuH_old.set_attrs("internal",
                           "ice thickness times effective viscosity (before an update)",
                           "Pa s m", ""); CHKERRQ(ierr);

  ierr = taud.create(grid, "taud", false); CHKERRQ(ierr);
  ierr = taud.set_attrs("diagnostic",
                        "X-component of the driving shear stress at the base of ice",
                        "Pa", "", 0); CHKERRQ(ierr);
  ierr = taud.set_attrs("diagnostic",
                        "Y-component of the driving shear stress at the base of ice",
                        "Pa", "", 1); CHKERRQ(ierr);

  ierr = velocity_old.create(grid, "velocity_old", true); CHKERRQ(ierr);
  ierr = velocity_old.set_attrs("internal",
                                "old SSA velocity field; used for re-trying with a different epsilon",
                                "m s-1", ""); CHKERRQ(ierr);

  // override velocity metadata
  ierr = velocity.set_name("bar_ssa"); CHKERRQ(ierr);
  ierr = velocity.set_attrs("internal_restart", "SSA model ice velocity in the X direction",
                            "m s-1", "", 0); CHKERRQ(ierr);

  ierr = velocity.set_attrs("internal_restart", "SSA model ice velocity in the Y direction",
                            "m s-1", "", 1); CHKERRQ(ierr);

  ierr = velocity.set_glaciological_units("m year-1"); CHKERRQ(ierr);

  // mimic IceGrid::createDA() with TRANSPOSE :
  PetscInt dof=2, stencil_width=1;
  ierr = DACreate2d(grid.com, DA_XYPERIODIC, DA_STENCIL_BOX,
                    grid.My, grid.Mx,
                    grid.Ny, grid.Nx,
                    dof, stencil_width,
                    grid.procs_y, grid.procs_x,
                    &SSADA); CHKERRQ(ierr);

  ierr = DACreateGlobalVector(SSADA, &SSAX); CHKERRQ(ierr);
  
  return 0;
}


PetscErrorCode SSA::deallocate() {
  PetscErrorCode ierr;

  if (SSAX != PETSC_NULL) {
    ierr = VecDestroy(SSAX); CHKERRQ(ierr);
  }

  if (SSADA != PETSC_NULL) {
    ierr = DADestroy(SSADA);CHKERRQ(ierr);
  }
  
  return 0;
}


//! \brief Update the SSA solution.
PetscErrorCode SSA::update(bool fast) {
  PetscErrorCode ierr;

  if (fast)
    return 0;

  grid.profiler->begin(event_ssa);

  ierr = solve(); CHKERRQ(ierr); 

  ierr = compute_basal_frictional_heating(basal_frictional_heating); CHKERRQ(ierr);
  ierr = compute_D2(D2); CHKERRQ(ierr);

  ierr = compute_maximum_velocity(); CHKERRQ(ierr);

  grid.profiler->end(event_ssa);

  return 0;
}


//! \brief Compute the D2 term (for the strain heating computation).
/*!
  Documented in [\ref BBssasliding].
 */
PetscErrorCode SSA::compute_D2(IceModelVec2S &result) {
  PetscErrorCode ierr;
  PetscReal dx = grid.dx, dy = grid.dy;

  ierr = velocity.begin_access(); CHKERRQ(ierr);
  ierr = result.begin_access(); CHKERRQ(ierr);
  ierr = mask->begin_access(); CHKERRQ(ierr);
  
  for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
      if (mask->value(i,j) == MASK_DRAGGING_SHEET) {
        const PetscScalar 
          u_x   = (velocity(i+1,j).u - velocity(i-1,j).u)/(2*dx),
          u_y   = (velocity(i,j+1).u - velocity(i,j-1).u)/(2*dy),
          v_x   = (velocity(i+1,j).v - velocity(i-1,j).v)/(2*dx),
          v_y   = (velocity(i,j+1).v - velocity(i,j-1).v)/(2*dy),
          D2ssa = PetscSqr(u_x) + PetscSqr(v_y) + u_x * v_y
          + PetscSqr(0.5*(u_y + v_x));

        result(i,j) = D2ssa;
      } else {
        result(i,j) = 0.0;
      }
    }
  }

  ierr = mask->end_access(); CHKERRQ(ierr);
  ierr = result.end_access(); CHKERRQ(ierr);
  ierr = velocity.end_access(); CHKERRQ(ierr);
  return 0;
}


//! \brief Compute the basal frictional heating.
/*!
  Ice shelves have zero basal friction heating.
 */
PetscErrorCode SSA::compute_basal_frictional_heating(IceModelVec2S &result) {
  PetscErrorCode ierr;

  ierr = velocity.begin_access(); CHKERRQ(ierr);
  ierr = result.begin_access(); CHKERRQ(ierr);
  ierr = tauc->begin_access(); CHKERRQ(ierr);
  ierr = mask->begin_access(); CHKERRQ(ierr);
  
  for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
      if (mask->is_floating(i,j)) {
        result(i,j) = 0.0;
      } else {
        const PetscScalar 
          C = basal.drag((*tauc)(i,j), velocity(i,j).u, velocity(i,j).v),
	  basal_stress_x = - C * velocity(i,j).u,
	  basal_stress_y = - C * velocity(i,j).v;
        result(i,j) = - basal_stress_x * velocity(i,j).u - basal_stress_y * velocity(i,j).v;
      }
    }
  }

  ierr = mask->end_access(); CHKERRQ(ierr);
  ierr = tauc->end_access(); CHKERRQ(ierr);
  ierr = result.end_access(); CHKERRQ(ierr);
  ierr = velocity.end_access(); CHKERRQ(ierr);
  return 0;
}


//! \brief Computes vertically-averaged ice hardness on the staggered grid.
PetscErrorCode SSA::compute_hardav_staggered(IceModelVec2Stag &result) {
  PetscErrorCode ierr;
  PetscScalar *E, *E_ij, *E_offset;

  E = new PetscScalar[grid.Mz];

  ierr = thickness->begin_access(); CHKERRQ(ierr);
  ierr = enthalpy->begin_access(); CHKERRQ(ierr);
  ierr = result.begin_access(); CHKERRQ(ierr);
  
  for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
      ierr = enthalpy->getInternalColumn(i,j,&E_ij); CHKERRQ(ierr);
      for (PetscInt o=0; o<2; o++) {
        const PetscInt oi = 1-o, oj=o;  
        const PetscScalar H = 0.5 * ((*thickness)(i,j) + (*thickness)(i+oi,j+oj));

        if (H == 0) {
          result(i,j,o) = -1e6; // an obviously impossible value
          continue;
        }

        ierr = enthalpy->getInternalColumn(i+oi,j+oj,&E_offset); CHKERRQ(ierr);
        // build a column of enthalpy values a the current location:
        for (int k = 0; k < grid.Mz; ++k) {
          E[k] = 0.5 * (E_ij[k] + E_offset[k]);
        }
        
        result(i,j,o) = ice.averagedHardness_from_enth(H, grid.kBelowHeight(H),
                                                       grid.zlevels, E); CHKERRQ(ierr); 
      } // o
    }   // j
  }     // i

  ierr = result.end_access(); CHKERRQ(ierr);
  ierr = enthalpy->end_access(); CHKERRQ(ierr);
  ierr = thickness->end_access(); CHKERRQ(ierr);

  delete [] E;
  return 0;
}


//! \brief Compute the norm of nuH and the change in nuH.
/*!
Verification and PST experiments
suggest that an \f$L^1\f$ criterion for convergence is best.  For verification
there seems to be little difference, presumably because the solutions are smooth
and the norms are roughly equivalent on a subspace of smooth functions.  For PST,
the \f$L^1\f$ criterion gives faster runs with essentially the same results.
Presumably that is because rapid (temporal and spatial) variation in 
\f$\bar\nu H\f$ occurs at margins, occupying very few horizontal grid cells.
For the significant (e.g.~in terms of flux) parts of the flow, it is o.k. to ignore
a bit of bad behavior at these few places, and \f$L^1\f$ ignores it more than
\f$L^2\f$ (much less \f$L^\infty\f$, which might not work at all).
 */
PetscErrorCode SSA::compute_nuH_norm(PetscReal &norm, PetscReal &norm_change) {
  PetscErrorCode ierr;

  PetscReal nuNorm[2], nuChange[2];

  const PetscScalar area = grid.dx * grid.dy;
#define MY_NORM     NORM_1

  // Test for change in nu
  ierr = nuH_old.add(-1, nuH); CHKERRQ(ierr);

  ierr = nuH_old.norm_all(MY_NORM, nuChange[0], nuChange[1]); CHKERRQ(ierr);
  ierr =     nuH.norm_all(MY_NORM, nuNorm[0],   nuNorm[1]);   CHKERRQ(ierr);

  nuChange[0] *= area;
  nuChange[1] *= area;
  nuNorm[0] *= area;
  nuNorm[1] *= area;

  norm_change = sqrt(PetscSqr(nuChange[0]) + PetscSqr(nuChange[1]));
  norm = sqrt(PetscSqr(nuNorm[0]) + PetscSqr(nuNorm[1]));
  return 0;
}


//! \brief Compute the product of ice thickness and effective viscosity (on the
//! staggered grid).
/*! 
In PISM the product \f$\nu H\f$ can be
  - constant, or
  - can be computed with a constant ice hardness \f$\bar B\f$ (temperature-independent)
    but with dependence of the viscosity on the strain rates, or 
  - it can depend on the strain rates \e and have a vertically-averaged ice
    hardness depending on temperature or enthalpy.

The flow law in ice stream and ice shelf regions must, for now, be a 
(temperature-dependent) Glen law.  This is the only flow law we know how to
convert to ``viscosity form''.  More general forms like Goldsby-Kohlstedt are
not yet inverted.

The viscosity form of a Glen law is
   \f[ \nu(T^*,D) = \frac{1}{2} B(T^*) D^{(1/n)-1}\, D_{ij} \f]
where 
   \f[  D_{ij} = \frac{1}{2} \left(\frac{\partial U_i}{\partial x_j} +
                                   \frac{\partial U_j}{\partial x_i}\right) \f]
is the strain rate tensor and \f$B\f$ is an ice hardness related to 
the ice softness \f$A(T^*)\f$ by
   \f[ B(T^*)=A(T^*)^{-1/n}  \f]
in the case of a temperature dependent Glen-type law.  (Here \f$T^*\f$ is the
pressure-adjusted temperature.)

The effective viscosity is then
   \f[ \nu = \frac{\bar B}{2} \left[\left(\frac{\partial u}{\partial x}\right)^2 + 
                               \left(\frac{\partial v}{\partial y}\right)^2 + 
                               \frac{\partial u}{\partial x} \frac{\partial v}{\partial y} + 
                               \frac{1}{4} \left(\frac{\partial u}{\partial y}
                                                 + \frac{\partial v}{\partial x}\right)^2
                               \right]^{(1-n)/(2n)}  \f]
where in the temperature-dependent case
   \f[ \bar B = \frac{1}{H}\,\int_b^h B(T^*)\,dz\f]
This integral is approximately computed by the trapezoid rule.

In fact the entire effective viscosity is regularized by adding a constant.
The regularization constant \f$\epsilon\f$ is an argument to this procedure.

Also we put \f$\bar\nu H = \f$\c constantNuHForSSA anywhere the ice is thinner
than \c min_thickness_SSA.  The geometry is not changed, but this has the effect 
of producing a shelf extension in ice free ocean, which affects the driving stress
and the force balance at the calving front.
 */
PetscErrorCode SSA::compute_nuH_staggered(IceModelVec2Stag &result, PetscReal epsilon) {
  PetscErrorCode ierr;
  PISMVector2 **uv;

  ierr = result.begin_access(); CHKERRQ(ierr);
  ierr = velocity.get_array(uv); CHKERRQ(ierr);
  ierr = hardness.begin_access(); CHKERRQ(ierr);
  ierr = thickness->begin_access(); CHKERRQ(ierr);

  const PetscScalar dx = grid.dx, dy = grid.dy;

  for (PetscInt o=0; o<2; ++o) {
    const PetscInt oi = 1 - o, oj=o;
    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {

        const PetscScalar H = 0.5 * ((*thickness)(i,j) + (*thickness)(i+oi,j+oj));

        if (H < strength_extension.get_min_thickness()) {
          // Extends strength of SSA (i.e. nuH coeff) into the ice free region.
          //  Does not add or subtract ice mass.
          result(i,j,o) = strength_extension.get_notional_strength();
          continue;
        }

        PetscScalar u_x, u_y, v_x, v_y;
        // Check the offset to determine how to differentiate velocity
        if (o == 0) {
          u_x = (uv[i+1][j].u - uv[i][j].u) / dx;
          u_y = (uv[i][j+1].u + uv[i+1][j+1].u - uv[i][j-1].u - uv[i+1][j-1].u) / (4*dy);
          v_x = (uv[i+1][j].v - uv[i][j].v) / dx;
          v_y = (uv[i][j+1].v + uv[i+1][j+1].v - uv[i][j-1].v - uv[i+1][j-1].v) / (4*dy);
        } else {
          u_x = (uv[i+1][j].u + uv[i+1][j+1].u - uv[i-1][j].u - uv[i-1][j+1].u) / (4*dx);
          u_y = (uv[i][j+1].u - uv[i][j].u) / dy;
          v_x = (uv[i+1][j].v + uv[i+1][j+1].v - uv[i-1][j].v - uv[i-1][j+1].v) / (4*dx);
          v_y = (uv[i][j+1].v - uv[i][j].v) / dy;
        }

        result(i,j,o) = H * ice.effectiveViscosity(hardness(i,j,o), u_x, u_y, v_x, v_y);

        if (! finite(result(i,j,o)) || false) {
          ierr = PetscPrintf(grid.com, "nuH[%d][%d][%d] = %e\n", o, i, j, result(i,j,o));
          CHKERRQ(ierr); 
          ierr = PetscPrintf(grid.com, "  u_x, u_y, v_x, v_y = %e, %e, %e, %e\n", 
                             u_x, u_y, v_x, v_y);
          CHKERRQ(ierr);
        }
          
        // We ensure that nuH is bounded below by a positive constant.
        result(i,j,o) += epsilon;
      } // j
    } // i
  } // o

  ierr = thickness->end_access(); CHKERRQ(ierr);
  ierr = hardness.end_access(); CHKERRQ(ierr);
  ierr = result.end_access(); CHKERRQ(ierr);
  ierr = velocity.end_access(); CHKERRQ(ierr);

  // Some communication
  ierr = result.beginGhostComm(); CHKERRQ(ierr);
  ierr = result.endGhostComm(); CHKERRQ(ierr);
  return 0;
}


//! \brief Compute the driving stress.
/*!
Computes the driving stress at the base of the ice:
   \f[ \tau_d = - \rho g H \nabla h \f]

If configuration parameter \c surface_gradient_method = \c eta then the surface gradient
\f$\nabla h\f$ is computed by the gradient of the
transformed variable  \f$\eta= H^{(2n+2)/n}\f$ (frequently, \f$\eta= H^{8/3}\f$).
Because this quantity is more regular at ice sheet margins, we get a 
better surface gradient.  When the thickness at a grid point is very small
(below \c minThickEtaTransform in the procedure), the formula is slightly 
modified to give a lower driving stress.

In floating parts the surface gradient is always computed by the \c mahaffy formula.
 
Saves it in user supplied Vecs \c vtaudx and \c vtaudy, which are treated 
as global.  (I.e. we do not communicate ghosts.)
 */
PetscErrorCode SSA::compute_driving_stress(IceModelVec2V &result) {
  PetscErrorCode ierr;

  IceModelVec2S &thk = *thickness; // to improve readability (below)

  const PetscScalar n       = ice.exponent(), // frequently n = 3
                    etapow  = (2.0 * n + 2.0)/n,  // = 8/3 if n = 3
                    invpow  = 1.0 / etapow,  // = 3/8
                    dinvpow = (- n - 2.0) / (2.0 * n + 2.0); // = -5/8
  const PetscScalar minThickEtaTransform = 5.0; // m
  const PetscScalar dx=grid.dx, dy=grid.dy;

  bool compute_surf_grad_inward_ssa = config.get_flag("compute_surf_grad_inward_ssa");
  PetscReal standard_gravity = config.get("standard_gravity");
  bool use_eta = (config.get_string("surface_gradient_method") == "eta");

  ierr =   surface->begin_access();    CHKERRQ(ierr);
  ierr =       bed->begin_access();  CHKERRQ(ierr);
  ierr =      mask->begin_access();  CHKERRQ(ierr);
  ierr =        thk.begin_access();  CHKERRQ(ierr);

  ierr = result.begin_access(); CHKERRQ(ierr);

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      const PetscScalar pressure = ice.rho * standard_gravity * thk(i,j);
      if (pressure <= 0.0) {
        result(i,j).u = 0.0;
        result(i,j).v = 0.0;
      } else {
        PetscScalar h_x = 0.0, h_y = 0.0;
        // FIXME: we need to handle grid periodicity correctly.
        if (mask->is_grounded(i,j) && (use_eta == true)) {
	        // in grounded case, differentiate eta = H^{8/3} by chain rule
          if (thk(i,j) > 0.0) {
            const PetscScalar myH = (thk(i,j) < minThickEtaTransform ?
                                     minThickEtaTransform : thk(i,j));
            const PetscScalar eta = pow(myH, etapow), factor = invpow * pow(eta, dinvpow);
            h_x = factor * (pow(thk(i+1,j),etapow) - pow(thk(i-1,j),etapow)) / (2*dx);
            h_y = factor * (pow(thk(i,j+1),etapow) - pow(thk(i,j-1),etapow)) / (2*dy);
          }
          // now add bed slope to get actual h_x,h_y
          // FIXME: there is no reason to assume user's bed is periodized
          h_x += bed->diff_x(i,j);
          h_y += bed->diff_y(i,j);
        } else {  // floating or eta transformation is not used
          if (compute_surf_grad_inward_ssa) {
            h_x = surface->diff_x_p(i,j);
            h_y = surface->diff_y_p(i,j);
          } else {
            h_x = surface->diff_x(i,j);
            h_y = surface->diff_y(i,j);
          }
        }

        result(i,j).u = - pressure * h_x;
        result(i,j).v = - pressure * h_y;
      }
    }
  }

  ierr =        thk.end_access(); CHKERRQ(ierr);
  ierr =       bed->end_access(); CHKERRQ(ierr);
  ierr =   surface->end_access(); CHKERRQ(ierr);
  ierr =      mask->end_access(); CHKERRQ(ierr);
  ierr =     result.end_access(); CHKERRQ(ierr);
  return 0;
}


//! \brief Compute maximum ice velocity; uses the mask to ignore values in
//!  ice-free areas.
PetscErrorCode SSA::compute_maximum_velocity() {
  PetscErrorCode ierr;
  PetscReal my_max_u = 0.0,
    my_max_v = 0.0;

  bool do_ocean_kill = config.get_flag("ocean_kill"),
    floating_ice_killed = config.get_flag("floating_ice_killed");

  ierr = velocity.begin_access(); CHKERRQ(ierr);
  ierr = mask->begin_access(); CHKERRQ(ierr);
  
  for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
      // the following conditionals, both -ocean_kill and -float_kill, are also applied in 
      //   IceModel::massContExplicitStep() when zeroing thickness
      const bool ignorableOcean = ( do_ocean_kill && (mask->value(i,j) == MASK_OCEAN_AT_TIME_0) )
	|| ( floating_ice_killed && mask->is_floating(i,j) );
  
      if (!ignorableOcean) {
        my_max_u = PetscMax(my_max_u, PetscAbs(velocity(i,j).u));
        my_max_v = PetscMax(my_max_v, PetscAbs(velocity(i,j).v));
      }
    }
  }

  ierr = mask->end_access(); CHKERRQ(ierr);
  ierr = velocity.end_access(); CHKERRQ(ierr);

  ierr = PetscGlobalMax(&my_max_u, &max_u, grid.com); CHKERRQ(ierr); 
  ierr = PetscGlobalMax(&my_max_v, &max_v, grid.com); CHKERRQ(ierr); 
  return 0;
}


//! \brief Update the nuH viewer.
/*!
 * FIXME FIXME FIXME 
 */
PetscErrorCode SSA::update_nuH_viewers() {
  PetscErrorCode ierr;
  IceModelVec2S tmp;

  return 0;

  ierr = tmp.create(grid, "nuH", false); CHKERRQ(ierr);

  ierr = nuH.begin_access(); CHKERRQ(ierr);
  ierr = tmp.begin_access(); CHKERRQ(ierr);
  
  for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
      PetscReal avg_nuH = 0.5 * (nuH(i,j,0) + nuH(i,j,1));
        if (avg_nuH > 1.0e14) {
          tmp(i,j) = log10(avg_nuH);
        } else {
          tmp(i,j) = 14.0;
        }
    }
  }

  ierr = tmp.end_access(); CHKERRQ(ierr);
  ierr = nuH.end_access(); CHKERRQ(ierr);

  ierr = tmp.view(300); CHKERRQ(ierr);

  return 0;
}


PetscErrorCode SSA::stdout_report(string &result) {
  result = stdout_ssa;
  return 0;
}


//! \brief Set the initial guess of the SSA velocity.
PetscErrorCode SSA::set_initial_guess(IceModelVec2V &guess) {
  PetscErrorCode ierr;
  ierr = velocity.copy_from(guess); CHKERRQ(ierr);
  return 0;
}


void SSA::add_vars_to_output(string /*keyword*/, set<string> &result) {
  result.insert("velbar_ssa");
}


PetscErrorCode SSA::define_variables(set<string> vars, const NCTool &nc, nc_type nctype) {
  PetscErrorCode ierr;

  if (set_contains(vars, "velbar_ssa")) {
    ierr = velocity.define(nc, nctype); CHKERRQ(ierr);
  }

  return 0;
}


PetscErrorCode SSA::write_variables(set<string> vars, string filename) {
  PetscErrorCode ierr;

  if (set_contains(vars, "velbar_ssa")) {
    ierr = velocity.write(filename.c_str()); CHKERRQ(ierr);
  }

  return 0;
}


PetscErrorCode SSA::write_model_state(string filename) {
  PetscErrorCode ierr;

  ierr = velocity.write(filename.c_str()); CHKERRQ(ierr);

  return 0;
}

