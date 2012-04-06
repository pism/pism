// Copyright (C) 2004--2012 Constantine Khroulev, Ed Bueler and Jed Brown
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

#include "SSAFD.hh"
#include "Mask.hh"
#include "basal_resistance.hh"
#include "pism_options.hh"
#include "flowlaw_factory.hh"

SSA *SSAFDFactory(IceGrid &g, IceBasalResistancePlasticLaw &b,
                  EnthalpyConverter &ec, const NCConfigVariable &c)
{
  return new SSAFD(g,b,ec,c);
}


//! \brief Allocate objects specific to the SSAFD object.
/*!
Because the FD implementation of the SSA uses Picard iteration, a PETSc KSP
and Mat are used directly.  In particular we set up \f$A\f$
(Mat SSAStiffnessMatrix) and a \f$b\f$ (= Vec SSARHS) and iteratively solve
linear systems
  \f[ A x = b \f]
where \f$x\f$ (= Vec SSAX).  A PETSc SNES object is never created.
 */
PetscErrorCode SSAFD::allocate_fd() {
  PetscErrorCode ierr;

  // note SSADA and SSAX are allocated in SSA::allocate()
  ierr = VecDuplicate(SSAX, &SSARHS); CHKERRQ(ierr);

  ierr = DMGetMatrix(SSADA, MATAIJ, &SSAStiffnessMatrix); CHKERRQ(ierr);

  ierr = KSPCreate(grid.com, &SSAKSP); CHKERRQ(ierr);
  // the default PC type somehow is ILU, which now fails (?) while block jacobi
  //   seems to work; runtime options can override (see test J in vfnow.py)
  PC pc;
  ierr = KSPGetPC(SSAKSP,&pc); CHKERRQ(ierr);
  ierr = PCSetType(pc,PCBJACOBI); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(SSAKSP); CHKERRQ(ierr);

  const PetscScalar power = 1.0 / flow_law->exponent();
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

  scaling = 1.0e9;  // comparable to typical beta for an ice stream;

  // The nuH viewer:
  view_nuh = false;
  nuh_viewer_size = 300;
  nuh_viewer = PETSC_NULL;

  dump_system_matlab = false;

  return 0;
}

//! \brief De-allocate SIAFD internal objects.
PetscErrorCode SSAFD::deallocate_fd() {
  PetscErrorCode ierr;

  if (SSAKSP != PETSC_NULL) {
    ierr = KSPDestroy(&SSAKSP); CHKERRQ(ierr);
  }

  if (SSAStiffnessMatrix != PETSC_NULL) {
    ierr = MatDestroy(&SSAStiffnessMatrix); CHKERRQ(ierr);
  }

  if (SSARHS != PETSC_NULL) {
    ierr = VecDestroy(&SSARHS); CHKERRQ(ierr);
  }

  return 0;
}

PetscErrorCode SSAFD::init(PISMVars &vars) {
  PetscErrorCode ierr;
  ierr = SSA::init(vars); CHKERRQ(ierr);

  // The FD solver does not support direct specification of a driving stress;
  // a surface elevation must be explicitly given.
  if(surface == NULL) {
    SETERRQ(grid.com, 1, "The finite difference SSA solver requires a surface elevation.\
An explicit driving stress was specified instead and cannot be used.");
  }

  ierr = verbPrintf(2,grid.com,
                    "  [using the KSP-based finite difference implementation]\n"); CHKERRQ(ierr);

  ierr = PetscOptionsBegin(grid.com, "", "SSAFD options", ""); CHKERRQ(ierr);
  {
    bool flag;
    ierr = PISMOptionsInt("-ssa_nuh_viewer_size", "nuH viewer size",
                          nuh_viewer_size, flag); CHKERRQ(ierr);
    ierr = PISMOptionsIsSet("-ssa_view_nuh", "Enable the SSAFD nuH runtime viewer",
                            view_nuh); CHKERRQ(ierr);
  }
  ierr = PetscOptionsEnd(); CHKERRQ(ierr);

  if (config.get_flag("calving_front_stress_boundary_condition")) {
    ierr = verbPrintf(2,grid.com,
      "  using PISM-PIK calving-front stress boundary condition ...\n"); CHKERRQ(ierr);
  }

  // option to save linear system in Matlab-readable ASCII format at end of each
  // numerical solution of SSA equations; can be given with or without filename prefix
  // (i.e. "-ssa_matlab " or "-ssa_matlab foo" are both legal; in former case get
  // "pism_SSA_[year].m" if "pism_SSA" is default prefix, and in latter case get "foo_[year].m")
  string tempPrefix;
  ierr = PISMOptionsIsSet("-ssafd_matlab", "Save linear system in Matlab-readable ASCII format",
			  dump_system_matlab); CHKERRQ(ierr);

  return 0;
}

//! \brief Computes the right-hand side ("rhs") of the linear problem for the
//! Picard iteration and finite-difference implementation of the SSA equations.
/*!
The right side of the SSA equations is just the driving stress term
   \f[ - \rho g H \nabla h. \f]
The basal stress is put on the left side of the system.  This method builds the
discrete approximation of the right side.  For more about the discretization
of the SSA equations, see comments for assemble_matrix().

The values of the driving stress on the i,j grid come from a call to
compute_driving_stress().

In the case of Dirichlet boundary conditions, the entries on the right-hand side
come from known velocity values.  The fields vel_bc and bc_locations are used for
this.
 */
PetscErrorCode SSAFD::assemble_rhs(Vec rhs) {
  PetscErrorCode ierr;
  PISMVector2 **rhs_uv;
  const double dx = grid.dx, dy = grid.dy;

  Mask M;

  const double standard_gravity = config.get("standard_gravity"),
    ocean_rho = config.get("sea_water_density"),
    ice_rho = config.get("ice_density");
  const bool use_cfbc = config.get_flag("calving_front_stress_boundary_condition");

  ierr = VecSet(rhs, 0.0); CHKERRQ(ierr);

  // get driving stress components
  ierr = compute_driving_stress(taud); CHKERRQ(ierr);

  ierr = taud.begin_access(); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(SSADA, rhs, &rhs_uv); CHKERRQ(ierr);

  bool bedrock_boundary = config.get_flag("ssa_dirichlet_bc");

  if (vel_bc && bc_locations) {
    ierr = vel_bc->begin_access(); CHKERRQ(ierr);
    ierr = bc_locations->begin_access(); CHKERRQ(ierr);
  }

  if (use_cfbc) {
    ierr = thickness->begin_access(); CHKERRQ(ierr);
    ierr = bed->begin_access(); CHKERRQ(ierr);
    ierr = mask->begin_access(); CHKERRQ(ierr);
  }

  for (PetscInt i = grid.xs; i < grid.xs + grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys + grid.ym; ++j) {

      if (vel_bc && (bc_locations->as_int(i, j) == 1)) {
        rhs_uv[i][j].u = scaling * (*vel_bc)(i, j).u;
        rhs_uv[i][j].v = scaling * (*vel_bc)(i, j).v;
        continue;
      }

      if (use_cfbc) {
        PetscScalar H_ij = (*thickness)(i,j);
        PetscInt M_ij = mask->as_int(i,j),
          M_e = mask->as_int(i + 1,j),
          M_w = mask->as_int(i - 1,j),
          M_n = mask->as_int(i,j + 1),
          M_s = mask->as_int(i,j - 1);

        // Note: this sets velocities at both ice-free ocean and ice-free
        // bedrock to zero. This means that we need to set boundary conditions
        // at both ice/ice-free-ocean and ice/ice-free-bedrock interfaces below
        // to be consistent.
        if (M.ice_free(M_ij)) {
          rhs_uv[i][j].u = 0.0;
          rhs_uv[i][j].v = 0.0;
          continue;
        }

        if (is_marginal(i, j, bedrock_boundary)) {
          PetscInt aMM = 1, aPP = 1, bMM = 1, bPP = 1;
          // direct neighbors
	  if (bedrock_boundary) {
            if (M.ice_free_ocean(M_e)) aPP = 0;
            if (M.ice_free_ocean(M_w)) aMM = 0;
            if (M.ice_free_ocean(M_n)) bPP = 0;
            if (M.ice_free_ocean(M_s)) bMM = 0;}
	  else {
            if (M.ice_free(M_e)) aPP = 0;
            if (M.ice_free(M_w)) aMM = 0;
            if (M.ice_free(M_n)) bPP = 0;
            if (M.ice_free(M_s)) bMM = 0;
          }

          const double ice_pressure = ice_rho * standard_gravity * H_ij,
                       H_ij2        = H_ij*H_ij;
          double ocean_pressure,
                 h_ij = 0.0,
                 tdx  = taud(i,j).u,
                 tdy  = taud(i,j).v;

          if ((*bed)(i,j) < (sea_level - (ice_rho / ocean_rho) * H_ij)) {
            //calving front boundary condition for floating shelf
            ocean_pressure = 0.5 * ice_rho * standard_gravity * (1 - (ice_rho / ocean_rho))*H_ij2;
            // this is not really the ocean_pressure, but the difference between
            // ocean_pressure and isotrop.normal stresses (=pressure) from within
            // the ice
	    h_ij = (1.0 - ice_rho / ocean_rho) * H_ij;
						
	    // what is the force balance of an iceshelf facing a bedrock wall?! 
	    // this is not relevant as long as we ask only for ice_free_ocean neighbors
	    //if ((aPP==0 && (*bed)(i+1,j)>h_ij) || (aMM==0 && (*bed)(i-1,j)>h_ij) ||
	    //    (bPP==0 && (*bed)(i,j+1)>h_ij) || (bMM==0 && (*bed)(i,j-1)>h_ij)){
	    //  ocean_pressure = 0.0; 
	    //}

          } else {
            if( (*bed)(i,j) >= sea_level) {
              // boundary condition for a "cliff" (grounded ice next to
              // ice-free ocean) or grounded margin.
              ocean_pressure = 0.5 * ice_rho * standard_gravity * H_ij2;
              // this is not 'zero' because the isotrop.normal stresses
              // (=pressure) from within the ice figures on RHS
	      h_ij = H_ij;
            } else {
              // boundary condition for marine terminating glacier
              ocean_pressure = 0.5 * ice_rho * standard_gravity *
                (H_ij2 - (ocean_rho / ice_rho)*(sea_level - (*bed)(i,j))*(sea_level - (*bed)(i,j)));
	      h_ij = H_ij + (*bed)(i,j) - sea_level;
            }
          }

          //here we take the direct gradient at the boundary (not centered)

          if (aPP == 0 && aMM == 1) tdx = ice_pressure*h_ij / dx;
          else if (aMM == 0 && aPP == 1) tdx = -ice_pressure*h_ij / dx;
          else if (aPP == 0 && aMM == 0) tdx = 0; //in case of some kind of ice nose, or ice bridge

          if (bPP == 0 && bMM == 1) tdy = ice_pressure*h_ij / dy;
          else if (bMM == 0 && bPP == 1) tdy = -ice_pressure*h_ij / dy;
          else if (bPP == 0 && bMM == 0) tdy = 0;

          // Note that if the current cell is "marginal" but not a CFBC
          // location, the following two lines are equaivalent to the "usual
          // case" below.
          rhs_uv[i][j].u = tdx - (aMM - aPP)*ocean_pressure / dx;
          rhs_uv[i][j].v = tdy - (bMM - bPP)*ocean_pressure / dy;

          continue;
        } // end of "if (is_marginal(i, j))"

        // If we reached this point, then CFBC are enabled, but we are in the
        // interior of a sheet or shelf. See "usual case" below.

      }   // end of "if (use_cfbc)"

      // usual case: use already computed driving stress
      rhs_uv[i][j].u = taud(i, j).u;
      rhs_uv[i][j].v = taud(i, j).v;
    }
  }

  if (use_cfbc) {
    ierr = thickness->end_access(); CHKERRQ(ierr);
    ierr = bed->end_access(); CHKERRQ(ierr);
    ierr = mask->end_access(); CHKERRQ(ierr);
  }

  if (vel_bc && bc_locations) {
    ierr = bc_locations->end_access(); CHKERRQ(ierr);
    ierr = vel_bc->end_access(); CHKERRQ(ierr);
  }

  ierr = taud.end_access(); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(SSADA, rhs, &rhs_uv); CHKERRQ(ierr);

  return 0;
}


//! \brief Assemble the left-hand side matrix for the KSP-based, Picard iteration,
//! and finite difference implementation of the SSA equations.
/*!
Recall the SSA equations are
\f{align*}
 - 2 \left[\nu H \left(2 u_x + v_y\right)\right]_x
        - \left[\nu H \left(u_y + v_x\right)\right]_y
        - \tau_{(b)1}  &= - \rho g H h_x, \\
   - \left[\nu H \left(u_y + v_x\right)\right]_x
        - 2 \left[\nu H \left(u_x + 2 v_y\right)\right]_y
        - \tau_{(b)2}  &= - \rho g H h_y,
\f}
where \f$u\f$ is the \f$x\f$-component of the velocity and \f$v\f$ is the
\f$y\f$-component of the velocity.

The coefficient \f$\nu\f$ is the vertically-averaged effective viscosity.
(The product \f$\nu H\f$ is computed by compute_nuH_staggered().)
The Picard iteration idea is that, to solve the nonlinear equations in which
the effective viscosity depends on the velocity, we freeze the effective
viscosity using its value at the current estimate of the velocity and we solve
the linear equations which come from this viscosity.  In abstract symbols, the
Picard iteration replaces the above nonlinear SSA equations by a sequence of
linear problems
	\f[ A(U^{(k)}) U^{(k+1)} = b \f]
where \f$A(U)\f$ is the matrix from discretizing the SSA equations supposing
the viscosity is a function of known velocities \f$U\f$, and where \f$U^{(k)}\f$
denotes the \f$k\f$th iterate in the process.  The current method assembles \f$A(U)\f$.

For ice shelves \f$\tau_{(b)i} = 0\f$ [\ref MacAyealetal].
For ice streams with a basal till modelled as a plastic material,
\f$\tau_{(b)i} = - \tau_c u_i/|\mathbf{u}|\f$ where
\f$\mathbf{u} = (u,v)\f$, \f$|\mathbf{u}| = \left(u^2 + v^2\right)^{1/2}\f$, where
\f$\tau_c(t,x,y)\f$ is the yield stress of the till [\ref SchoofStream].
More generally, ice streams can be modeled with a pseudo-plastic basal till;
see IceModel::initBasalTillModel() and IceModel::updateYieldStressUsingBasalWater()
and reference [\ref BKAJS].  The pseudo-plastic till model includes all power law
sliding relations and the linearly-viscous model for sliding,
\f$\tau_{(b)i} = - \beta u_i\f$ where \f$\beta\f$ is the basal drag
(friction) parameter [\ref MacAyeal].  In any case, PISM assumes that the basal shear
stress can be factored this way, <i>even if the coefficient depends on the
velocity</i>, \f$\beta(u,v)\f$.  Such factoring is possible even in the case of
(regularized) plastic till.  This scalar coefficient \f$\beta\f$ is what is
returned by IceBasalResistancePlasticLaw::drag().

Note that the basal shear stress appears on the \em left side of the linear
system we actually solve.  We believe this is crucial, because of its effect on
the spectrum of the linear approximations of each stage.  The effect on spectrum
is clearest in the linearly-viscous till case but
there seems to be an analogous effect in the plastic till case.

This method assembles the matrix for the left side of the above SSA equations.
The numerical method is finite difference.  Suppose we use difference notation
\f$\delta_{+x}f^{i,j} = f^{i+1,j}-f^{i,j}\f$,
\f$\delta_{-x}f^{i,j} = f^{i,j}-f^{i-1,j}\f$, and
\f$\Delta_{x}f^{i,j} = f^{i+1,j}-f^{i-1,j}\f$, and corresponding notation for
\f$y\f$ differences, and that we write \f$N = \nu H\f$ then the first of the
two "concrete" SSA equations above has this discretization:
\f{align*}
- &2 \frac{N^{i+\frac{1}{2},j}}{\Delta x} \left[2\frac{\delta_{+x}u^{i,j}}{\Delta x} + \frac{\Delta_{y} v^{i+1,j} + \Delta_{y} v^{i,j}}{4 \Delta y}\right] + 2 \frac{N^{i-\frac{1}{2},j}}{\Delta x} \left[2\frac{\delta_{-x}u^{i,j}}{\Delta x} + \frac{\Delta_y v^{i,j} + \Delta_y v^{i-1,j}}{4 \Delta y}\right] \\
&\qquad- \frac{N^{i,j+\frac{1}{2}}}{\Delta y} \left[\frac{\delta_{+y} u^{i,j}}{\Delta y} + \frac{\Delta_x v^{i,j+1} + \Delta_x v^{i,j}}{4 \Delta x}\right] + \frac{N^{i,j-\frac{1}{2}}}{\Delta y} \left[\frac{\delta_{-y}u^{i,j}}{\Delta y} + \frac{\Delta_x v^{i,j} + \Delta_x v^{i,j-1}}{4 \Delta x}\right] - \tau_{(b)1}^{i,j} = - \rho g H^{i,j} \frac{\Delta_x h^{i,j}}{2\Delta x}.
\f}
As a picture, see Figure \ref ssastencil.

\image html ssastencil.png "\b ssastencil:  Stencil for our finite difference discretization of the first of the two scalar SSA equations.  Triangles show staggered grid points where N = nu * H is evaluated.  Circles and squares show where u and v are approximated, respectively."
\anchor ssastencil

It follows immediately that the matrix we assemble in the current method has
13 nonzeros entries per row because, for this first SSA equation, there are 5
grid values of \f$u\f$ and 8 grid values of \f$v\f$ used in this scheme.  For
the second equation we also have 13 nonzeros per row.

FIXME:  document use of DAGetMatrix and MatStencil and MatSetValuesStencil

*/
PetscErrorCode SSAFD::assemble_matrix(bool include_basal_shear, Mat A) {
  PetscErrorCode  ierr;

  const PetscScalar   dx=grid.dx, dy=grid.dy;
  const PetscScalar   beta_ice_free_bedrock = config.get("beta_ice_free_bedrock");
  const bool use_cfbc = config.get_flag("calving_front_stress_boundary_condition");

  // shortcut:
  IceModelVec2V &vel = velocity;

  ierr = MatZeroEntries(A); CHKERRQ(ierr);

  /* matrix assembly loop */

  ierr = nuH.begin_access(); CHKERRQ(ierr);
  ierr = tauc->begin_access(); CHKERRQ(ierr);
  ierr = vel.begin_access(); CHKERRQ(ierr);
  ierr = mask->begin_access(); CHKERRQ(ierr);

  Mask M;

  const bool bedrock_boundary = config.get_flag("ssa_dirichlet_bc");

  if (vel_bc && bc_locations) {
    ierr = bc_locations->begin_access(); CHKERRQ(ierr);
  }

  // handles friction of the ice cell along ice-free bedrock margins when bedrock higher than ice surface (in simplified setups)
  bool nuBedrockSet=config.get_flag("nuBedrockSet");
  if (nuBedrockSet) {
    ierr =    thickness->begin_access();  CHKERRQ(ierr);
    ierr =    bed->begin_access();        CHKERRQ(ierr);
    ierr =    surface->begin_access();    CHKERRQ(ierr);
  }
  PetscScalar nuBedrock=config.get("nuBedrock");
  PetscScalar HminFrozen=0.0;

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {

      // Handle the easy case: provided Dirichlet boundary conditions
      if (vel_bc && bc_locations && bc_locations->as_int(i,j) == 1) {
        // set diagonal entry to one (scaled); RHS entry will be known velocity;
        ierr = set_diagonal_matrix_entry(A, i, j, scaling); CHKERRQ(ierr);
        continue;
      }

      /* Provide shorthand for the following staggered coefficients  nu H:
       *      c_n
       *  c_w     c_e
       *      c_s
       */
      // const
      PetscScalar c_w = nuH(i-1,j,0);
      PetscScalar c_e = nuH(i,j,0);
      PetscScalar c_s = nuH(i,j-1,1);
      PetscScalar c_n = nuH(i,j,1);

      if (nuBedrockSet){
       // if option is set, the viscosity at ice-bedrock boundary layer will
       // be prescribed and is a temperature-independent free (user determined) parameter

	// direct neighbors
	PetscInt  M_e = mask->as_int(i + 1,j),
	          M_w = mask->as_int(i - 1,j),
	          M_n = mask->as_int(i,j + 1),
		  M_s = mask->as_int(i,j - 1);

        if ((*thickness)(i,j) > HminFrozen) {  
	  if ((*bed)(i-1,j) > (*surface)(i,j) && M.ice_free_land(M_w)) {
	    c_w = nuBedrock * 0.5 * ((*thickness)(i,j)+(*thickness)(i-1,j));	    
	  }
	  if ((*bed)(i+1,j) > (*surface)(i,j) && M.ice_free_land(M_e)) {
	   c_e = nuBedrock * 0.5 * ((*thickness)(i,j)+(*thickness)(i+1,j));
	  }
	  if ((*bed)(i,j+1) > (*surface)(i,j) && M.ice_free_land(M_n)) {
	    c_n = nuBedrock * 0.5 * ((*thickness)(i,j)+(*thickness)(i,j+1));
  	  }
	  if ((*bed)(i,j-1) > (*surface)(i,j) && M.ice_free_land(M_s)) {
	    c_s = nuBedrock * 0.5 * ((*thickness)(i,j)+(*thickness)(i+1,j));
	  }
        }
      }

      // We use DAGetMatrix to obtain the SSA matrix, which means that all 18
      // non-zeros get allocated, even though we use only 13 (or 14). The
      // remaining 5 (or 4) coefficients are zeros, but we set them anyway,
      // because this makes the code easier to understand.
      const PetscInt sten = 18;
      MatStencil row, col[sten];

      PetscInt aMn = 1, aPn = 1, aMM = 1, aPP = 1, aMs = 1, aPs = 1;
      PetscInt bPw = 1, bPP = 1, bPe = 1, bMw = 1, bMM = 1, bMe = 1;

      PetscInt M_ij = mask->as_int(i,j);

      if (use_cfbc) {
        int
          // direct neighbors
          M_e = mask->as_int(i + 1,j),
          M_w = mask->as_int(i - 1,j),
          M_n = mask->as_int(i,j + 1),
          M_s = mask->as_int(i,j - 1),
          // "diagonal" neighbors
          M_ne = mask->as_int(i + 1,j + 1),
          M_se = mask->as_int(i + 1,j - 1),
          M_nw = mask->as_int(i - 1,j + 1),
          M_sw = mask->as_int(i - 1,j - 1);

        // Note: this sets velocities at both ice-free ocean and ice-free
        // bedrock to zero. This means that we need to set boundary conditions
        // at both ice/ice-free-ocean and ice/ice-free-bedrock interfaces below
        // to be consistent.
        if (M.ice_free(M_ij)) {
          ierr = set_diagonal_matrix_entry(A, i, j, scaling); CHKERRQ(ierr);
          continue;
        }

        if (is_marginal(i, j, bedrock_boundary)) {
          // If at least one of the following four conditions is "true", we're
          // at a CFBC location.
	  if (bedrock_boundary) {

            if (M.ice_free_ocean(M_e)) aPP = 0;
            if (M.ice_free_ocean(M_w)) aMM = 0;
            if (M.ice_free_ocean(M_n)) bPP = 0;
            if (M.ice_free_ocean(M_s)) bMM = 0;

            // decide whether to use centered or one-sided differences
            if (M.ice_free_ocean(M_n) || M.ice_free_ocean(M_ne)) aPn = 0;
            if (M.ice_free_ocean(M_e) || M.ice_free_ocean(M_ne)) bPe = 0;
            if (M.ice_free_ocean(M_e) || M.ice_free_ocean(M_se)) bMe = 0;
            if (M.ice_free_ocean(M_s) || M.ice_free_ocean(M_se)) aPs = 0;
            if (M.ice_free_ocean(M_s) || M.ice_free_ocean(M_sw)) aMs = 0;
            if (M.ice_free_ocean(M_w) || M.ice_free_ocean(M_sw)) bMw = 0;
            if (M.ice_free_ocean(M_w) || M.ice_free_ocean(M_nw)) bPw = 0;
            if (M.ice_free_ocean(M_n) || M.ice_free_ocean(M_nw)) aMn = 0;}

	  else {

            if (M.ice_free(M_e)) aPP = 0;
            if (M.ice_free(M_w)) aMM = 0;
            if (M.ice_free(M_n)) bPP = 0;
            if (M.ice_free(M_s)) bMM = 0;

            // decide whether to use centered or one-sided differences
            if (M.ice_free(M_n) || M.ice_free(M_ne)) aPn = 0;
            if (M.ice_free(M_e) || M.ice_free(M_ne)) bPe = 0;
            if (M.ice_free(M_e) || M.ice_free(M_se)) bMe = 0;
            if (M.ice_free(M_s) || M.ice_free(M_se)) aPs = 0;
            if (M.ice_free(M_s) || M.ice_free(M_sw)) aMs = 0;
            if (M.ice_free(M_w) || M.ice_free(M_sw)) bMw = 0;
            if (M.ice_free(M_w) || M.ice_free(M_nw)) bPw = 0;
            if (M.ice_free(M_n) || M.ice_free(M_nw)) aMn = 0;				}
           }
      } // end of "if (use_cfbc)"

      /* begin Maxima-generated code */
      const PetscReal dx2 = dx*dx, dy2 = dy*dy, d4 = 4*dx*dy, d2 = 2*dx*dy;

      /* Coefficients of the discretization of the first equation; u first, then v. */
      PetscReal eq1[] = {
        0,  -c_n*bPP/dy2,  0,
        -4*c_w*aMM/dx2,  (c_n*bPP+c_s*bMM)/dy2+(4*c_e*aPP+4*c_w*aMM)/dx2,  -4*c_e*aPP/dx2,
        0,  -c_s*bMM/dy2,  0,
        c_w*aMM*bPw/d2+c_n*aMn*bPP/d4,  (c_n*aPn*bPP-c_n*aMn*bPP)/d4+(c_w*aMM*bPP-c_e*aPP*bPP)/d2,  -c_e*aPP*bPe/d2-c_n*aPn*bPP/d4,
        (c_w*aMM*bMw-c_w*aMM*bPw)/d2+(c_n*aMM*bPP-c_s*aMM*bMM)/d4,  (c_n*aPP*bPP-c_n*aMM*bPP-c_s*aPP*bMM+c_s*aMM*bMM)/d4+(c_e*aPP*bPP-c_w*aMM*bPP-c_e*aPP*bMM+c_w*aMM*bMM)/d2,  (c_e*aPP*bPe-c_e*aPP*bMe)/d2+(c_s*aPP*bMM-c_n*aPP*bPP)/d4,
        -c_w*aMM*bMw/d2-c_s*aMs*bMM/d4,  (c_s*aMs*bMM-c_s*aPs*bMM)/d4+(c_e*aPP*bMM-c_w*aMM*bMM)/d2,  c_e*aPP*bMe/d2+c_s*aPs*bMM/d4,
      };

      /* Coefficients of the discretization of the second equation; u first, then v. */
      PetscReal eq2[] = {
        c_w*aMM*bPw/d4+c_n*aMn*bPP/d2,  (c_n*aPn*bPP-c_n*aMn*bPP)/d2+(c_w*aMM*bPP-c_e*aPP*bPP)/d4,  -c_e*aPP*bPe/d4-c_n*aPn*bPP/d2,
        (c_w*aMM*bMw-c_w*aMM*bPw)/d4+(c_n*aMM*bPP-c_s*aMM*bMM)/d2,  (c_n*aPP*bPP-c_n*aMM*bPP-c_s*aPP*bMM+c_s*aMM*bMM)/d2+(c_e*aPP*bPP-c_w*aMM*bPP-c_e*aPP*bMM+c_w*aMM*bMM)/d4,  (c_e*aPP*bPe-c_e*aPP*bMe)/d4+(c_s*aPP*bMM-c_n*aPP*bPP)/d2,
        -c_w*aMM*bMw/d4-c_s*aMs*bMM/d2,  (c_s*aMs*bMM-c_s*aPs*bMM)/d2+(c_e*aPP*bMM-c_w*aMM*bMM)/d4,  c_e*aPP*bMe/d4+c_s*aPs*bMM/d2,
        0,  -4*c_n*bPP/dy2,  0,
        -c_w*aMM/dx2,  (4*c_n*bPP+4*c_s*bMM)/dy2+(c_e*aPP+c_w*aMM)/dx2,  -c_e*aPP/dx2,
        0,  -4*c_s*bMM/dy2,  0,
      };

      /* i indices */
      const PetscInt I[] = {
        i-1,  i,  i+1,
        i-1,  i,  i+1,
        i-1,  i,  i+1,
        i-1,  i,  i+1,
        i-1,  i,  i+1,
        i-1,  i,  i+1,
      };

      /* j indices */
      const PetscInt J[] = {
        j+1,  j+1,  j+1,
        j,  j,  j,
        j-1,  j-1,  j-1,
        j+1,  j+1,  j+1,
        j,  j,  j,
        j-1,  j-1,  j-1,
      };

      /* component indices */
      const PetscInt C[] = {
        0,  0,  0,
        0,  0,  0,
        0,  0,  0,
        1,  1,  1,
        1,  1,  1,
        1,  1,  1,
      };
      /* end Maxima-generated code */

      /* Dragging ice experiences friction at the bed determined by the
       *    IceBasalResistancePlasticLaw::drag() methods.  These may be a plastic,
       *    pseudo-plastic, or linear friction law.  Dragging is done implicitly
       *    (i.e. on left side of SSA eqns).  */
      PetscReal beta = 0.0;
      if (include_basal_shear) {
        if (M.grounded_ice(M_ij)) {
          beta = basal.drag((*tauc)(i,j), vel(i,j).u, vel(i,j).v);
        } else if (M.ice_free_land(M_ij)) {
          // apply drag even in this case, to help with margins; note ice free
          // areas already have a strength extension
          beta = beta_ice_free_bedrock;
        }
      }

      // add beta to diagonal entries
      eq1[4]  += beta;
      eq2[13] += beta;

      // build equations: NOTE TRANSPOSE
      row.j = i; row.i = j;
      for (PetscInt m = 0; m < sten; m++) {
        col[m].j = I[m]; col[m].i = J[m]; col[m].c = C[m];
      }

      // set coefficients of the first equation:
      row.c = 0;
      ierr = MatSetValuesStencil(A, 1, &row, sten, col, eq1, INSERT_VALUES); CHKERRQ(ierr);

      // set coefficients of the second equation:
      row.c = 1;
      ierr = MatSetValuesStencil(A, 1, &row, sten, col, eq2, INSERT_VALUES); CHKERRQ(ierr);
    }
  }

  if (vel_bc && bc_locations) {
    ierr = bc_locations->end_access(); CHKERRQ(ierr);
  }

  ierr = mask->end_access(); CHKERRQ(ierr);
  ierr = vel.end_access(); CHKERRQ(ierr);
  ierr = tauc->end_access(); CHKERRQ(ierr);
  ierr = nuH.end_access(); CHKERRQ(ierr);

  if (nuBedrockSet) {
  	ierr =    thickness->end_access();    CHKERRQ(ierr);
  	ierr =  		bed->end_access();  CHKERRQ(ierr);
  	ierr =    	surface->end_access();    CHKERRQ(ierr);
  }

  ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
#if (PISM_DEBUG==1)
  ierr = MatSetOption(A,MAT_NEW_NONZERO_LOCATION_ERR,PETSC_TRUE);CHKERRQ(ierr);
#endif

  return 0;
}


//! \brief Compute the vertically-averaged horizontal velocity from the shallow
//! shelf approximation.
/*!
This is the main procedure in the SSAFD.  It manages the nonlinear solve process
and the Picard iteration.

The outer loop (over index \c k) is the nonlinear iteration.  In this loop the effective
viscosity is computed by compute_nuH_staggered() and then the linear system is
set up and solved.

Specifically, we call the PETSc procedure KSPSolve() to solve the linear system.
Solving the linear system is also a loop, an iteration, but it occurs
inside KSPSolve().  The user has full control of the KSP solve through the PETSc
interface.  The default choicess for KSP type <tt>-ksp_type</tt> and preconditioner type
<tt>-pc_type</tt> are GMRES(30) for the former and block Jacobi with ILU on the
blocks for the latter.  The defaults usually work.  These choices are important
but poorly understood.  The eigenvalues of the linearized
SSA vary with ice sheet geometry and temperature in ways that are not well-studied.
Nonetheless these eigenvalues determine the convergence of
this (inner) linear iteration.  A well-chosen preconditioner can put the eigenvalues
in the right place so that the KSP can converge quickly.

The preconditioner will behave differently on different numbers of
processors.  If the user wants the results of SSA calculations to be
independent of the number of processors, then <tt>-pc_type none</tt> could
be used, but performance will be poor.

If you want to test different KSP methods, it may be helpful to see how many
iterations were necessary.  Use <tt>-ksp_monitor</tt>.
Initial testing implies that CGS takes roughly half the iterations of
GMRES(30), but is not significantly faster because the iterations are each
roughly twice as slow.  The outputs of PETSc options <tt>-ksp_monitor_singular_value</tt>,
<tt>-ksp_compute_eigenvalues</tt> and <tt>-ksp_plot_eigenvalues -draw_pause N</tt>
(the last holds plots for N seconds) may be useful to diagnose.

The outer loop terminates when the effective viscosity times thickness is no longer changing
much, according to the tolerance set by the option <tt>-ssa_rtol</tt>.  The
outer loop also terminates when a maximum number of iterations is exceeded.
We save the velocity from the last time step in order to have a better estimate
of the effective viscosity than the u=v=0 result.

In truth there is an "outer outer" loop (over index \c l).  This attempts to
over-regularize the effective viscosity if the nonlinear iteration (the "outer"
loop over \c k) is not converging with the default regularization.  The same
over-regularization is attempted if the KSP object reports that it has not
converged.

(An alternative for recovery in the KSP diverged case, suggested by Jed, is to
revert to a direct linear solve, either for the whole domain (not scalable) or
on the subdomains.  This recovery alternative requires a more nontrivial choice
but it may be worthwhile.  Note the user can already do <tt>-pc_type asm
-sub_pc_type lu</tt> at the command line, forcing subdomain direct solves.)
 */
PetscErrorCode SSAFD::solve() {
  PetscErrorCode ierr;
  Mat A = SSAStiffnessMatrix; // solve  A SSAX = SSARHS
  PetscReal   norm, normChange;
  PetscInt    ksp_iterations, ksp_iterations_total = 0, outer_iterations;
  KSPConvergedReason  reason;

  stdout_ssa.clear();

  PetscReal ssaRelativeTolerance = config.get("ssafd_relative_convergence"),
            epsilon              = config.get("epsilon_ssa");
  PetscInt ssaMaxIterations = static_cast<PetscInt>(config.get("max_iterations_ssafd"));
  // this has no units; epsilon goes up by this ratio when previous value failed
  const PetscScalar DEFAULT_EPSILON_MULTIPLIER_SSA = 4.0;

  ierr = velocity.copy_to(velocity_old); CHKERRQ(ierr);

  // computation of RHS only needs to be done once; does not depend on
  // solution; but matrix changes under nonlinear iteration (loop over k below)
  ierr = assemble_rhs(SSARHS); CHKERRQ(ierr);

  ierr = compute_hardav_staggered(hardness); CHKERRQ(ierr);

  for (PetscInt l=0; ; ++l) { // iterate with increasing regularization parameter
    ierr = compute_nuH_staggered(nuH, epsilon); CHKERRQ(ierr);

    ierr = update_nuH_viewers(); CHKERRQ(ierr);
    // iterate on effective viscosity: "outer nonlinear iteration":
    for (PetscInt k = 0; k < ssaMaxIterations; ++k) {
      if (getVerbosityLevel() > 2) {
        char tempstr[50] = "";  snprintf(tempstr,50, "  %d,%2d:", l, k);
        stdout_ssa += tempstr;
      }

      // in preparation of measuring change of effective viscosity:
      ierr = nuH.copy_to(nuH_old); CHKERRQ(ierr);

      do {
        // assemble (or re-assemble) matrix, which depends on updated viscosity
        ierr = assemble_matrix(true, A); CHKERRQ(ierr);
        if (getVerbosityLevel() > 2)
          stdout_ssa += "A:";

        // call PETSc to solve linear system by iterative method; "inner iteration"
        ierr = KSPSetOperators(SSAKSP, A, A, SAME_NONZERO_PATTERN); CHKERRQ(ierr);
        ierr = KSPSolve(SSAKSP, SSARHS, SSAX); CHKERRQ(ierr); // SOLVE

        // check if diverged; report to standard out about iteration
        ierr = KSPGetConvergedReason(SSAKSP, &reason); CHKERRQ(ierr);
        if (reason < 0) {
          // KSP diverged
          ierr = verbPrintf(1,grid.com,
              "\nPISM WARNING:  KSPSolve() reports 'diverged'; reason = %d = '%s'\n",
              reason,KSPConvergedReasons[reason]); CHKERRQ(ierr);
          // write a file with a fixed filename; avoid zillions of files
          char filename[PETSC_MAX_PATH_LEN] = "SSAFD_kspdivergederror.petsc";
          ierr = verbPrintf(1,grid.com,
              "  writing linear system to PETSc binary file %s ...\n",filename);
              CHKERRQ(ierr);
          PetscViewer    viewer;
          ierr = PetscViewerBinaryOpen(grid.com, filename, FILE_MODE_WRITE, 
                                       &viewer); CHKERRQ(ierr);
          ierr = MatView(A,viewer); CHKERRQ(ierr);
          ierr = VecView(SSARHS,viewer); CHKERRQ(ierr);
          ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
          // attempt recovery by same mechanism as for outer iteration
          //   (FIXME: could force direct solve on subdomains?)
          if (epsilon <= 0.0) {
            ierr = verbPrintf(1,grid.com,
                "KSP diverged AND epsilon<=0.0.\n  ENDING ...\n"); CHKERRQ(ierr);
            PISMEnd();
          }
          ierr = verbPrintf(1,grid.com,
              "  KSP diverged with epsilon=%8.2e.  Retrying with epsilon multiplied by %8.2e.\n\n",
            epsilon, DEFAULT_EPSILON_MULTIPLIER_SSA); CHKERRQ(ierr);
          epsilon *= DEFAULT_EPSILON_MULTIPLIER_SSA;
          // recovery requires recompute on nuH  (FIXME: could be implemented by max(eps,nuH) here?)
          ierr = compute_nuH_staggered(nuH, epsilon); CHKERRQ(ierr);
        } else {
          // report on KSP success; the "inner" iteration is done
          ierr = KSPGetIterationNumber(SSAKSP, &ksp_iterations); CHKERRQ(ierr);
          ksp_iterations_total += ksp_iterations;
          if (getVerbosityLevel() > 2) {
            char tempstr[50] = "";  snprintf(tempstr,50, "S:%d,%d: ", ksp_iterations, reason);
            stdout_ssa += tempstr;
          }
        }

      } while (reason < 0);  // keep trying till KSP converged

      // Communicate so that we have stencil width for evaluation of effective
      //   viscosity on next "outer" iteration (and geometry etc. if done):
      ierr = velocity.copy_from(SSAX); CHKERRQ(ierr);

      ierr = velocity.beginGhostComm(); CHKERRQ(ierr);
      ierr = velocity.endGhostComm(); CHKERRQ(ierr);

      // update viscosity and check for viscosity convergence
      ierr = compute_nuH_staggered(nuH, epsilon); CHKERRQ(ierr);
      ierr = update_nuH_viewers(); CHKERRQ(ierr);
      ierr = compute_nuH_norm(norm, normChange); CHKERRQ(ierr);
      if (getVerbosityLevel() > 2) {
        char tempstr[100] = "";
        snprintf(tempstr,100, "|nu|_2, |Delta nu|_2/|nu|_2 = %10.3e %10.3e\n",
                         norm, normChange/norm);
        stdout_ssa += tempstr;
      }

      if (getVerbosityLevel() > 2) { // assume that high verbosity shows interest
                                     //   in immediate feedback about SSA iterations
        ierr = verbPrintf(2,grid.com, stdout_ssa.c_str()); CHKERRQ(ierr);
        stdout_ssa.clear();
      }

      outer_iterations = k + 1;
      if (norm == 0 || normChange / norm < ssaRelativeTolerance) goto done;

    } // end of the "outer loop" (index: k)

    if (epsilon > 0.0) {
       ierr = verbPrintf(1,grid.com,
	 "WARNING: Effective viscosity not converged after %d outer iterations\n"
	 "\twith epsilon=%8.2e. Retrying with epsilon multiplied by %8.2e.\n",
	 ssaMaxIterations, epsilon, DEFAULT_EPSILON_MULTIPLIER_SSA); CHKERRQ(ierr);
       ierr = velocity.copy_from(velocity_old); CHKERRQ(ierr);
       epsilon *= DEFAULT_EPSILON_MULTIPLIER_SSA;
    } else {
       ierr = verbPrintf(1,grid.com,
         "Effective viscosity not converged after %d iterations AND epsilon<=0.0.\n"
         "  ENDING ...\n", ssaMaxIterations); CHKERRQ(ierr);
       PISMEnd();
    }

  } // end of the "outer outer loop" (index: l)

  done:

  if (getVerbosityLevel() > 2) {
    char tempstr[100] = "";
    snprintf(tempstr, 100, "... =%5d outer iterations, ~%3.1f KSP iterations each\n",
             outer_iterations, ((double) ksp_iterations_total) / outer_iterations);
    stdout_ssa += tempstr;
  } else if (getVerbosityLevel() == 2) {
    // at default verbosity, just record last normchange and iterations
    char tempstr[100] = "";
    snprintf(tempstr, 100, "%5d outer iterations, ~%3.1f KSP iterations each\n",
             outer_iterations, ((double) ksp_iterations_total) / outer_iterations);
    stdout_ssa += tempstr;
  }
  if (getVerbosityLevel() >= 2)
    stdout_ssa = "  SSA: " + stdout_ssa;

  if (config.get_flag("write_ssa_system_to_matlab")) {
    ierr = writeSSAsystemMatlab(); CHKERRQ(ierr);
  }

  if (config.get_flag("scalebrutalSet")){
    const PetscScalar sliding_scale_brutalFactor = config.get("sliding_scale_brutal");
    ierr = velocity.scale(sliding_scale_brutalFactor); CHKERRQ(ierr);

    ierr = velocity.beginGhostComm(); CHKERRQ(ierr);
    ierr = velocity.endGhostComm(); CHKERRQ(ierr);
  }

  return 0;
}


//! \brief Write the SSA system to an .m (MATLAB) file (for debugging).
PetscErrorCode SSAFD::writeSSAsystemMatlab() {
  PetscErrorCode ierr;
  PetscViewer    viewer;
  string prefix = "pism_SSAFD", file_name;
  char           yearappend[PETSC_MAX_PATH_LEN];

  IceModelVec2S component;
  ierr = component.create(grid, "temp_storage", false); CHKERRQ(ierr);

  bool flag;
  ierr = PISMOptionsString("ssafd_matlab",
                           "Save the linear system to an ASCII .m file. Sets the file prefix.",
                           prefix, flag); CHKERRQ(ierr);

  snprintf(yearappend, PETSC_MAX_PATH_LEN, "_y%.0f.m", grid.time->seconds_to_years(grid.time->current()));
  file_name = prefix + string(yearappend);

  ierr = verbPrintf(2, grid.com,
                    "writing Matlab-readable file for SSAFD system A xsoln = rhs to file `%s' ...\n",
                    file_name.c_str()); CHKERRQ(ierr);
  ierr = PetscViewerCreate(grid.com, &viewer);CHKERRQ(ierr);
  ierr = PetscViewerSetType(viewer, PETSCVIEWERASCII);CHKERRQ(ierr);
  ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);CHKERRQ(ierr);
  ierr = PetscViewerFileSetName(viewer, file_name.c_str());CHKERRQ(ierr);

  // get the command which started the run
  string cmdstr = pism_args_string();

  // save linear system; gives system A xsoln = rhs at last (nonlinear) iteration of SSA
  ierr = PetscViewerASCIIPrintf(viewer,
    "%% A PISM linear system report for the finite difference implementation\n"
    "%% of the SSA stress balance from this run:\n"
    "%%   '%s'\n"
    "%% Writes matrix A (sparse), and vectors uv and rhs, for the linear\n"
    "%% system which was solved at the last step of the nonlinear iteration:\n"
    "%%    A * uv = rhs.\n"
    "%% Also writes the year, the coordinates x,y, their gridded versions\n"
    "%% xx,yy, and the thickness (thk) and surface elevation (usurf).\n"
    "%% Also writes i-offsetvalues of vertically-integrated viscosity\n"
    "%% (nuH_0 = nu * H), and j-offset version of same thing (nuH_1 = nu * H);\n"
    "%% these are on the staggered grid.\n",
    cmdstr.c_str());  CHKERRQ(ierr);

  ierr = PetscViewerASCIIPrintf(viewer,"\n\necho off\n");  CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) SSAStiffnessMatrix,"A"); CHKERRQ(ierr);
  ierr = MatView(SSAStiffnessMatrix, viewer);CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"clear zzz\n\n");  CHKERRQ(ierr);

  ierr = PetscObjectSetName((PetscObject) SSARHS,"rhs"); CHKERRQ(ierr);
  ierr = VecView(SSARHS, viewer);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) SSAX,"uv"); CHKERRQ(ierr);
  ierr = VecView(SSAX, viewer);CHKERRQ(ierr);

  // save coordinates (for viewing, primarily)
  ierr = PetscViewerASCIIPrintf(viewer,"\nyear=%10.6f;\n",grid.time->seconds_to_years(grid.time->current()));  CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
            "x=%12.3f + (0:%d)*%12.3f;\n"
            "y=%12.3f + (0:%d)*%12.3f;\n",
            -grid.Lx,grid.Mx-1,grid.dx,-grid.Ly,grid.My-1,grid.dy); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"[xx,yy]=meshgrid(x,y);\n");  CHKERRQ(ierr);

  // also save thickness and effective viscosity
  ierr = thickness->view_matlab(viewer); CHKERRQ(ierr);
  ierr = surface->view_matlab(viewer); CHKERRQ(ierr);

  ierr = nuH.get_component(0, component); CHKERRQ(ierr);
  ierr = component.set_name("nuH_0"); CHKERRQ(ierr);
  ierr = component.set_attr("long_name",
    "effective viscosity times thickness (i offset) at current time step"); CHKERRQ(ierr);
  ierr = component.view_matlab(viewer); CHKERRQ(ierr);

  ierr = nuH.get_component(0, component); CHKERRQ(ierr);
  ierr = component.set_name("nuH_1"); CHKERRQ(ierr);
  ierr = component.set_attr("long_name",
    "effective viscosity times thickness (j offset) at current time step"); CHKERRQ(ierr);
  ierr = component.view_matlab(viewer); CHKERRQ(ierr);

  ierr = PetscViewerASCIIPrintf(viewer,"echo on\n");  CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);

  return 0;
}

//! \brief Compute the norm of nu H and the change in nu H.
/*!
Verification and PST experiments
suggest that an \f$L^1\f$ criterion for convergence is best.  For verification
there seems to be little difference, presumably because the solutions are smooth
and the norms are roughly equivalent on a subspace of smooth functions.  For PST,
the \f$L^1\f$ criterion gives faster runs with essentially the same results.
Presumably that is because rapid (temporal and spatial) variation in
\f$\nu H\f$ occurs at margins, occupying very few horizontal grid cells.
For the significant (e.g.~in terms of flux) parts of the flow, it is o.k. to ignore
a bit of bad behavior at these few places, and \f$L^1\f$ ignores it more than
\f$L^2\f$ (much less \f$L^\infty\f$, which might not work at all).
 */
PetscErrorCode SSAFD::compute_nuH_norm(PetscReal &norm, PetscReal &norm_change) {
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

//! \brief Computes vertically-averaged ice hardness on the staggered grid.
PetscErrorCode SSAFD::compute_hardav_staggered(IceModelVec2Stag &result) {
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

        result(i,j,o) = flow_law->averaged_hardness(H, grid.kBelowHeight(H),
                                                       &grid.zlevels[0], E); CHKERRQ(ierr);
      } // o
    }   // j
  }     // i

  ierr = result.end_access(); CHKERRQ(ierr);
  ierr = enthalpy->end_access(); CHKERRQ(ierr);
  ierr = thickness->end_access(); CHKERRQ(ierr);

  delete [] E;
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
invert to ``viscosity form''.  More general forms like Goldsby-Kohlstedt are
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

The result of this procedure is \f$\nu H\f$, not just \f$\nu\f$, this it is
a vertical integral, not average, of viscosity.

The resulting effective viscosity times thickness is regularized by ensuring that
its minimum is at least \f$\epsilon\f$.  This regularization constant is an argument.

In this implementation we set \f$\nu H\f$ to a constant anywhere the ice is
thinner than a certain minimum. See SSAStrengthExtension and compare how this
issue is handled when -cfbc is set.
*/
PetscErrorCode SSAFD::compute_nuH_staggered(IceModelVec2Stag &result, PetscReal epsilon) {
  PetscErrorCode ierr;
  PISMVector2 **uv;

  ierr = result.begin_access(); CHKERRQ(ierr);
  ierr = velocity.get_array(uv); CHKERRQ(ierr);
  ierr = hardness.begin_access(); CHKERRQ(ierr);
  ierr = thickness->begin_access(); CHKERRQ(ierr);

  PetscScalar ssa_enhancement_factor = flow_law->enhancement_factor(),
    n_glen = flow_law->exponent(),
    nu_enhancement_scaling = 1.0 / pow(ssa_enhancement_factor, 1.0/n_glen);

  const PetscScalar dx = grid.dx, dy = grid.dy;

  for (PetscInt o=0; o<2; ++o) {
    const PetscInt oi = 1 - o, oj=o;
    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {

        const PetscScalar H = 0.5 * ((*thickness)(i,j) + (*thickness)(i+oi,j+oj));

        if (H < strength_extension->get_min_thickness()) {
          // Extends strength of SSA (i.e. nuH coeff) into the ice free region.
          //  Does not add or subtract ice mass.
          result(i,j,o) = strength_extension->get_notional_strength();
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

        result(i,j,o) = H * flow_law->effective_viscosity(hardness(i,j,o), u_x, u_y, v_x, v_y);

        if (! finite(result(i,j,o)) || false) {
          ierr = PetscPrintf(grid.com, "nuH[%d][%d][%d] = %e\n", o, i, j, result(i,j,o));
          CHKERRQ(ierr);
          ierr = PetscPrintf(grid.com, "  u_x, u_y, v_x, v_y = %e, %e, %e, %e\n",
                             u_x, u_y, v_x, v_y);
          CHKERRQ(ierr);
        }

        // include the SSA enhancement factor; in most cases ssa_enhancement_factor is 1
        result(i,j,o) *= nu_enhancement_scaling;

        // We ensure that nuH is bounded below by a positive constant.
//OLD WAY before 4/5/12:       result(i,j,o) = PetscMax(epsilon,result(i,j,o));
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


//! Update the nuH viewer, which shows log10(nu H).
PetscErrorCode SSAFD::update_nuH_viewers() {
  PetscErrorCode ierr;
  IceModelVec2S tmp;

  if (!view_nuh) return 0;

  ierr = tmp.create(grid, "nuH", false); CHKERRQ(ierr);
  ierr = tmp.set_attrs("temporary",
                       "log10 of (viscosity * thickness)",
                       "Pa s m", ""); CHKERRQ(ierr);

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

  if (nuh_viewer == PETSC_NULL) {
    ierr = grid.create_viewer(nuh_viewer_size, "nuH", nuh_viewer); CHKERRQ(ierr);
  }

  ierr = tmp.view(nuh_viewer, PETSC_NULL); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode SSAFD::set_diagonal_matrix_entry(Mat A, int i, int j,
                                                PetscScalar value) {
  PetscErrorCode ierr;
  MatStencil row, col;
  row.j = i; row.i = j;
  col.j = i; col.i = j;

  row.c = 0; col.c = 0;
  ierr = MatSetValuesStencil(A, 1, &row, 1, &col, &value, INSERT_VALUES); CHKERRQ(ierr);

  row.c = 1; col.c = 1;
  ierr = MatSetValuesStencil(A, 1, &row, 1, &col, &value, INSERT_VALUES); CHKERRQ(ierr);

  return 0;
}

//! \brief Checks if a cell is near or at the ice front.
/*!
 * You need to call mask->begin_access() before and mask->end_access() after using this.
 *
 * Note that a cell is a CFBC location of one of four direct neighbors is ice-free.
 *
 * If one of the diagonal neighbors is ice-free we don't use the CFBC, but we
 * do need to compute weights used in the SSA discretization (see
 * assemble_matrix()) to avoid differentiating across interfaces between icy
 * and ice-free cells.
 *
 * This method ensures that checks in assemble_rhs() and assemble_matrix() are
 * consistent.
 */
bool SSAFD::is_marginal(int i, int j, bool ssa_dirichlet_bc) {
	
  const PetscInt M_ij = mask->as_int(i,j),
    // direct neighbors
    M_e = mask->as_int(i + 1,j),
    M_w = mask->as_int(i - 1,j),
    M_n = mask->as_int(i,j + 1),
    M_s = mask->as_int(i,j - 1),
    // "diagonal" neighbors
    M_ne = mask->as_int(i + 1,j + 1),
    M_se = mask->as_int(i + 1,j - 1),
    M_nw = mask->as_int(i - 1,j + 1),
    M_sw = mask->as_int(i - 1,j - 1);

  Mask M;

  if (ssa_dirichlet_bc) {
    return (!M.ice_free(M_ij)) &&
      (M.ice_free(M_e) || M.ice_free(M_w) || M.ice_free(M_n) || M.ice_free(M_s) ||
       M.ice_free(M_ne) || M.ice_free(M_se) || M.ice_free(M_nw) || M.ice_free(M_sw));}
  else {
    return (!M.ice_free(M_ij)) &&
      (M.ice_free_ocean(M_e) || M.ice_free_ocean(M_w) || M.ice_free_ocean(M_n) || M.ice_free_ocean(M_s) ||
       M.ice_free_ocean(M_ne) || M.ice_free_ocean(M_se) || M.ice_free_ocean(M_nw) || M.ice_free_ocean(M_sw));
  }
}

SSAFD_nuH::SSAFD_nuH(SSAFD *m, IceGrid &g, PISMVars &my_vars)
  : PISMDiag<SSAFD>(m, g, my_vars) {

  // set metadata:
  dof = 2;
  vars.resize(2);
  vars[0].init_2d("nuH[0]", grid);
  vars[1].init_2d("nuH[1]", grid);

  set_attrs("ice thickness times effective viscosity, i-offset", "",
            "Pa s m", "kPa s m", 0);
  set_attrs("ice thickness times effective viscosity, j-offset", "",
            "Pa s m", "kPa s m", 1);
}

PetscErrorCode SSAFD_nuH::compute(IceModelVec* &output) {
  PetscErrorCode ierr;

  IceModelVec2Stag *result = new IceModelVec2Stag;
  ierr = result->create(grid, "nuH", true); CHKERRQ(ierr);
  ierr = result->set_metadata(vars[0], 0); CHKERRQ(ierr);
  ierr = result->set_metadata(vars[1], 1); CHKERRQ(ierr);
  result->write_in_glaciological_units = true;

  ierr = model->nuH.copy_to(*result); CHKERRQ(ierr);

  output = result;
  return 0;
}

void SSAFD::get_diagnostics(map<string, PISMDiagnostic*> &dict) {
  SSA::get_diagnostics(dict);

  dict["nuH"] = new SSAFD_nuH(this, grid, *variables);
}

