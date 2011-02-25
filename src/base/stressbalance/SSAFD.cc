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

#include "SSAFD.hh"


//! \brief Allocate objects specific to the SSAFD object.
/*!
Because the FD implementation of the FD uses Picard iteration, a PETSc KSP
and Mat are used directly.  In particular we set up a \f$A =\f$
(Mat SSAStiffnessMatrix) and a \f$b = \f$ (Vec SSARHS) and iteratively solve
linear systems
  \f[ A x = b \f]
where \f$x = \f$ (Vec SSAX).  A PETSc SNES object is never created.
 */
PetscErrorCode SSAFD::allocate_fd() {
  PetscErrorCode ierr;

  // note SSADA and SSAX are allocated in SSA::allocate()
  ierr = VecDuplicate(SSAX, &SSARHS); CHKERRQ(ierr);

  ierr = DAGetMatrix(SSADA, MATMPIAIJ, &SSAStiffnessMatrix); CHKERRQ(ierr);

  ierr = KSPCreate(grid.com, &SSAKSP); CHKERRQ(ierr);
  // the default PC type somehow is ILU, which now fails (?) while block jacobi
  //   seems to work; runtime options can override (see test J in vfnow.py)
  PC pc;
  ierr = KSPGetPC(SSAKSP,&pc); CHKERRQ(ierr);
  ierr = PCSetType(pc,PCBJACOBI); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(SSAKSP); CHKERRQ(ierr);

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

  return 0;
}

PetscErrorCode SSAFD::deallocate_fd() {
  PetscErrorCode ierr;

  if (SSAKSP != PETSC_NULL) {
    ierr = KSPDestroy(SSAKSP); CHKERRQ(ierr);
  }

  if (SSAStiffnessMatrix != PETSC_NULL) {
    ierr = MatDestroy(SSAStiffnessMatrix); CHKERRQ(ierr);
  }

  if (SSARHS != PETSC_NULL) {
    ierr = VecDestroy(SSARHS); CHKERRQ(ierr);
  }
  
  return 0;
}


//! \brief Computes the right-hand side ("rhs") of the linear problem for the
//! finite-difference implementation of the SSA equations.
/*! 
The right side of the SSA equations is just the driving stress term
   \f[ - \rho g H \nabla h. \f]
The basal stress is put on the left side of the system.  This method builds the
discrete approximation of the right side.  For more about the discretization
of the SSA equations, see comments for assembleSSAMatrix().

The values of the driving stress on the i,j grid come from a call to
computeDrivingStress().

Grid points with mask value MASK_SHEET correspond to the trivial equations
   \f[ \bar u_{ij} = \frac{uvbar(i-1,j,0) + uvbar(i,j,0)}{2}, \f]
and similarly for \f$\bar v_{ij}\f$.  That is, the vertically-averaged
horizontal velocity is already known for these points because it was either 
computed (on the staggered grid) using the SIA or was set by the -ssaBC
mechanism.
 */
PetscErrorCode SSAFD::assemble_rhs(Vec rhs) {
  PetscErrorCode ierr;
  PISMVector2     **rhs_uv;

  // next constant not too sensitive, but must match value in assembleSSAMatrix():
  const PetscScalar   scaling = 1.0e9;  // comparable to typical beta for an ice stream;

  ierr = VecSet(rhs, 0.0); CHKERRQ(ierr);

  // get driving stress components
  ierr = compute_driving_stress(taud); CHKERRQ(ierr);

  ierr = taud.begin_access(); CHKERRQ(ierr);
  ierr = DAVecGetArray(SSADA,rhs,&rhs_uv); CHKERRQ(ierr);

  if (vel_bc && bc_locations) {
    ierr = vel_bc->begin_access(); CHKERRQ(ierr);
    ierr = bc_locations->begin_access(); CHKERRQ(ierr);
  }

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (vel_bc && (bc_locations->value(i,j) == MASK_BC)) {
        rhs_uv[i][j].u = scaling * (*vel_bc)(i,j).u;
        rhs_uv[i][j].v = scaling * (*vel_bc)(i,j).v;
      } else {
	// usual case: use already computed driving stress
        rhs_uv[i][j].u = taud(i,j).u;
        rhs_uv[i][j].v = taud(i,j).v;
      }
    }
  }

  if (vel_bc) {
    ierr = bc_locations->end_access(); CHKERRQ(ierr);
    ierr = vel_bc->end_access(); CHKERRQ(ierr);
  }
  
  ierr = taud.end_access(); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(SSADA,rhs,&rhs_uv); CHKERRQ(ierr);

  ierr = VecAssemblyBegin(rhs); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(rhs); CHKERRQ(ierr);

  return 0;
}

//! \brief Assemble the left-hand side matrix for the finite difference implementation
//! of the SSA equations.
/*! 
The SSA equations are in their clearest divergence form
    \f[ - \frac{\partial T_{ij}}{\partial x_j} - \tau_{(b)i} = f_i \f]
where \f$i,j\f$ range over \f$x,y\f$, \f$T_{ij}\f$ is a depth-integrated viscous
stress tensor (%i.e. equation (2.6) in [\ref SchoofStream]; also [\ref Morland]).
These equations determine velocity in a more-or-less elliptic manner.
Here \f$\tau_{(b)i}\f$ are the components of the basal shear stress applied to
the base of the ice and \f$f_i\f$ is the driving shear stress,
    \f[ f_i = - \rho g H \frac{\partial h}{\partial x_i}. \f]
Here \f$H\f$ is the ice thickness and \f$h\f$ is the elevation of the surface of
the ice.

More concretely, the SSA equations are
\f{align*}
 - 2 \left[\bar\nu H \left(2 u_x + v_y\right)\right]_x
        - \left[\bar\nu H \left(u_y + v_x\right)\right]_y
        - \tau_{(b)1}  &= - \rho g H h_x, \\
   - \left[\bar\nu H \left(u_y + v_x\right)\right]_x
        - 2 \left[\bar\nu H \left(u_x + 2 v_y\right)\right]_y
        - \tau_{(b)2}  &= - \rho g H h_y, 
\f}
where \f$u\f$ is the \f$x\f$-component of the velocity and \f$v\f$ is the
\f$y\f$-component of the velocity.  Note \f$\bar\nu\f$ is the vertically-averaged
effective viscosity of the ice.

For ice shelves \f$\tau_{(b)i} = 0\f$ [\ref MacAyealetal].
For ice streams with a basal till modelled as a plastic material,
\f$\tau_{(b)i} = - \tau_c u_i/|\mathbf{u}|\f$ where 
\f$\mathbf{u} = (u,v)\f$, \f$|\mathbf{u}| = \left(u^2 + v^2\right)^{1/2}\f$.
Here \f$\tau_c(t,x,y)\f$ is the yield stress of the till [\ref SchoofStream].
More generally, ice streams can be modeled with a pseudo-plastic basal till;
see initBasalTillModel() and updateYieldStressUsingBasalWater() and
[\ref BKAJS].

The pseudo-plastic till model includes all power law sliding relations 
[\ref BKAJS], and in particular it includes modeling the basal till as a linearly-viscous 
material, \f$\tau_{(b)i} = - \beta u_i\f$ where \f$\beta\f$ is the basal drag
(friction) parameter [\ref MacAyeal].  PISM assumes that the basal shear
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
\f$y\f$ differences, and that we write \f$N = \bar\nu\f$ then the first of the 
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
  // next constant not too sensitive, but must match value in assembleSSARhs():
  const PetscScalar   scaling = 1.0e9;  // comparable to typical beta for an ice stream
  IceModelVec2V vel = velocity;         // a shortcut

  ierr = MatZeroEntries(A); CHKERRQ(ierr);

  PetscReal beta_shelves_drag_too = config.get("beta_shelves_drag_too");
  bool shelvesDragToo = config.get_flag("shelves_drag_too");

  /* matrix assembly loop */

  ierr = nuH.begin_access(); CHKERRQ(ierr);
  ierr = tauc->begin_access(); CHKERRQ(ierr);
  ierr = vel.begin_access(); CHKERRQ(ierr);

  if (bc_locations) {
    ierr = bc_locations->begin_access(); CHKERRQ(ierr);
  }

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      PismMask mask_value = MASK_UNKNOWN;

      if (bc_locations)
        mask_value = bc_locations->value(i,j);

      if (mask_value == MASK_BC) {
        // set diagonal entry to one; RHS entry will be known (e.g. SIA) velocity;
        //   this is where boundary value to SSA is set
        MatStencil  row, col;
        row.j = i; row.i = j; row.c = 0;
        col.j = i; col.i = j; col.c = 0;
        ierr = MatSetValuesStencil(A,1,&row,1,&col,&scaling,INSERT_VALUES); CHKERRQ(ierr);
        row.c = 1;
        col.c = 1;
        ierr = MatSetValuesStencil(A,1,&row,1,&col,&scaling,INSERT_VALUES); CHKERRQ(ierr);
      } else {
        const PetscScalar dx2 = dx*dx, d4 = dx*dy*4, dy2 = dy*dy;
        /* Provide shorthand for the following staggered coefficients  nu H:
        *      c11
        *  c00     c01
        *      c10
        * Note that the positive i (x) direction is right and the positive j (y)
        * direction is up. */
        const PetscScalar c00 = nuH(i-1,j,0);
        const PetscScalar c01 = nuH(i,j,0);
        const PetscScalar c10 = nuH(i,j-1,1);
        const PetscScalar c11 = nuH(i,j,1);

        const PetscInt sten = 13;
        MatStencil  row, col[sten];

        /* start with the values at the points */
        PetscScalar valU[] = {
          /*               */ -c11/dy2,
          (2*c00+c11)/d4,     2*(c00-c01)/d4,                 -(2*c01+c11)/d4,
          -4*c00/dx2,         4*(c01+c00)/dx2+(c11+c10)/dy2,  -4*c01/dx2,
          (c11-c10)/d4,                                       (c10-c11)/d4,
          /*               */ -c10/dy2,
          -(2*c00+c10)/d4,    2*(c01-c00)/d4,                 (2*c01+c10)/d4 };
        PetscScalar valV[] = {
          (2*c11+c00)/d4,     (c00-c01)/d4,                   -(2*c11+c01)/d4,
          /*               */ -4*c11/dy2,
          2*(c11-c10)/d4,                                     2*(c10-c11)/d4,
          -c00/dx2,           4*(c11+c10)/dy2+(c01+c00)/dx2,  -c01/dx2,
          -(2*c10+c00)/d4,    (c01-c00)/d4,                   (2*c10+c01)/d4,
          /*               */ -4*c10/dy2 };

        /* Dragging ice experiences friction at the bed determined by the
         *    basalDrag[x|y]() methods.  These may be a plastic, pseudo-plastic,
         *    or linear friction law according to basal->drag(), which gets called
         *    by basalDragx(),basalDragy().  */
        if (include_basal_shear && (mask_value == MASK_DRAGGING_SHEET)) {
          // Dragging is done implicitly (i.e. on left side of SSA eqns for u,v).
          valU[5] += basal.drag((*tauc)(i,j), vel(i,j).u, vel(i,j).v);
          valV[7] += basal.drag((*tauc)(i,j), vel(i,j).u, vel(i,j).v);
        }

        // make shelf drag a little bit if desired
        if (shelvesDragToo && (mask_value == MASK_FLOATING)) {
          //ierr = verbPrintf(1,grid.com,"... SHELF IS DRAGGING ..."); CHKERRQ(ierr);
          valU[5] += beta_shelves_drag_too;
          valV[7] += beta_shelves_drag_too;
        }

        // build "u" equation: NOTE TRANSPOSE
        row.j = i; row.i = j; row.c = 0;
        const PetscInt UI[] = {
          /*       */ i,
          i-1,        i,          i+1,
          i-1,        i,          i+1,
          i-1,                    i+1,
          /*       */ i,
          i-1,        i,          i+1};
        const PetscInt UJ[] = {
          /*       */ j+1,
          j+1,        j+1,        j+1,
          j,          j,          j,
          j,                      j,
          /*       */ j-1,
          j-1,        j-1,        j-1};
        const PetscInt UC[] = {
          /*       */ 0,
          1,          1,          1,
          0,          0,          0,
          1,                      1,
          /*       */ 0,
          1,          1,          1};
        for (PetscInt m=0; m<sten; m++) {
          col[m].j = UI[m]; col[m].i = UJ[m], col[m].c = UC[m];
        }
        ierr = MatSetValuesStencil(A,1,&row,sten,col,valU,INSERT_VALUES); CHKERRQ(ierr);

        // build "v" equation: NOTE TRANSPOSE
        row.j = i; row.i = j; row.c = 1;
        const PetscInt VI[] = {
          i-1,        i,          i+1,
          /*       */ i,
          i-1,                    i+1,
          i-1,        i,          i+1,
          i-1,        i,          i+1,
          /*       */ i};
        const PetscInt VJ[] = {
          j+1,        j+1,        j+1,
          /*       */ j+1,
          j,                      j,
          j,          j,          j,
          j-1,        j-1,        j-1,
          /*       */ j-1};
        const PetscInt VC[] = {
          0,          0,          0,
          /*       */ 1,
          0,                      0,
          1,          1,          1,
          0,          0,          0,
          /*       */ 1};
        for (PetscInt m=0; m<sten; m++) {
          col[m].j = VI[m]; col[m].i = VJ[m], col[m].c = VC[m];
        }
        ierr = MatSetValuesStencil(A,1,&row,sten,col,valV,INSERT_VALUES); CHKERRQ(ierr);

      }
    }
  }

  if (bc_locations) {
    ierr = bc_locations->end_access(); CHKERRQ(ierr);
  }

  ierr = vel.end_access(); CHKERRQ(ierr);  
  ierr = tauc->end_access(); CHKERRQ(ierr);  
  ierr = nuH.end_access(); CHKERRQ(ierr);

  ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

  return 0; 
}


//! \brief Compute the vertically-averaged horizontal velocity from the shallow
//! shelf approximation (SSA).
/*!
This is the main procedure implementing the SSA.  

The outer loop (over index \c k) is the nonlinear iteration.  In this loop the effective 
viscosity is computed by computeEffectiveViscosity() and then the linear system is 
set up and solved.

This procedure creates a PETSC \c KSP, it calls assembleSSAMatrix() and assembleSSARhs() 
to store the linear system in the \c KSP, and then calls the PETSc procedure KSPSolve()
to solve the linear system.  

Solving the linear system is also a loop, an iteration, but it occurs
inside KSPSolve().  This inner loop is controlled by PETSc but the user can set the option 
<tt>-ksp_rtol</tt>, in particular, and that tolerance controls when the inner loop terminates.

Note that <tt>-ksp_type</tt> can be used to choose the \c KSP.  This will set 
which type of linear iterative method is used.  The KSP choice is important
but poorly understood because the eigenvalues of the linearized SSA are not 
well understood.  Nonetheless these eigenvalues determine the convergence of 
this (inner) linear iteration.  The default KSP is GMRES(30).

Note that <tt>-pc_type</tt> will set which preconditioner to use; the default 
is ILU.  A well-chosen preconditioner can put the eigenvalues in the right
place so that the KSP can converge quickly.  The preconditioner is also
important because it will behave differently on different numbers of
processors.  If the user wants the results of SSA calculations to be 
independent of the number of processors, then <tt>-pc_type none</tt> should
be used, but performance will be poor.

If you want to test different KSP methods, it may be helpful to see how many 
iterations were necessary.  Use <tt>-ksp_monitor</tt>.
Initial testing implies that CGS takes roughly half the iterations of 
GMRES(30), but is not significantly faster because the iterations are each 
roughly twice as slow.  Furthermore, ILU and BJACOBI seem roughly equivalent
as preconditioners.

The outer loop terminates when the effective viscosity is no longer changing 
much, according to the tolerance set by the option <tt>-ssa_rtol</tt>.  (The 
outer loop also terminates when a maximum number of iterations is exceeded.)
We save the velocity from the last time step in order to have a better estimate
of the effective viscosity than the u=v=0 result.

In truth there is an "outer outer" loop (over index \c l).  This one manages an
attempt to over-regularize the effective viscosity so that the nonlinear
iteration (the "outer" loop over \c k) has a chance to converge if it doesn't 
converge with the default regularization.
 */
PetscErrorCode SSAFD::solve() {
  PetscErrorCode ierr;
  Mat A = SSAStiffnessMatrix; // solve  A SSAX = SSARHS
  PetscReal   norm, normChange;
  PetscInt    ksp_iterations, outer_iterations;
  KSPConvergedReason  reason;

  stdout_ssa = "";

  bool write_ssa_system_to_matlab = config.get_flag("write_ssa_system_to_matlab");
  PetscReal ssaRelativeTolerance = config.get("ssa_relative_convergence"),
            epsilon              = config.get("epsilon_ssa");

  PetscInt ssaMaxIterations = static_cast<PetscInt>(config.get("max_iterations_ssa"));
  
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

      // assemble (or re-assemble) matrix, which depends on updated viscosity
      ierr = assemble_matrix(true, A); CHKERRQ(ierr);
      if (getVerbosityLevel() > 2)
        stdout_ssa += "A:";

      // call PETSc to solve linear system by iterative method; "inner linear iteration"
      ierr = KSPSetOperators(SSAKSP, A, A, DIFFERENT_NONZERO_PATTERN); CHKERRQ(ierr);
      ierr = KSPSolve(SSAKSP, SSARHS, SSAX); CHKERRQ(ierr); // SOLVE

      // report to standard out about iteration
      ierr = KSPGetConvergedReason(SSAKSP, &reason); CHKERRQ(ierr);
      if (reason < 0) {
        ierr = verbPrintf(1,grid.com, 
            "\n\n\nPISM ERROR:  KSPSolve() reports 'diverged'; reason = %d = '%s';\n"
                  "  see PETSc man page for KSPGetConvergedReason();   ENDING ...\n\n",
            reason,KSPConvergedReasons[reason]); CHKERRQ(ierr);
        PISMEnd();
      }
      ierr = KSPGetIterationNumber(SSAKSP, &ksp_iterations); CHKERRQ(ierr);
      if (getVerbosityLevel() > 2) {
        char tempstr[50] = "";  snprintf(tempstr,50, "S:%d,%d: ", ksp_iterations, reason);
        stdout_ssa += tempstr;
      }

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

      outer_iterations = k + 1;
      if (norm == 0 || normChange / norm < ssaRelativeTolerance) goto done;

    } // end of the "outer loop" (index: k)

    if (epsilon > 0.0) {
       // this has no units; epsilon goes up by this ratio when previous value failed
       const PetscScalar DEFAULT_EPSILON_MULTIPLIER_SSA = 4.0;
       ierr = verbPrintf(1,grid.com,
			 "WARNING: Effective viscosity not converged after %d iterations\n"
			 "\twith epsilon=%8.2e. Retrying with epsilon * %8.2e.\n",
			 ssaMaxIterations, epsilon, DEFAULT_EPSILON_MULTIPLIER_SSA);
       CHKERRQ(ierr);

       ierr = velocity.copy_from(velocity_old); CHKERRQ(ierr);
       epsilon *= DEFAULT_EPSILON_MULTIPLIER_SSA;
    } else {
       SETERRQ1(1, 
         "Effective viscosity not converged after %d iterations; epsilon=0.0.\n"
         "  Stopping.\n", 
         ssaMaxIterations);
    }

  } // end of the "outer outer loop" (index: l)

  done:

  if (getVerbosityLevel() > 2) {
    char tempstr[50] = "";
    snprintf(tempstr,50, "... =%5d outer iterations\n", outer_iterations);
    stdout_ssa += tempstr;
  } else if (getVerbosityLevel() == 2) {
    // at default verbosity, just record last normchange and iterations
    char tempstr[50] = "";
    snprintf(tempstr,50, "%5d outer iterations\n", outer_iterations);
    stdout_ssa += tempstr;
  }
  if (getVerbosityLevel() >= 2)
    stdout_ssa = "  SSA: " + stdout_ssa;
  if (write_ssa_system_to_matlab) {
    ierr = writeSSAsystemMatlab(); CHKERRQ(ierr);
  }

  return 0;
}


//! \brief Write the SSA system to an .m (MATLAB) file (for debugging).
PetscErrorCode SSAFD::writeSSAsystemMatlab() {
  PetscErrorCode ierr;
  PetscViewer    viewer;
  char           file_name[PETSC_MAX_PATH_LEN], yearappend[PETSC_MAX_PATH_LEN];

  IceModelVec2S component;
  ierr = component.create(grid, "temp_storage", false); CHKERRQ(ierr);

  // FIXME: the file name prefix should be an option
  strcpy(file_name,"pism_SSAFD");
  snprintf(yearappend, PETSC_MAX_PATH_LEN, "_y%.0f.m", grid.year);
  strcat(file_name,yearappend);
  ierr = verbPrintf(2, grid.com, 
             "writing Matlab-readable file for SSAFD system A xsoln = rhs to file `%s' ...\n",
             file_name); CHKERRQ(ierr);
  ierr = PetscViewerCreate(grid.com, &viewer);CHKERRQ(ierr);
  ierr = PetscViewerSetType(viewer, PETSC_VIEWER_ASCII);CHKERRQ(ierr);
  ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);CHKERRQ(ierr);
  ierr = PetscViewerFileSetName(viewer, file_name);CHKERRQ(ierr);

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
  ierr = PetscViewerASCIIPrintf(viewer,"\nyear=%10.6f;\n",grid.year);  CHKERRQ(ierr);
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
  ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);

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
PetscErrorCode SSAFD::compute_nuH_staggered(IceModelVec2Stag &result, PetscReal epsilon) {
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

//! \brief Update the nuH viewer.
/*!
 * FIXME FIXME FIXME 
 */
PetscErrorCode SSAFD::update_nuH_viewers() {
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
