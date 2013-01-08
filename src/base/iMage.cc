// Copyright (C) 2004-2011, 2013 Jed Brown, Ed Bueler and Constantine Khroulev
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

#include <petscdmda.h>
#include "iceModelVec.hh"
#include "columnSystem.hh"
#include "iceModel.hh"
#include "PISMStressBalance.hh"
#include "IceGrid.hh"
#include "pism_options.hh"

//! Tridiagonal linear system for vertical column of age (pure advection) problem.
class ageSystemCtx : public columnSystemCtx {

public:
  ageSystemCtx(PetscInt my_Mz, string my_prefix);
  PetscErrorCode initAllColumns();

  PetscErrorCode solveThisColumn(PetscScalar **x, PetscErrorCode &pivoterrorindex);  

public:
  // constants which should be set before calling initForAllColumns()
  PetscScalar  dx,
               dy,
               dtAge,
               dzEQ;
  // pointers which should be set before calling initForAllColumns()
  PetscScalar  *u,
               *v,
               *w;
  IceModelVec3 *tau3;

protected: // used internally
  PetscScalar nuEQ;
  bool        initAllDone;
};


ageSystemCtx::ageSystemCtx(PetscInt my_Mz, string my_prefix)
      : columnSystemCtx(my_Mz, my_prefix) { // size of system is Mz
  initAllDone = false;
  // set values so we can check if init was called on all
  dx = -1.0;
  dy = -1.0;
  dtAge = -1.0;
  dzEQ = -1.0;
  u = NULL;
  v = NULL;
  w = NULL;
  tau3 = NULL;
}


PetscErrorCode ageSystemCtx::initAllColumns() {
  // check whether each parameter & pointer got set
  if (dx <= 0.0) { SETERRQ(PETSC_COMM_SELF, 2,"un-initialized dx in ageSystemCtx"); }
  if (dy <= 0.0) { SETERRQ(PETSC_COMM_SELF, 3,"un-initialized dy in ageSystemCtx"); }
  if (dtAge <= 0.0) { SETERRQ(PETSC_COMM_SELF, 4,"un-initialized dtAge in ageSystemCtx"); }
  if (dzEQ <= 0.0) { SETERRQ(PETSC_COMM_SELF, 5,"un-initialized dzEQ in ageSystemCtx"); }
  if (u == NULL) { SETERRQ(PETSC_COMM_SELF, 6,"un-initialized pointer u in ageSystemCtx"); }
  if (v == NULL) { SETERRQ(PETSC_COMM_SELF, 7,"un-initialized pointer v in ageSystemCtx"); }
  if (w == NULL) { SETERRQ(PETSC_COMM_SELF, 8,"un-initialized pointer w in ageSystemCtx"); }
  if (tau3 == NULL) { SETERRQ(PETSC_COMM_SELF, 9,"un-initialized pointer tau3 in ageSystemCtx"); }
  nuEQ = dtAge / dzEQ; // derived constant
  initAllDone = true;
  return 0;
}

//! Conservative first-order upwind scheme with implicit in the vertical: one column solve.
/*!
The PDE being solved is
    \f[ \frac{\partial \tau}{\partial t} + \frac{\partial}{\partial x}\left(u \tau\right) + \frac{\partial}{\partial y}\left(v \tau\right) + \frac{\partial}{\partial z}\left(w \tau\right) = 1. \f]
This PDE has the conservative form identified in the comments on IceModel::ageStep().

Let
    \f[ \mathcal{U}(x,y_{i+1/2}) = x \, \begin{Bmatrix} y_i, \quad x \ge 0 \\ y_{i+1}, \quad x \le 0 \end{Bmatrix}. \f]
Note that the two cases agree when \f$x=0\f$, so there is no conflict.  This is
part of the upwind rule, and \f$x\f$ will be the cell-boundary (finite volume sense)
value of the velocity.  Our discretization of the PDE uses this upwind notation 
to build an explicit scheme for the horizontal terms and an implicit scheme for
the vertical terms, as follows.

Let
    \f[ A_{i,j,k}^n \approx \tau(x_i,y_j,z_k) \f]
be the numerical approximation of the exact value on the grid.  The scheme is
\f{align*}{
  \frac{A_{ijk}^{n+1} - A_{ijk}^n}{\Delta t} &+ \frac{\mathcal{U}(u_{i+1/2},A_{i+1/2,j,k}^n) - \mathcal{U}(u_{i-1/2},A_{i-1/2,j,k}^n)}{\Delta x} + \frac{\mathcal{U}(v_{j+1/2},A_{i,j+1/2,k}^n) - \mathcal{U}(v_{j-1/2},A_{i,j-1/2,k}^n)}{\Delta y} \\
    &\qquad \qquad + \frac{\mathcal{U}(w_{k+1/2},A_{i,j,k+1/2}^{n+1}) - \mathcal{U}(w_{k-1/2},A_{i,j,k-1/2}^{n+1})}{\Delta z} = 1.
  \f}
Here velocity components \f$u,v,w\f$ are all evaluated at time \f$t_n\f$, so
\f$u_{i+1/2} = u_{i+1/2,j,k}^n\f$ in more detail, and so on for all the other
velocity values.  Note that this discrete form
is manifestly conservative, in that, for example, the same term at \f$u_{i+1/2}\f$
is used both in updating \f$A_{i,j,k}^{n+1}\f$ and \f$A_{i+1,j,k}^{n+1}\f$.

Rewritten as a system of equations in the vertical index, let
   \f[ \Phi_k = \Delta t - \frac{\Delta t}{\Delta x} \left[\mathcal{U}(u_{i+1/2},A_{i+1/2,j,k}^n) - \mathcal{U}(u_{i-1/2},A_{i-1/2,j,k}^n)\right] - \frac{\Delta t}{\Delta y} \left[\mathcal{U}(v_{j+1/2},A_{i,j+1/2,k}^n) - \mathcal{U}(v_{j-1/2},A_{i,j-1/2,k}^n)\right]. \f]
Let \f$\nu = \Delta t / \Delta z\f$ and for slight simplification denote \f$w_\pm = w_{k\pm 1/2}\f$.  The equation determining the unknown
new age values in the column, denoted \f$a_k = A_{i,j,k}^{n+1}\f$, is
   \f[ a_k + \nu \left[\mathcal{U}(w_+,a_{k+1/2}) - \mathcal{U}(w_-,a_{k-1/2})\right] = A_{ijk}^n + \Phi_k. \f]
This is perhaps easiest to understand as four cases:
\f{align*}{
w_+\ge 0, w_-\ge 0:  && (-\nu w_-) a_{k-1} + (1 + \nu w_+) a_k + (0) a_{k+1} &= A_{ijk}^n + \Phi_k, \\
w_+\ge 0, w_- < 0:   && (0) a_{k-1} + (1 + \nu w_+ - \nu w_-) a_k + (0) a_{k+1} &= A_{ijk}^n + \Phi_k, \\
w_+ < 0,  w_-\ge 0:  && (-\nu w_-) a_{k-1} + (1) a_k + (+\nu w_+) a_{k+1} &= A_{ijk}^n + \Phi_k, \\
w_+ < 0,  w_- < 0:   && (0) a_{k-1} + (1 - \nu w_-) a_k + (+\nu w_+) a_{k+1} &= A_{ijk}^n + \Phi_k.
 \f}
These equation form a tridiagonal system, i.e. the bandwidth does not exceed three.
In every case the on-diagonal coefficient is greater than or equal to one,
while the off-diagonal coefficients are always negative.  The coefficients
approximately sum to one in each case, but only up to errors of size \f$O(\Delta z)\f$.
These facts APPARENTLY imply that the method has a maximum principle \ref MortonMayers.


FIXME:  THE COMMENT ABOVE HAS BEEN UPDATED TO THE 'CONSERVATIVE' FORM, BUT THE
CODE STILL REFLECTS THE OLD SCHEME.

FIXME:  CARE MUST BE TAKEN TO MAINTAIN CONSERVATISM AT SURFACE.
 */
PetscErrorCode ageSystemCtx::solveThisColumn(PetscScalar **x, PetscErrorCode &pivoterrorindex) {
  PetscErrorCode ierr;
  if (!initAllDone) {  SETERRQ(PETSC_COMM_SELF, 2,
     "solveThisColumn() should only be called after initAllColumns() in ageSystemCtx"); }

  // set up system: 0 <= k < ks
  for (PetscInt k = 0; k < ks; k++) {
    planeStar<PetscScalar> ss;  // note ss.ij = tau[k]
    ierr = tau3->getPlaneStar_fine(i,j,k,&ss); CHKERRQ(ierr);
    // do lowest-order upwinding, explicitly for horizontal
    rhs[k] =  (u[k] < 0) ? u[k] * (ss.e -  ss.ij) / dx
                         : u[k] * (ss.ij  - ss.w) / dx;
    rhs[k] += (v[k] < 0) ? v[k] * (ss.n -  ss.ij) / dy
                         : v[k] * (ss.ij  - ss.s) / dy;
    // note it is the age eqn: dage/dt = 1.0 and we have moved the hor.
    //   advection terms over to right:
    rhs[k] = ss.ij + dtAge * (1.0 - rhs[k]);

    // do lowest-order upwinding, *implicitly* for vertical
    PetscScalar AA = nuEQ * w[k];
    if (k > 0) {
      if (AA >= 0) { // upward velocity
        L[k] = - AA;
        D[k] = 1.0 + AA;
        U[k] = 0.0;
      } else { // downward velocity; note  -AA >= 0
        L[k] = 0.0;
        D[k] = 1.0 - AA;
        U[k] = + AA;
      }
    } else { // k == 0 case
      // note L[0] not an allocated location
      if (AA > 0) { // if strictly upward velocity apply boundary condition:
                    // age = 0 because ice is being added to base
        D[0] = 1.0;
        U[0] = 0.0;
        rhs[0] = 0.0;
      } else { // downward velocity; note  -AA >= 0
        D[0] = 1.0 - AA;
        U[0] = + AA;
        // keep rhs[0] as is
      }
    }
  }  // done "set up system: 0 <= k < ks"
      
  // surface b.c. at ks
  if (ks>0) {
    L[ks] = 0;
    D[ks] = 1.0;   // ignore U[ks]
    rhs[ks] = 0.0;  // age zero at surface
  }

  // solve it
  pivoterrorindex = solveTridiagonalSystem(ks+1,x);
  return 0;
}


//! Take a semi-implicit time-step for the age equation.
/*!
Let \f$\tau(t,x,y,z)\f$ be the age of the ice.  Denote the three-dimensional
velocity field within the ice fluid as \f$(u,v,w)\f$.  The age equation
is \f$d\tau/dt = 1\f$, that is, ice may move but it gets one year older in one
year.  Thus
    \f[ \frac{\partial \tau}{\partial t} + u \frac{\partial \tau}{\partial x}
        + v \frac{\partial \tau}{\partial y} + w \frac{\partial \tau}{\partial z} = 1 \f]
This equation is purely advective and hyperbolic.  The right-hand side is "1" as
long as age \f$\tau\f$ and time \f$t\f$ are measured in the same units.
Because the velocity field is incompressible, \f$\nabla \cdot (u,v,w) = 0\f$,
we can rewrite the equation as
    \f[ \frac{\partial \tau}{\partial t} + \nabla \left( (u,v,w) \tau \right) = 1 \f]
There is a conservative first-order numerical method; see ageSystemCtx::solveThisColumn().

The boundary condition is that when the ice falls as snow it has age zero.  
That is, \f$\tau(t,x,y,h(t,x,y)) = 0\f$ in accumulation areas.  There is no 
boundary condition elsewhere on the ice upper surface, as the characteristics
go outward in the ablation zone.  If the velocity in the bottom cell of ice
is upward (\f$w>0\f$) then we also apply a zero age boundary condition,
\f$\tau(t,x,y,0) = 0\f$.  This is the case where ice freezes on at the base,
either grounded basal ice freezing on stored water in till, or marine basal ice.
(Note that the water that is frozen-on as ice might be quite "old" in the sense
that its most recent time in the atmosphere was long ago; this comment is
relevant to any analysis which relates isotope ratios to PISM's modeled age.)

The numerical method is a conservative form of first-order upwinding, but the
vertical advection term is computed implicitly.  Thus there is no CFL-type
stability condition from the vertical velocity; CFL is only for the horizontal
velocity.  We use a finely-spaced, equally-spaced vertical grid in the
calculation.  Note that the IceModelVec3 methods getValColumn...() and
setValColumn..() interpolate back and forth between this fine grid and
the storage grid.  The storage grid may or may not be equally-spaced.  See
ageSystemCtx::solveThisColumn() for the actual method.
 */
PetscErrorCode IceModel::ageStep() {
  PetscErrorCode  ierr;

  // set up fine grid in ice
  PetscInt    fMz = grid.Mz_fine;
  PetscScalar fdz = grid.dz_fine;

  PetscScalar *x;  
  x = new PetscScalar[fMz]; // space for solution

  bool viewOneColumn;
  ierr = PISMOptionsIsSet("-view_sys", viewOneColumn); CHKERRQ(ierr);

  ageSystemCtx system(fMz, "age"); // linear system to solve in each column
  system.dx    = grid.dx;
  system.dy    = grid.dy;
  system.dtAge = dt_TempAge;
  system.dzEQ  = fdz;
  // pointers to values in current column
  system.u     = new PetscScalar[fMz];
  system.v     = new PetscScalar[fMz];
  system.w     = new PetscScalar[fMz];
  // system needs access to tau3 for planeStar()
  system.tau3  = &tau3;
  // this checks that all needed constants and pointers got set
  ierr = system.initAllColumns(); CHKERRQ(ierr);

  IceModelVec3 *u3, *v3, *w3;
  ierr = stress_balance->get_3d_velocity(u3, v3, w3); CHKERRQ(ierr); 

  ierr = vH.begin_access(); CHKERRQ(ierr);
  ierr = tau3.begin_access(); CHKERRQ(ierr);
  ierr = u3->begin_access(); CHKERRQ(ierr);
  ierr = v3->begin_access(); CHKERRQ(ierr);
  ierr = w3->begin_access(); CHKERRQ(ierr);
  ierr = vWork3d.begin_access(); CHKERRQ(ierr);

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      // this should *not* be replaced by a call to grid.kBelowHeight()
      const PetscInt  fks = static_cast<PetscInt>(floor(vH(i,j)/fdz));

      if (fks == 0) { // if no ice, set the entire column to zero age
        ierr = vWork3d.setColumn(i,j,0.0); CHKERRQ(ierr);
      } else { // general case: solve advection PDE; start by getting 3D velocity ...

	ierr = u3->getValColumn(i,j,fks,system.u); CHKERRQ(ierr);
	ierr = v3->getValColumn(i,j,fks,system.v); CHKERRQ(ierr);
	ierr = w3->getValColumn(i,j,fks,system.w); CHKERRQ(ierr);

        ierr = system.setIndicesAndClearThisColumn(i,j,fks); CHKERRQ(ierr);

        // solve the system for this column; call checks that params set
        PetscErrorCode pivoterr;
        ierr = system.solveThisColumn(&x,pivoterr); CHKERRQ(ierr);

        if (pivoterr != 0) {
          ierr = PetscPrintf(PETSC_COMM_SELF,
            "\n\ntridiagonal solve of ageSystemCtx in ageStep() FAILED at (%d,%d)\n"
                " with zero pivot position %d; viewing system to m-file ... \n",
            i, j, pivoterr); CHKERRQ(ierr);
          ierr = system.reportColumnZeroPivotErrorMFile(pivoterr); CHKERRQ(ierr);
          SETERRQ(grid.com, 1,"PISM ERROR in ageStep()\n");
        }
        if (viewOneColumn && issounding(i,j)) {
          ierr = PetscPrintf(PETSC_COMM_SELF,
            "\n\nin ageStep(): viewing ageSystemCtx at (i,j)=(%d,%d) to m-file ... \n\n",
            i, j); CHKERRQ(ierr);
          ierr = system.viewColumnInfoMFile(x, fMz); CHKERRQ(ierr);
        }

        // x[k] contains age for k=0,...,ks, but set age of ice above (and at) surface to zero years
        for (PetscInt k=fks+1; k<fMz; k++) {
          x[k] = 0.0;
        }

        // put solution in IceModelVec3
        ierr = vWork3d.setValColumnPL(i,j,x); CHKERRQ(ierr);
      }
    }
  }

  ierr = vH.end_access(); CHKERRQ(ierr);
  ierr = tau3.end_access();  CHKERRQ(ierr);
  ierr = u3->end_access();  CHKERRQ(ierr);
  ierr = v3->end_access();  CHKERRQ(ierr);
  ierr = w3->end_access();  CHKERRQ(ierr);
  ierr = vWork3d.end_access();  CHKERRQ(ierr);

  delete [] x;
  delete [] system.u;  delete [] system.v;  delete [] system.w;

  ierr = vWork3d.update_ghosts(tau3); CHKERRQ(ierr);

  return 0;
}

