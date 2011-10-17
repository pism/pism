// Copyright (C) 2004-2011 Jed Brown, Ed Bueler and Constantine Khroulev
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

#include <petscda.h>
#include "iceModelVec.hh"
#include "columnSystem.hh"
#include "iceModel.hh"
#include "PISMStressBalance.hh"
#include "IceGrid.hh"

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
  if (dx <= 0.0) { SETERRQ(2,"un-initialized dx in ageSystemCtx"); }
  if (dy <= 0.0) { SETERRQ(3,"un-initialized dy in ageSystemCtx"); }
  if (dtAge <= 0.0) { SETERRQ(4,"un-initialized dtAge in ageSystemCtx"); }
  if (dzEQ <= 0.0) { SETERRQ(5,"un-initialized dzEQ in ageSystemCtx"); }
  if (u == NULL) { SETERRQ(6,"un-initialized pointer u in ageSystemCtx"); }
  if (v == NULL) { SETERRQ(7,"un-initialized pointer v in ageSystemCtx"); }
  if (w == NULL) { SETERRQ(8,"un-initialized pointer w in ageSystemCtx"); }
  if (tau3 == NULL) { SETERRQ(9,"un-initialized pointer tau3 in ageSystemCtx"); }
  nuEQ = dtAge / dzEQ; // derived constant
  initAllDone = true;
  return 0;
}


PetscErrorCode ageSystemCtx::solveThisColumn(PetscScalar **x, PetscErrorCode &pivoterrorindex) {
  PetscErrorCode ierr;
  if (!initAllDone) {  SETERRQ(2,
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
The age equation is\f$d\tau/dt = 1\f$, that is,
    \f[ \frac{\partial \tau}{\partial t} + u \frac{\partial \tau}{\partial x}
        + v \frac{\partial \tau}{\partial y} + w \frac{\partial \tau}{\partial z} = 1\f]
where \f$\tau(t,x,y,z)\f$ is the age of the ice and \f$(u,v,w)\f$  is the three dimensional
velocity field.  This equation is purely advective.  And it is hyperbolic.

The boundary condition is that when the ice falls as snow it has age zero.  
That is, \f$\tau(t,x,y,h(t,x,y)) = 0\f$ in accumulation areas, while there is no 
boundary condition elsewhere, as the characteristics go outward in the ablation zone.
(Some more numerical-analytic attention to this is worthwhile.)

If the velocity in the bottom cell of ice is upward (\code (w[i][j][0] > 0 \endcode)
then we also apply an age = 0 boundary condition.  This is the case where ice freezes
on at the base, either grounded basal ice freezing on stored water in till, or marine basal ice.

The numerical method is first-order upwind but the vertical advection term is computed
implicitly.  (Thus there is no CFL-type stability condition for that part.)

We use a finely-spaced, equally-spaced vertical grid in the calculation.  Note that the IceModelVec3 
methods getValColumn...() and setValColumn..() interpolate back and forth between the grid 
on which calculation is done and the storage grid.  Thus the storage grid can be either 
equally spaced or not.
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
          SETERRQ(1,"PISM ERROR in ageStep()\n");
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

  ierr = tau3.beginGhostCommTransfer(vWork3d); CHKERRQ(ierr);
  ierr = tau3.endGhostCommTransfer(vWork3d); CHKERRQ(ierr);

  return 0;
}

