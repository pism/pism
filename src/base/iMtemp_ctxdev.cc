// Copyright (C) 2004-2009 Jed Brown, Ed Bueler and Constantine Khroulev
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

#include <petsc.h>
#include "iceModelVec.hh"
#include "iceModel.hh"
#include "columnSystem.hh"


//! Take a semi-implicit time-step for the age equation.  Also check the horizontal CFL for advection.
/*!
The age equation is\f$d\tau/dt = 1\f$, that is,
    \f[ \frac{\partial \tau}{\partial t} + u \frac{\partial \tau}{\partial x}
        + v \frac{\partial \tau}{\partial y} + w \frac{\partial \tau}{\partial z} = 1\f]
where \f$\tau(t,x,y,z)\f$ is the age of the ice and \f$(u,v,w)\f$  is the three dimensional
velocity field.  This equation is hyperbolic (purely advective).  
The boundary condition is that when the ice fell as snow it had age zero.  
That is, \f$\tau(t,x,y,h(t,x,y)) = 0\f$ in accumulation areas, while there is no 
boundary condition elsewhere (as the characteristics go outward elsewhere).

If the velocity in the bottom cell of ice is upward (w[i][j][0]>0) then we apply
an age=0 boundary condition.  This is the case where ice freezes on at the base,
either grounded basal ice or marine basal ice.

A related matter:  By default, when computing the grain size for the 
Goldsby-Kohlstedt flow law, the age \f$\tau\f$ is not used.  Instead a pseudo age 
is computed by updateGrainSizeNow().  If you want the age computed by this routine 
to be used for the grain size estimation, 
from the Vostok core relation as in grainSizeVostok(), add option 
<tt>-real_age_grainsize</tt>.

\latexonly\index{BOMBPROOF!implementation for age equation}\endlatexonly
The numerical method is first-order upwind but the vertical advection term is computed
implicitly.  Thus there is no CFL-type stability condition for that part.

We use equally-spaced vertical grid in the calculation.  Note that the IceModelVec3 
methods getValColumn() and setValColumn() interpolate back and forth between the grid 
on which calculation is done and the storage grid.  Thus the storage grid can be either 
equally spaced or not.
 */
#if 1
PetscErrorCode IceModel::ageStep(PetscScalar* CFLviol) {
  PetscErrorCode  ierr;

  PetscInt    Mz, dummyMbz;
  PetscScalar dzEQ, dummydz, *zlevEQ, *dummylev;

  ierr = getMzMbzForTempAge(Mz, dummyMbz); CHKERRQ(ierr);

  zlevEQ = new PetscScalar[Mz];
  dummylev = new PetscScalar[dummyMbz];

  ierr = getVertLevsForTempAge(Mz, dummyMbz, dzEQ, dummydz, zlevEQ, dummylev);
     CHKERRQ(ierr);

  delete [] dummylev;

  PetscScalar **H, *tau, *u, *v, *w;
  tau = new PetscScalar[Mz];
  u = new PetscScalar[Mz];
  v = new PetscScalar[Mz];
  w = new PetscScalar[Mz];

  PetscScalar *x;  
  x = new PetscScalar[Mz]; // space for solution

  ageSystemCtx system(Mz); // linear system to solve in each column

  const PetscScalar cflx = grid.dx / dtTempAge,
                    cfly = grid.dy / dtTempAge,
                    nuEQ = dtTempAge / dzEQ;

  // set constants within system of equations which are independent of column
  ierr = system.ageSetConstants(grid.dx,grid.dy,dtTempAge,nuEQ,zlevEQ); CHKERRQ(ierr);

  ierr = vH.get_array(H); CHKERRQ(ierr);
  ierr = tau3.begin_access(); CHKERRQ(ierr);
  ierr = u3.begin_access(); CHKERRQ(ierr);
  ierr = v3.begin_access(); CHKERRQ(ierr);
  ierr = w3.begin_access(); CHKERRQ(ierr);
  ierr = taunew3.begin_access(); CHKERRQ(ierr);

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      // re task #4218: keep this check despite creation of the grid extension
      // mechanism; check here is slightly different; later this one can be
      // removed if never active.

      // this should *not* be replaced by a call to grid.kBelowHeightEQ()
      const PetscInt  ks = static_cast<PetscInt>(floor(H[i][j]/dzEQ));
      if (ks > Mz-1) {
        SETERRQ3(1,
           "ageStep() ERROR: ks = %d too high in ice column;\n"
           "  H[i][j] = %5.4f exceeds Lz = %5.4f\n",
           ks, H[i][j], grid.Lz);
      }

      if (ks == 0) { // if no ice, set the entire column to zero age
                     // and ignore the velocities in that column
        ierr = taunew3.setColumn(i,j,0.0); CHKERRQ(ierr);
      } else { // general case
        ierr = tau3.getValColumnQUAD(i,j,Mz,zlevEQ,tau); CHKERRQ(ierr);
        ierr = u3.getValColumnQUAD(i,j,Mz,zlevEQ,u); CHKERRQ(ierr);
        ierr = v3.getValColumnQUAD(i,j,Mz,zlevEQ,v); CHKERRQ(ierr);
        ierr = w3.getValColumnQUAD(i,j,Mz,zlevEQ,w); CHKERRQ(ierr);

        // age evolution is pure advection (so provides check on temp calculation):
        //   check horizontal CFL conditions at each point
        for (PetscInt k=0; k<ks; k++) {
          if (PetscAbs(u[k]) > cflx)  *CFLviol += 1.0;
          if (PetscAbs(v[k]) > cfly)  *CFLviol += 1.0;
        }

        // set up and solve the system for this column
        ierr = system.ageColumnWiseSetUpAndSolve(i,j,ks,u,v,w,tau3,&x);
        if (ierr > 0) {
          SETERRQ3(2,
            "Tridiagonal solve failed at (%d,%d) with zero pivot position %d.\n",
            i, j, ierr);
        } else { CHKERRQ(ierr); }

        // x[k] contains age for k=0,...,ks
        for (PetscInt k=ks+1; k<Mz; k++) {
          x[k] = 0.0;  // age of ice above (and at) surface is zero years
        }
        
        // put solution in IceModelVec3
        ierr = taunew3.setValColumnPL(i,j,Mz,zlevEQ,x); CHKERRQ(ierr);
      }
    }
  }

  ierr = vH.end_access(); CHKERRQ(ierr);
  ierr = tau3.end_access();  CHKERRQ(ierr);
  ierr = u3.end_access();  CHKERRQ(ierr);
  ierr = v3.end_access();  CHKERRQ(ierr);
  ierr = w3.end_access();  CHKERRQ(ierr);
  ierr = taunew3.end_access();  CHKERRQ(ierr);

  delete [] x;

  delete [] tau;  delete [] u;  delete [] v;  delete [] w;

  delete [] zlevEQ;

  return 0;
}
#endif

