// Copyright (C) 2004-2008 Jed Brown and Ed Bueler
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

#include <petscfix.h>
#include <petscda.h>
#include "grid.hh"


const PetscScalar IceParam::DEFAULT_ICEPARAM_Lx   = 1500.0e3,
                  IceParam::DEFAULT_ICEPARAM_Ly   = 1500.0e3,
                  IceParam::DEFAULT_ICEPARAM_Lz   = 4000.0, 
                  IceParam::DEFAULT_ICEPARAM_year = 0.0;
const PetscInt    IceParam::DEFAULT_ICEPARAM_Mx   = 61,
                  IceParam::DEFAULT_ICEPARAM_My   = 61,
                  IceParam::DEFAULT_ICEPARAM_Mz   = 31,
                  IceParam::DEFAULT_ICEPARAM_Mbz  = 1;

IceParam::IceParam() {
  Lx = DEFAULT_ICEPARAM_Lx;
  Ly = DEFAULT_ICEPARAM_Ly;
  Lz = DEFAULT_ICEPARAM_Lz;
  year = DEFAULT_ICEPARAM_year;
  Mx = DEFAULT_ICEPARAM_Mx;
  My = DEFAULT_ICEPARAM_My;
  Mz = DEFAULT_ICEPARAM_Mz;
  Mbz = DEFAULT_ICEPARAM_Mbz;
  dx = 2.0 * Lx / (Mx - 1);
  dy = 2.0 * Ly / (My - 1);
//  dz = Lz / (Mz - 1);
  Lbz = (Mbz - 1) * (Lz / (Mz - 1));
}


IceGrid::IceGrid(MPI_Comm c,
                 PetscMPIInt r,
                 PetscMPIInt s):
  com(c), rank(r), size(s), createDA_done(PETSC_FALSE) { 
  p = new IceParam;
  zlevels = new PetscScalar[p->Mz];
  zblevels = new PetscScalar[p->Mbz];
  zlevelsEQ = new PetscScalar[p->Mz];
  zblevelsEQ = new PetscScalar[p->Mbz];

  setLevelsFromLsMs();
  equally_spaced = PETSC_TRUE;
}


IceGrid::IceGrid(MPI_Comm c,
                 PetscMPIInt r,
                 PetscMPIInt s,
                 IceParam *prm):
  com(c), rank(r), size(s), p(prm), createDA_done(PETSC_FALSE) { }


IceGrid::~IceGrid() {
  if (destroyDA() != 0) {
    PetscPrintf(com, "IceGrid destructor: invalid destroyDA() return; ENDING\n");
    PetscEnd();
  }
  delete p;
  p = PETSC_NULL;
  delete [] zlevels;
  zlevels = PETSC_NULL;
  delete [] zblevels;
  zblevels = PETSC_NULL;
  delete [] zlevelsEQ;
  zlevelsEQ = PETSC_NULL;
  delete [] zblevelsEQ;
  zlevelsEQ = PETSC_NULL;
}


PetscErrorCode  IceGrid::setLevelsFromLsMs() {

  // note lengths of these arrays can change at a rescale
  delete [] zlevels;
  delete [] zblevels;
  delete [] zlevelsEQ;
  delete [] zblevelsEQ;
  zlevels = new PetscScalar[p->Mz];
  zblevels = new PetscScalar[p->Mbz];
  zlevelsEQ = new PetscScalar[p->Mz];
  zblevelsEQ = new PetscScalar[p->Mbz];

  dzEQ = p->Lz / (p->Mz - 1);
  dzbEQ = p->Lbz / (p->Mbz - 1);

// EQUAL:
  for (PetscInt k=0; k < p->Mz; k++) {
    zlevelsEQ[k] = dzEQ * ((PetscScalar) k);
    zlevels[k] = zlevelsEQ[k];
  }
  for (PetscInt kb=0; kb < p->Mbz; kb++) {
    zblevelsEQ[kb] = - p->Lbz + dzbEQ * ((PetscScalar) kb);
    zblevels[kb] = zblevelsEQ[kb];
  }
  return 0;
}



//! Create the PETSc DAs for the grid specified in *p (a pointer to IceParam).
PetscErrorCode IceGrid::createDA() {
  PetscErrorCode ierr;

  if (createDA_done == PETSC_TRUE) {
    ierr = destroyDA(); CHKERRQ(ierr);
  }

  ierr = DACreate2d(com, DA_XYPERIODIC, DA_STENCIL_BOX,
                    p->My, p->Mx, PETSC_DECIDE, PETSC_DECIDE, 1, 1,
                    PETSC_NULL, PETSC_NULL, &da2); CHKERRQ(ierr);

/* see IceModelVec3:
  PetscInt    M, N, m, n;
  ierr = DAGetInfo(da2, PETSC_NULL, &N, &M, PETSC_NULL, &n, &m, PETSC_NULL,
                   PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL); CHKERRQ(ierr);
  ierr = DACreate3d(com, DA_YZPERIODIC, DA_STENCIL_STAR, p->Mz, N, M, 1, n, m, 1, 1,
                    PETSC_NULL, PETSC_NULL, PETSC_NULL, &da3); CHKERRQ(ierr);
  ierr = DACreate3d(com, DA_YZPERIODIC, DA_STENCIL_STAR, p->Mbz, N, M, 1, n, m, 1, 1,
                    PETSC_NULL, PETSC_NULL, PETSC_NULL, &da3b); CHKERRQ(ierr);
*/

  DALocalInfo info;
  ierr = DAGetLocalInfo(da2, &info); CHKERRQ(ierr);
  xs = info.ys; xm = info.ym;
  ys = info.xs; ym = info.xm;

  ierr = rescale(p->Lx,p->Ly,p->Lz); CHKERRQ(ierr);
  createDA_done = PETSC_TRUE;
  return 0;
}


PetscErrorCode IceGrid::destroyDA() {
  PetscErrorCode ierr;

  ierr = DADestroy(da2); CHKERRQ(ierr);
  return 0;
}


bool IceGrid::equalVertSpacing() {
  return (equally_spaced == PETSC_TRUE);
}

  
PetscInt IceGrid::kBelowHeight(const PetscScalar height) {
  if (height < 0.0 - 1.0e-6) {
    SETERRQ1(1,"height = %5.4f is below base of ice (height must be nonnegative)\n",
             height);
  }
  if (height > p->Lz + 1.0e-6) {
    SETERRQ2(2,"height = %20.18f is above top of computational grid Lz = %20.18f\n",
             height,p->Lz);
  }
  PetscInt mcurr = 0;
  while (zlevels[mcurr+1] < height) {
    mcurr++;
  }
  return mcurr;
}


PetscInt IceGrid::kBelowHeightEQ(const PetscScalar height) {
  if (height < 0.0 - 1.0e-6) {
    SETERRQ1(1,"height = %5.4f is below base of ice (height must be nonnegative)\n",
             height);
  }
  const PetscInt ks = static_cast<PetscInt>(floor(height/dzEQ));
  if (ks >= p->Mz) {
    SETERRQ2(2,"height = %20.18f is more than dzEQ above top of computational grid Lz = %20.18f\n",
             height,p->Lz);
  }
  return ks;
}


//! Rescale IceGrid based on new values for \c Lx, \c Ly, \c Lz (default version).
/*! 
This method computes \c dx, \c dy, \c dz, and \c Lbz based on the current values of \c Mx,
\c My, \c Mz, \c Mbz and the input values of \c Lx, \c Ly, and \c Lz.  Note that \c dz
is the vertical spacing both in the ice and in bedrock.  Thus <tt>Lbz = Mbz * dz</tt> 
determines \c Lbz.

See the comment for <tt>rescale(lx,ly,lz,truelyPeriodic)</tt>.
 */
PetscErrorCode IceGrid::rescale(const PetscScalar lx, const PetscScalar ly, 
                                const PetscScalar lz) {
  PetscErrorCode ierr;
  
  ierr = rescale(lx,ly,lz,PETSC_FALSE); CHKERRQ(ierr);
  return 0;
}


//! Rescale IceGrid based on new values for \c Lx, \c Ly,\c Lz, but optionally allowing for periodicity.
/*! 
The grid used in PISM, in particular the PETSc DAs used here, are periodic in x and y.
This means that the ghosted values <tt> foo[i+1][j], foo[i-1][j], foo[i][j+1], foo[i][j-1]</tt>
for all 2D Vecs, and similarly in the x and y directions for 3D Vecs, are always available.
That is, they are available even if i,j is a point at the edge of the grid.  On the other
hand, by default, \c dx  is the full width  <tt>2 * Lx</tt>  divided by  <tt>Mx - 1</tt>.
This means that we conceive of the computational domain as starting at the <tt>i = 0</tt>
grid location and ending at the  <tt>i = Mx - 1</tt>  grid location, in particular.  
This idea is not quite compatible with the periodic nature of the grid.

The upshot is that if one computes in a truely periodic way then the gap between the  
<tt>i = 0</tt>  and  <tt>i = Mx - 1</tt>  grid points should \em also have width  \c dx.  
Thus we compute  <tt>dx = 2 * Lx / Mx</tt>.
*/
PetscErrorCode IceGrid::rescale(const PetscScalar lx, const PetscScalar ly, 
                                const PetscScalar lz, const PetscTruth truelyPeriodic) {
  PetscErrorCode ierr;

  if (lz<=0) {
    SETERRQ(1, "rescale: error, lz must be positive\n");
  }

  p->Lz = lz;
  if (p->Mz <= 1) {
    SETERRQ(2, "rescale: error, Mz must be at least 2 (i.e. computational domain must have thickness)\n");
  }
  p->Lbz = ( p->Lz / (p->Mz - 1) ) * (p->Mbz - 1);  // make dz in bed equal to dz in ice; gives Lbz=0 if Mbz=1
  
  ierr = setLevelsFromLsMs(); CHKERRQ(ierr);

  ierr = rescale_using_zlevels(lx, ly, truelyPeriodic); CHKERRQ(ierr);
  return 0;
}


//! Rescale IceGrid based on new values for \c Lx, \c Ly but assuming zlevels and zblevels already have correct values.
PetscErrorCode IceGrid::rescale_using_zlevels(const PetscScalar lx, const PetscScalar ly, 
                                              const PetscTruth truelyPeriodic) {
  PetscErrorCode ierr;

  if (lx<=0) {
    SETERRQ(1, "rescale: error, lx must be positive\n");
  }
  if (ly<=0) {
    SETERRQ(2, "rescale: error, ly must be positive\n");
  }

  p->Lx = lx; p->Ly = ly;
  if (truelyPeriodic == PETSC_TRUE) {
    p->dx = 2.0 * p->Lx / (p->Mx);
    p->dy = 2.0 * p->Ly / (p->My);
  } else {
    p->dx = 2.0 * p->Lx / (p->Mx - 1);
    p->dy = 2.0 * p->Ly / (p->My - 1);
  }

  if ( (!isIncreasing(p->Mz, zlevels)) || (PetscAbs(zlevels[0]) > 1.0e-10) ) {
    SETERRQ(3, "rescale: zlevels invalid; must be strictly increasing and start with z=0\n");
  }
  zlevels[0] = 0.0;
  if ( (!isIncreasing(p->Mbz, zblevels)) || (PetscAbs(zblevels[p->Mbz-1]) > 1.0e-10) ) {
    SETERRQ(3, "rescale: zlevels invalid; must be strictly increasing and start with z=0\n");
  }
  zblevels[p->Mbz-1] = 0.0;
  
  p->Lz = zlevels[p->Mz-1];
  p->Lbz = - zblevels[0];

  dzEQ = p->Lz / (p->Mz - 1);
  if (p->Mbz > 1) {
    dzbEQ = p->Lbz / (p->Mbz - 1);
  } else {
    dzbEQ = dzEQ;  // FIXME: to avoid problem in temperatureStep()
  }

  delete [] zlevelsEQ;
  delete [] zblevelsEQ;
  zlevelsEQ = new PetscScalar[p->Mz];
  zblevelsEQ = new PetscScalar[p->Mbz];
  for (PetscInt k=0; k < p->Mz; k++) {
    zlevelsEQ[k] = dzEQ * ((PetscScalar) k);
  }
  for (PetscInt kb=0; kb < p->Mbz; kb++) {
    zblevelsEQ[kb] = - p->Lbz + dzbEQ * ((PetscScalar) kb);
  }

  // it is not known if the following has any effect:
  ierr = DASetUniformCoordinates(da2, -p->Ly, p->Ly, -p->Lx, p->Lx,
                                 PETSC_NULL, PETSC_NULL); CHKERRQ(ierr);
  return 0;
}


bool IceGrid::isIncreasing(const PetscInt len, PetscScalar *vals) {
  for (PetscInt k = 0; k < len-1; k++) {
    if (vals[k] >= vals[k+1])  return false;
  }
  return true;
}

