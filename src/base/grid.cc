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
#include "iceModelpreamble.hh"  // just for verbPrintf()

const PetscScalar IceGrid::DEFAULT_ICEPARAM_Lx   = 1500.0e3,
                  IceGrid::DEFAULT_ICEPARAM_Ly   = 1500.0e3,
                  IceGrid::DEFAULT_ICEPARAM_Lz   = 4000.0, 
                  IceGrid::DEFAULT_ICEPARAM_year = 0.0;
const PetscInt    IceGrid::DEFAULT_ICEPARAM_Mx   = 61,
                  IceGrid::DEFAULT_ICEPARAM_My   = 61,
                  IceGrid::DEFAULT_ICEPARAM_Mz   = 31,
                  IceGrid::DEFAULT_ICEPARAM_Mbz  = 1;


IceGrid::IceGrid(MPI_Comm c,
                 PetscMPIInt r,
                 PetscMPIInt s):
  com(c), rank(r), size(s), createDA_done(PETSC_FALSE) { 

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
  Lbz = (Mbz - 1) * (Lz / (Mz - 1));

  zlevels = new PetscScalar[Mz];
  zblevels = new PetscScalar[Mbz];

  spacing_type = 1;
  setVertLevels();
}


IceGrid::~IceGrid() {
  if (destroyDA() != 0) {
    PetscPrintf(com, "IceGrid destructor: invalid destroyDA() return; ENDING\n");
    PetscEnd();
  }
  delete [] zlevels;
  zlevels = PETSC_NULL;
  delete [] zblevels;
  zblevels = PETSC_NULL;
}


PetscErrorCode IceGrid::chooseEquallySpacedVertical() {
  spacing_type = 1;
  return 0;
}


PetscErrorCode IceGrid::chooseChebyshevSpacedVertical() {
  spacing_type = 2;
  return 0;
}


PetscErrorCode IceGrid::chooseQuadraticSpacedVertical() {
  spacing_type = 3;
  return 0;
}


//! Set the vertical levels according to values in Mz, Lz, Mbz, and the spacing_type flag.
/*!
This procedure is only called when a grid is determined from scratch, e.g. by a derived class or when
bootstrapping from 2D data only, but not when reading a model state input file which will have its own grid.

Uses Mz, Lz, Mbz, and spacing_type.  Note that \c spacing_type cannot be zero at entry to this routine;
it must be 1, 2, or 3.

Sets dzMIN, dzMAX, and Lbz.  Sets and re-allocates zlevels[] and zblevels[].
 */
PetscErrorCode  IceGrid::setVertLevels() {
  
  if (spacing_type == 0) {
    SETERRQ(1,"IceGrid::setVertLevels():  spacing_type cannot be zero when running setVertLevels();\n"
               "  was rescale_...() called twice or was choose...VertSpacing() not called?\n");
  } else if ((spacing_type < 0) || (spacing_type > 3)) {
    SETERRQ1(2,"IceGrid::setVertLevels():  spacing_type = %d is invalid\n",spacing_type);
  }
  
  if (Mz < 2) {
    SETERRQ(3,"IceGrid::setVertLevels():  Mz must be at least 2 for setVertLevels to work\n");
  }

  // note lengths of these arrays can change at a rescale
  delete [] zlevels;
  zlevels = new PetscScalar[Mz];

  if (spacing_type == 1) { 
    dzMIN = Lz / ((PetscScalar) Mz - 1);
    dzMAX = dzMIN;
    // equally-spaced in ice
    for (PetscInt k=0; k < Mz - 1; k++) {
      zlevels[k] = dzMIN * ((PetscScalar) k);
    }
    zlevels[Mz - 1] = Lz;  // make sure it is exactly equal
  } else if (spacing_type == 2) { 
    // Spaced according to the Chebyshev extreme points in the interval [0,1], with 1 being 
    // the base of the ice.  That is, we take original Chebyshev extreme points x_j = cos(pi*j/N), 
    // with N = 2 Mz - 1 but choosing only j=0,1,...,Mz-1 (and not j=0,1,...,N).  These satisfy
    // 0 <= x_j <= 1 and they are clustered near x=1.  Then we flip and scale:
    //   z_j = Lz * (1 - x_z).
    // Note that the smallest spacing is a factor proportional to Mz smaller than dzEQ
    // (i.e. z_1 = C Lz / (Mz^2) while dzEQ = Lz/(Mz-1)).   Near the top the spacing is 
    // to very good approximation equal to pi * dzEQ.
    for (PetscInt k=0; k < Mz - 1; k++) {
      zlevels[k] = Lz * ( 1.0 - cos((pi/2.0) * k / (Mz-1)) );
    }
    zlevels[Mz - 1] = Lz;  // make sure it is exactly equal
    dzMIN = zlevels[1] - zlevels[0];
    dzMAX = zlevels[Mz-1] - zlevels[Mz-2];
  } else if (spacing_type == 3) { 
    // this quadratic scheme is an attempt to be less extreme in the fineness near the base.
    const PetscScalar  ll = 3.0;
    for (PetscInt k=0; k < Mz - 1; k++) {
      const PetscScalar zeta = ((PetscScalar) k) / ((PetscScalar) Mz - 1);
      zlevels[k] = Lz * ( (zeta / ll) * (1.0 + (ll - 1.0) * zeta) );
    }
    zlevels[Mz - 1] = Lz;  // make sure it is exactly equal
    dzMIN = zlevels[1] - zlevels[0];
    dzMAX = zlevels[Mz-1] - zlevels[Mz-2];
  }

  delete [] zblevels;
  zblevels = new PetscScalar[Mbz];

  // equally-spaced in bedrock with spacing determined by Lz and Mz;
  //   also set
  if (Mbz == 1) {
    zblevels[0] = 0.0;
    Lbz = 0.0;
  } else if (Mbz > 1) {
    const PetscScalar dzEQ = Lz / ((PetscScalar) Mz - 1);
    Lbz = dzEQ * ((PetscScalar) Mbz - 1);
    for (PetscInt kb=0; kb < Mbz; kb++) {
      zblevels[kb] = -Lbz + dzEQ * ((PetscScalar) kb);
    }
  } else {
    SETERRQ(4,"Mbz must be at least one");
  }

  return 0;
}


PetscErrorCode IceGrid::printVertLevels(const int verbosity) {
  PetscErrorCode ierr;
  ierr = verbPrintf(verbosity,com,  "    printing vertical levels in ice (Mz=%d,Lz=%5.4f): ",Mz,Lz);
     CHKERRQ(ierr);
  for (PetscInt k=0; k < Mz; k++) {
    ierr = verbPrintf(verbosity,com," %5.4f,",zlevels[k]); CHKERRQ(ierr);
  }
  ierr = verbPrintf(verbosity,com,"\n    printing vertical levels in bedrock (Mbz=%d,Lbz=%5.4f): ",Mbz,Lbz);
     CHKERRQ(ierr);
  for (PetscInt kb=0; kb < Mbz; kb++) {
    ierr = verbPrintf(verbosity,com," %5.4f,",zblevels[kb]); CHKERRQ(ierr);
  }
  ierr = verbPrintf(verbosity,com,"\n"); CHKERRQ(ierr);
  return 0;
}


bool IceGrid::isEqualVertSpacing() {
  return (spacing_type == 1);
}

  
PetscInt IceGrid::kBelowHeight(const PetscScalar height) {
  if (height < 0.0 - 1.0e-6) {
    PetscPrintf(com, 
       "IceGrid kBelowHeight(): height = %5.4f is below base of ice (height must be nonnegative)\n",
       height);
    PetscEnd();
  }
  if (height > Lz + 1.0e-6) {
    PetscPrintf(com, 
       "IceGrid kBelowHeight(): height = %5.4f is above top of computational grid Lz = %5.4f\n",
       height,Lz);
    PetscEnd();
  }
  PetscInt mcurr = 0;
  while (zlevels[mcurr+1] < height) {
    mcurr++;
  }
  return mcurr;
}


//! Rescale IceGrid based on new values for \c Lx, \c Ly, \c Lz (default version).
/*! 
This method computes \c dx, \c dy, \c dz, and \c Lbz based on the current values of \c Mx,
\c My, \c Mz, \c Mbz and the input values of \c lx, \c ly, and \c lz.  It also sets \c Lx = \c lx, etc.

Uses <tt>spacing_type</tt> to set the vertical levels.  In particular, one of the procedures 
chooseEquallySpacedVertical(), chooseQuadraticSpacedVertical(), or chooseChebyshevSpacedVertical() 
should be called before this one.

If vertical levels are already determined, use method rescale_using_zlevels().

See the comment for rescale_and_set_zlevels(<tt>lx,ly,lz,truelyPeriodic</tt>).
 */
PetscErrorCode IceGrid::rescale_and_set_zlevels(const PetscScalar lx, const PetscScalar ly, 
                                                const PetscScalar lz) {
  PetscErrorCode ierr;
  
  ierr = rescale_and_set_zlevels(lx,ly,lz,PETSC_FALSE); CHKERRQ(ierr);
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
PetscErrorCode IceGrid::rescale_and_set_zlevels(const PetscScalar lx, const PetscScalar ly, 
                                                const PetscScalar lz, const PetscTruth truelyPeriodic) {
  PetscErrorCode ierr;

  if (lz<=0) {
    SETERRQ(1, "rescale: lz must be positive\n");
  }
  if (Mz <= 1) {
    SETERRQ(2, "rescale: Mz must be at least 2 (ice part of computational domain must have thickness)\n");
  }

  Lz = lz;  
  ierr = setVertLevels(); CHKERRQ(ierr);  // note Lbz is determined within setVertLevels()

  ierr = rescale_using_zlevels(lx, ly, truelyPeriodic); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceGrid::rescale_using_zlevels(const PetscScalar lx, const PetscScalar ly) {
  PetscErrorCode ierr;
  ierr = rescale_using_zlevels(lx,ly,PETSC_FALSE);
  return 0;
}


//! Rescale IceGrid based on new values for \c Lx, \c Ly but assuming zlevels and zblevels already have correct values.
/*! 
Before calling this, always make sure Mz, Mbz, zlevels[], and zblevels[] are correctly set
and consistent.  In particular, they need to be strictly increasing and zlevels[0] = zblevels[Mbz-1] = 0.0
must be true.
 */
PetscErrorCode IceGrid::rescale_using_zlevels(const PetscScalar lx, const PetscScalar ly, 
                                              const PetscTruth truelyPeriodic) {
  PetscErrorCode ierr;

  if (lx<=0) {
    SETERRQ(1, "rescale: lx must be positive\n");
  }
  if (ly<=0) {
    SETERRQ(2, "rescale: ly must be positive\n");
  }

  Lx = lx; Ly = ly;
  if (truelyPeriodic == PETSC_TRUE) {
    dx = 2.0 * Lx / (Mx);
    dy = 2.0 * Ly / (My);
  } else {
    dx = 2.0 * Lx / (Mx - 1);
    dy = 2.0 * Ly / (My - 1);
  }

  if ( (!isIncreasing(Mz, zlevels)) || (PetscAbs(zlevels[0]) > 1.0e-10) ) {
    SETERRQ(3, "rescale: zlevels invalid; must be strictly increasing and start with z=0\n");
  }
  zlevels[0] = 0.0;
  if ( (!isIncreasing(Mbz, zblevels)) || (PetscAbs(zblevels[Mbz-1]) > 1.0e-10) ) {
    SETERRQ(3, "rescale: zlevels invalid; must be strictly increasing and start with z=0\n");
  }
  zblevels[Mbz-1] = 0.0;
  
  Lz = zlevels[Mz-1];
  Lbz = - zblevels[0];

  ierr = get_dzMIN_dzMAX_spacingtype(); CHKERRQ(ierr);
  
  // it is not known if the following has any effect:
  ierr = DASetUniformCoordinates(da2, -Ly, Ly, -Lx, Lx,
                                 PETSC_NULL, PETSC_NULL); CHKERRQ(ierr);
  return 0;
}


bool IceGrid::isIncreasing(const PetscInt len, PetscScalar *vals) {
  for (PetscInt k = 0; k < len-1; k++) {
    if (vals[k] >= vals[k+1])  return false;
  }
  return true;
}


//! From given vertical grid zlevels[], determine \c dzMIN, \c dzMAX and determine whether grid is equally-spaced.
/*!
The standard for equal spacing is \f$10^{-8}\f$ m max difference between \c dzMIN and \c dzMAX.
 */
PetscErrorCode IceGrid::get_dzMIN_dzMAX_spacingtype() {
  dzMIN = Lz; 
  dzMAX = 0.0;
  for (PetscInt k = 0; k < Mz - 1; k++) {
    const PetscScalar mydz = zlevels[k+1] - zlevels[k];
    dzMIN = PetscMin(mydz,dzMIN);
    dzMAX = PetscMax(mydz,dzMAX);
  }
  if (PetscAbs(dzMAX - dzMIN) <= 1.0e-8) {
    spacing_type = 1;
  } else {
    spacing_type = 0;
  }
  return 0;
}


//! Create the PETSc DA \c da2 for the horizontal grid.  Determine how the horizontal grid is divided among processors.
/*!
This procedure should only be called 
after the parameters describing the horizontal computational box (Lx,Ly) and the parameters for 
the horizontal grid (Mx,My) are already determined.  In particular, the input file (either 
\c -if or \c -bif) and user options (like \c -Mx) must have already been read to determine 
the parameters, and any conflicts must have been resolved.
 */
PetscErrorCode IceGrid::createDA() {
  PetscErrorCode ierr;

  if (createDA_done == PETSC_TRUE) {
    ierr = destroyDA(); CHKERRQ(ierr);
  }

  ierr = DACreate2d(com, DA_XYPERIODIC, DA_STENCIL_BOX,
                    My, Mx, PETSC_DECIDE, PETSC_DECIDE, 1, 1,
                    PETSC_NULL, PETSC_NULL, &da2); CHKERRQ(ierr);

/* see IceModelVec3; we no longer need da3 and da3b
  PetscInt    M, N, m, n;
  ierr = DAGetInfo(da2, PETSC_NULL, &N, &M, PETSC_NULL, &n, &m, PETSC_NULL,
                   PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL); CHKERRQ(ierr);
  ierr = DACreate3d(com, DA_YZPERIODIC, DA_STENCIL_STAR, Mz, N, M, 1, n, m, 1, 1,
                    PETSC_NULL, PETSC_NULL, PETSC_NULL, &da3); CHKERRQ(ierr);
  ierr = DACreate3d(com, DA_YZPERIODIC, DA_STENCIL_STAR, Mbz, N, M, 1, n, m, 1, 1,
                    PETSC_NULL, PETSC_NULL, PETSC_NULL, &da3b); CHKERRQ(ierr);
*/

  DALocalInfo info;
  ierr = DAGetLocalInfo(da2, &info); CHKERRQ(ierr);
  xs = info.ys; xm = info.ym;
  ys = info.xs; ym = info.xm;

  createDA_done = PETSC_TRUE;
  return 0;
}


PetscErrorCode IceGrid::destroyDA() {
  PetscErrorCode ierr;

  ierr = DADestroy(da2); CHKERRQ(ierr);
  return 0;
}

