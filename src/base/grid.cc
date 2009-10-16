// Copyright (C) 2004-2009 Jed Brown and Ed Bueler
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
#include "pism_const.hh"
#include "nc_util.hh"
#include "../udunits/udunits.h"

//! Use equally-spaced vertical by default.
const SpacingType IceGrid::DEFAULT_SPACING_TYPE = EQUAL;

//! Use non-periodic grid by default.
const Periodicity IceGrid::DEFAULT_PERIODICITY = NONE;

//! Quadratic spacing in the vertical is 4 times finer near the base than equal spacing.
const PetscScalar IceGrid::DEFAULT_QUADZ_LAMBDA  = 4.0;

//! Default computational box is 3000 km x 3000 km (= 2 \c Lx x 2 \c Ly) in horizontal.
const PetscScalar IceGrid::DEFAULT_ICEPARAM_Lx   = 1500.0e3;

//! Default computational box is 3000 km x 3000 km (= 2 \c Lx x 2 \c Ly) in horizontal.
const PetscScalar IceGrid::DEFAULT_ICEPARAM_Ly   = 1500.0e3;

//! Default computational box for ice is 4000 m high.
const PetscScalar IceGrid::DEFAULT_ICEPARAM_Lz   = 4000.0;

//! Start at year zero by default.
const PetscScalar IceGrid::DEFAULT_ICEPARAM_year = 0.0;


//! Default grid is 61 x 61 in horizontal.
const PetscInt    IceGrid::DEFAULT_ICEPARAM_Mx   = 61;

//! Default grid is 61 x 61 in horizontal.
const PetscInt    IceGrid::DEFAULT_ICEPARAM_My   = 61;

//! Default grid has 31 levels in the vertical.
const PetscInt    IceGrid::DEFAULT_ICEPARAM_Mz   = 31;

//! Default grid has no grid for bedrock; top bedrock level duplicates bottom of ice.
const PetscInt    IceGrid::DEFAULT_ICEPARAM_Mbz  = 1;


IceGrid::IceGrid(MPI_Comm c,
                 PetscMPIInt r,
                 PetscMPIInt s):
  com(c), rank(r), size(s) { 

  // The grid in symmetric with respect to zero by default.
  x0 = 0.0;
  y0 = 0.0;

  vertical_spacing = DEFAULT_SPACING_TYPE;
  periodicity = DEFAULT_PERIODICITY;
  Lx = DEFAULT_ICEPARAM_Lx;
  Ly = DEFAULT_ICEPARAM_Ly;
  Lz = DEFAULT_ICEPARAM_Lz;
  year = DEFAULT_ICEPARAM_year;
  Mx = DEFAULT_ICEPARAM_Mx;
  My = DEFAULT_ICEPARAM_My;
  Mz = DEFAULT_ICEPARAM_Mz;
  Mbz = DEFAULT_ICEPARAM_Mbz;

  da2 = PETSC_NULL;
  zlevels = NULL;
  zblevels = NULL;

  compute_vertical_levels();
  compute_horizontal_spacing();
}


IceGrid::~IceGrid() {
  if (da2 != PETSC_NULL) {
    DADestroy(da2);
  }
  delete [] zlevels;
  delete [] zblevels;
}

//! Set the vertical levels according to values in Mz, Lz, Mbz, and the vertical_spacing data member.
/*!
Sets \c dzMIN, \c dzMAX, and \c Lbz.  Sets and re-allocates \c zlevels[] 
and \c zblevels[].

Uses \c Mz, \c Lz, \c Mbz, and \c spacing_type.  Note that \c vertical_spacing
cannot be UNKNOWN.

This procedure is only called when a grid is determined from scratch, e.g. 
by a derived class or when bootstrapping from 2D data only, but not when 
reading a model state input file (which will have its own grid,
which may not even be a grid created by this routine).
  - When \c vertical_spacing == EQUAL, the vertical grid in the ice is equally spaced:
    <tt>zlevels[k] = k dzMIN</tt> where <tt>dzMIN = Lz / (Mz - 1)</tt>.  
    In this case <tt>dzMIN = dzMAX</tt>.
  - When \c vertical_spacing == CHEBYSHEV, the vertical grid in the ice is Chebyshev spaced.  
    Note that (generally speaking) the \f$N+1\f$ \e Chebyshev \e extreme 
    \e points are \f$x_j = \cos(j \pi/N)\f$ for \f$j=0,1,\dots,N\f$. 
    (See \lo\cite{Trefethen}\elo.)
    These are concentrated at either end of the interval \f$[-1,1]\f$.  In our 
    case we want points concentrated near zero, and we use only half of the 
    Chebyshev points because we don't need concentration near the top
    of the computational box.  So we take the original Chebyshev extreme points
    \f$x_j\f$ with \f$N= 2\, \mathtt{Mz} - 1\f$ but we choose only 
    \f$j=0,1,\dots,\mathtt{Mz}-1\f$.  These points satisfy 
    \f$0 \le x_j \le 1\f$ and they are clustered near \f$x=1\f$.  Then we 
    flip and scale: \f$z_j = \mathtt{Lz} (1 - x_j)\f$.  The smallest spacing 
    is a factor proportional to \c Mz smaller than the equal spacing.  That is, 
    \f$z_1 = C \mathtt{Lz} / (\mathtt{Mz}^2)\f$ while 
    \f$dzEQ = \mathtt{Lz}/(\mathtt{Mz}-1)\f$.   Near the top the spacing is, 
    to good approximation, equal to \f$\pi \mathtt{dzEQ}\f$; the actual top 
    space is recorded as \c dzMAX.
  - When \c vertical_spacing == QUADRATIC, the spacing is a quadratic function.  The intent 
    is that the spacing is smaller near the base than near the top, but that 
    the effect is less extreme than the Chebyshev case.  In particular, if 
    \f$\zeta_k = k / (\mathtt{Mz} - 1)\f$ then <tt>zlevels[k] = Lz * 
    ( (\f$\zeta_k\f$ / \f$\lambda\f$) * (1.0 + (\f$\lambda\f$ - 1.0) 
    * \f$\zeta_k\f$) )</tt> where \f$\lambda\f$ = 4.  The value \f$\lambda\f$ 
    indicates the slope of the quadratic function as it leaves the base.  
    Thus a value of \f$\lambda\f$ = 4 makes the spacing about four times finer 
    at the base than equal spacing would be.
      
In all cases the spacing is equal in the bedrock: <tt>zblevels[k] = </tt> where 
<tt>dzEQ = Lz / (Mz - 1)</tt> and <tt>Lbz = dzEQ * (Mbz - 1)</tt>.  That is, 
both the \c zblevels[] and the depth into the bedrock to which the 
computational box extends (\c Lbz) are set according to the value
of \c Mbz and according to the equal spacing which would apply to the ice.
 */
PetscErrorCode  IceGrid::compute_vertical_levels() {
  
  if (vertical_spacing == UNKNOWN) {
    SETERRQ(1,"IceGrid::compute_vertical_levels(): vertical_spacing can not be UNKNOWN");
  }
  
  if (Mz < 2) {
    SETERRQ(2,"IceGrid::compute_vertical_levels(): Mz must be at least 2.");
  }

  if (Mbz < 1) {
    SETERRQ(3, "IceGrid::compute_vertical_levesl(): Mbz must be at least 1.");
  }

  if (Lz <= 0) {
    SETERRQ(4, "IceGrid::compute_vertical_levels(): Lz must be positive.");
  }

  // Fill the levels in the ice:
  delete [] zlevels;
  zlevels = new PetscScalar[Mz];

  switch (vertical_spacing) {
  case EQUAL: {
    dzMIN = Lz / ((PetscScalar) Mz - 1);
    dzMAX = dzMIN;
    
    // Equal spacing
    for (PetscInt k=0; k < Mz - 1; k++) {
      zlevels[k] = dzMIN * ((PetscScalar) k);
    }
    zlevels[Mz - 1] = Lz;  // make sure it is exactly equal
    break;
  }
  case CHEBYSHEV: {
    // Spaced according to the Chebyshev extreme points in the interval [0,1],
    //   with 1 flipped to be the base of the ice, and stretched.
    for (PetscInt k=0; k < Mz - 1; k++) {
      zlevels[k] = Lz * ( 1.0 - cos((pi/2.0) * k / (Mz-1)) );
    }
    zlevels[Mz - 1] = Lz;  // make sure it is exactly equal
    dzMIN = zlevels[1] - zlevels[0];
    dzMAX = zlevels[Mz-1] - zlevels[Mz-2];
    break;
  }
  case QUADRATIC: {
    // this quadratic scheme is an attempt to be less extreme in the fineness near the base.
    const PetscScalar  lam = DEFAULT_QUADZ_LAMBDA;  
    for (PetscInt k=0; k < Mz - 1; k++) {
      const PetscScalar zeta = ((PetscScalar) k) / ((PetscScalar) Mz - 1);
      zlevels[k] = Lz * ( (zeta / lam) * (1.0 + (lam - 1.0) * zeta) );
    }
    zlevels[Mz - 1] = Lz;  // make sure it is exactly equal
    dzMIN = zlevels[1] - zlevels[0];
    dzMAX = zlevels[Mz-1] - zlevels[Mz-2];
    break;
  }
  default:
    SETERRQ(5, "Can't happen.");
  }

  // Fill the bedrock levels:
  delete [] zblevels;
  zblevels = new PetscScalar[Mbz];

  // equally-spaced in bedrock with spacing determined by Lz, Mz and vertical_spacing
  if (Mbz == 1) {
    zblevels[0] = 0.0;
    Lbz = 0.0;
  } else {
    const PetscScalar dzEQ = dzMIN;
    Lbz = dzEQ * ((PetscScalar) Mbz - 1);
    for (PetscInt kb=0; kb < Mbz; kb++) {
      zblevels[kb] = -Lbz + dzEQ * ((PetscScalar) kb);
    }
  }

  return 0;
}


//! Print the vertical levels in \c zlevels[] and \c zblevels[] to standard out.
PetscErrorCode IceGrid::printVertLevels(int verbosity) {
  PetscErrorCode ierr;
  ierr = verbPrintf(verbosity,com,
     "    vertical levels in ice (Mz=%d,Lz=%5.4f): ",Mz,Lz); CHKERRQ(ierr);
  for (PetscInt k=0; k < Mz; k++) {
    ierr = verbPrintf(verbosity,com," %5.4f,",zlevels[k]); CHKERRQ(ierr);
  }
  ierr = verbPrintf(verbosity,com,
     "\n    vertical levels in bedrock (Mbz=%d,Lbz=%5.4f): ",Mbz,Lbz); CHKERRQ(ierr);
  for (PetscInt kb=0; kb < Mbz; kb++) {
    ierr = verbPrintf(verbosity,com," %5.4f,",zblevels[kb]); CHKERRQ(ierr);
  }
  ierr = verbPrintf(verbosity,com,"\n"); CHKERRQ(ierr);
  return 0;
}

//! Return the index \c k into \c zlevels[] so that <tt>zlevels[k] <= height < zlevels[k+1]</tt> and <tt>k < Mz</tt>.
PetscInt IceGrid::kBelowHeight(PetscScalar height) {
  if (height < 0.0 - 1.0e-6) {
    PetscPrintf(com, 
       "IceGrid kBelowHeight(): height = %5.4f is below base of ice (height must be non-negative)\n",
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
//  while ((zlevels[mcurr+1] <= height) && (mcurr+1 < Mz)) {
  while (zlevels[mcurr+1] < height) {
    mcurr++;
  }
  return mcurr;
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
    vertical_spacing = EQUAL;
  } else {
    vertical_spacing = UNKNOWN;
  }
  return 0;
}


//! Create the PETSc DA \c da2 for the horizontal grid.  Determine how the horizontal grid is divided among processors.
/*!
  This procedure should only be called after the parameters describing the
  horizontal computational box (Lx,Ly) and the parameters for the horizontal
  grid (Mx,My) are already determined. In particular, the input file (either \c
  -i or \c -boot_from) and user options (like \c -Mx) must have already been
  read to determine the parameters, and any conflicts must have been resolved.
 */
PetscErrorCode IceGrid::createDA() {
  PetscErrorCode ierr;

  if (da2 != PETSC_NULL) {
    SETERRQ(1, "IceGrid::createDA(): da2 != PETSC_NULL");
  }

  // this line contains the fundamental transpose: My,Mx instead of "Mx,My";
  // this transpose allows us to index the Vecs by "[i][j]" but display them
  //   correctly in PETSc viewers; changing this to untransposed should *only*
  //   change that viewer behavior
  ierr = DACreate2d(com, DA_XYPERIODIC, DA_STENCIL_BOX,
                    My, Mx, PETSC_DECIDE, PETSC_DECIDE, 1, 1,
                    PETSC_NULL, PETSC_NULL, &da2); CHKERRQ(ierr);

  /* see IceModelVec3 for creation of three-dimensional DAs for each 3D Vec */

  DALocalInfo info;
  ierr = DAGetLocalInfo(da2, &info); CHKERRQ(ierr);
  // this continues the fundamental transpose
  xs = info.ys; xm = info.ym;
  ys = info.xs; ym = info.xm;

  return 0;
}

//! Sets grid vertical levels, Mz, Mbz, Lz and Lbz. Checks input for consistency.
PetscErrorCode IceGrid::set_vertical_levels(int new_Mz, int new_Mbz, double *new_zlevels, double *new_zblevels) {
  if (new_Mz < 2) {
    SETERRQ(1, "IceGrid::set_vertical_levels(): Mz has to be at least 2.");
  }
  if (new_Mbz < 1) {
    SETERRQ(2, "IceGrid::set_vertical_levels(): Mbz has to be at least 1.");
  }

  if ( (!isIncreasing(new_Mz, new_zlevels)) || (PetscAbs(new_zlevels[0]) > 1.0e-10) ) {
    SETERRQ(3, "IceGrid::set_vertical_levels(): invalid zlevels; must be strictly increasing and start with z=0.");
  }

  if ( (!isIncreasing(new_Mbz, new_zblevels)) || (PetscAbs(new_zblevels[new_Mbz-1]) > 1.0e-10) ) {
    SETERRQ(3, "rescale: zblevels invalid; must be strictly increasing and end with z=0\n");
  }

  Mz  =  new_Mz;
  Mbz =  new_Mbz;
  Lz  =  new_zlevels[Mz - 1];
  Lbz = -new_zblevels[0];

  // Fill the levels in the ice:
  delete[] zlevels;
  zlevels = new PetscScalar[Mz];
  for (int j = 0; j < Mz; j++)
    zlevels[j] = (PetscScalar) new_zlevels[j];
  zlevels[0] = 0.0;		// make sure zlevels start with zero

  // Fill the bedrock levels:
  delete[] zblevels;
  zblevels = new PetscScalar[Mbz];
  for (int j = 0; j < Mbz; j++)
    zblevels[j] = (PetscScalar) new_zblevels[j];
  zblevels[Mbz-1] = 0.0;	// make sure zblevels end with zero

  get_dzMIN_dzMAX_spacingtype();

  return 0;
}

//! Compute horizontal spacing parameters \c dx and \c dy using \c Mx, \c My, \c Lx, \c Ly and periodicity.
/*! 
The grid used in PISM, in particular the PETSc DAs used here, are periodic in x and y.
This means that the ghosted values <tt> foo[i+1][j], foo[i-1][j], foo[i][j+1], foo[i][j-1]</tt>
for all 2D Vecs, and similarly in the x and y directions for 3D Vecs, are always available.
That is, they are available even if i,j is a point at the edge of the grid.  On the other
hand, by default, \c dx  is the full width  <tt>2 * Lx</tt>  divided by  <tt>Mx - 1</tt>.
This means that we conceive of the computational domain as starting at the <tt>i = 0</tt>
grid location and ending at the  <tt>i = Mx - 1</tt>  grid location, in particular.  
This idea is not quite compatible with the periodic nature of the grid.

The upshot is that if one computes in a truly periodic way then the gap between the  
<tt>i = 0</tt>  and  <tt>i = Mx - 1</tt>  grid points should \em also have width  \c dx.  
Thus we compute  <tt>dx = 2 * Lx / Mx</tt>.
 */
PetscErrorCode IceGrid::compute_horizontal_spacing() {
  if (Mx < 3) {
    SETERRQ(1, "IceGrid::set_horizontal_dims(): Mx has to be at least 3.");
  }

  if (My < 3) {
    SETERRQ(2, "IceGrid::set_horizontal_dims(): My has to be at least 3.");
  }

  if (Lx <= 0) {
    SETERRQ(3, "IceGrid::set_horizontal_dims(): Lx has to be positive.");
  }

  if (Ly <= 0) {
    SETERRQ(3, "IceGrid::set_horizontal_dims(): Ly has to be positive.");
  }

  if (periodicity & X_PERIODIC) {
    dx = 2.0 * Lx / Mx;
  } else {
    dx = 2.0 * Lx / (Mx - 1);
  }

  if (periodicity & Y_PERIODIC) {
    dy = 2.0 * Ly / My;
  } else {
    dy = 2.0 * Ly / (My - 1);
  }

  return 0;
}

//! \brief Computes values of x and y corresponding to the computational grid,
//! with accounting for periodicity.
/*! This method allocates arrays \c x and \c y, and they have to be freed
  (using delete[]) by the caller.
 */
PetscErrorCode IceGrid::compute_horizontal_coordinates(double* &x, double* &y) {
  PetscErrorCode ierr;

  // make sure that dx and dy are set correctly:
  ierr = compute_horizontal_spacing(); CHKERRQ(ierr);

  x = new double[Mx];
  y = new double[My];

  double x_min = -Lx + x0,
    x_max = Lx + x0,
    y_min = -Ly + y0,
    y_max = Ly + y0;

  if (periodicity & X_PERIODIC) {
    x_max -= dx;
  }

  if (periodicity & Y_PERIODIC) {
    y_max -= dy;
  }

  for (int i = 0; i < Mx; ++i)
    x[i] = x_min + i * dx;
  x[Mx - 1] = x_max;

  for (int i = 0; i < My; ++i)
    y[i] = y_min + i * dy;
  y[My - 1] = y_max;

  return 0;
}
