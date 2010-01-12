// Copyright (C) 2004-2010 Jed Brown, Ed Bueler and Constantine Khroulev
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
const SpacingType IceGrid::DEFAULT_ICE_SPACING_TYPE = QUADRATIC;

const SpacingType IceGrid::DEFAULT_BED_SPACING_TYPE = QUADRATIC;

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

//! Default thickness of the bedrock layer.
const PetscScalar IceGrid::DEFAULT_ICEPARAM_Lbz   = 0.0;

//! Start at year zero by default.
const PetscScalar IceGrid::DEFAULT_ICEPARAM_start_year = 0.0;

//! Run for 100 years by default
const PetscScalar IceGrid::DEFAULT_ICEPARAM_run_length = 1000.0;

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

  ice_vertical_spacing = DEFAULT_ICE_SPACING_TYPE;
  bed_vertical_spacing = DEFAULT_BED_SPACING_TYPE;
  periodicity = DEFAULT_PERIODICITY;
  Lx = DEFAULT_ICEPARAM_Lx;
  Ly = DEFAULT_ICEPARAM_Ly;
  Lz = DEFAULT_ICEPARAM_Lz;
  Lbz = DEFAULT_ICEPARAM_Lbz;
  start_year = DEFAULT_ICEPARAM_start_year;
  year = start_year;
  end_year = start_year + DEFAULT_ICEPARAM_run_length;
  Mx = DEFAULT_ICEPARAM_Mx;
  My = DEFAULT_ICEPARAM_My;
  Mz = DEFAULT_ICEPARAM_Mz;
  Mbz = DEFAULT_ICEPARAM_Mbz;

  initial_Mz = 0;		// will be set to a correct value in
				// IceModel::check_maximum_thickness()

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


//! Compute vertical levels in both ice and bedrock.
PetscErrorCode IceGrid::compute_vertical_levels() {
  PetscErrorCode ierr;

  ierr = compute_ice_vertical_levels(); CHKERRQ(ierr);
  ierr = compute_bed_vertical_levels(); CHKERRQ(ierr);

  return 0;
}

//! \brief Set the vertical levels in the ice according to values in Mz, Lz,
//! Mbz, and the ice_vertical_spacing data member.
/*!
Sets \c dzMIN, \c dzMAX, and \c Lbz.  Sets and re-allocates \c zlevels[] 
and \c zblevels[].

Uses \c Mz, \c Lz, \c Mbz, and \c spacing_type.  Note that \c ice_vertical_spacing
cannot be UNKNOWN.

This procedure is only called when a grid is determined from scratch, %e.g. 
by a derived class or when bootstrapping from 2D data only, but not when 
reading a model state input file (which will have its own grid,
which may not even be a grid created by this routine).
  - When \c vertical_spacing == EQUAL, the vertical grid in the ice is equally spaced:
    <tt>zlevels[k] = k dzMIN</tt> where <tt>dzMIN = Lz / (Mz - 1)</tt>.  
    In this case <tt>dzMIN = dzMAX</tt>.
  - When \c vertical_spacing == QUADRATIC, the spacing is a quadratic function.  The intent 
    is that the spacing is smaller near the base than near the top.  In particular, if 
    \f$\zeta_k = k / (\mathtt{Mz} - 1)\f$ then <tt>zlevels[k] = Lz * 
    ( (\f$\zeta_k\f$ / \f$\lambda\f$) * (1.0 + (\f$\lambda\f$ - 1.0) 
    * \f$\zeta_k\f$) )</tt> where \f$\lambda\f$ = 4.  The value \f$\lambda\f$ 
    indicates the slope of the quadratic function as it leaves the base.  
    Thus a value of \f$\lambda\f$ = 4 makes the spacing about four times finer 
    at the base than equal spacing would be.
 */
PetscErrorCode  IceGrid::compute_ice_vertical_levels() {
  
  if (Mz < 2) {
    SETERRQ(2,"IceGrid::compute_ice_vertical_levels(): Mz must be at least 2.");
  }

  if (Lz <= 0) {
    SETERRQ(4, "IceGrid::compute_ice_vertical_levels(): Lz must be positive.");
  }

  // Fill the levels in the ice:
  delete [] zlevels;
  zlevels = new PetscScalar[Mz];

  switch (ice_vertical_spacing) {
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
    SETERRQ(1,"IceGrid::compute_ice_vertical_levels(): ice_vertical_spacing can not be UNKNOWN.");
  }

  return 0;
}

//! Compute vertical levels in the bedrock.
/*!
  PISM supports equal and quadratic spacing of bedrock levels.

  
  
  Please see IceGrid::compute_ice_vertical_levels() for more.
 */
PetscErrorCode IceGrid::compute_bed_vertical_levels() {

  if (Mbz < 1) {
    SETERRQ(3, "IceGrid::compute_bed_vertical_levels(): Mbz must be at least 1.");
  }

  if (Lbz < 0) {
    SETERRQ(4, "IceGrid::compute_bed_vertical_levels(): Lbz must be zero or positive.");
  }

  if ((Lbz < 1e-9) && (Mbz > 1)) {
    SETERRQ(5, "IceGrid::compute_bed_vertical_levels(): Lbz must be positive if Mbz > 1.");
  }

  // Fill the bedrock levels:
  delete [] zblevels;
  zblevels = new PetscScalar[Mbz];

  if (Mbz == 1) {		// we have no bedrock
    zblevels[0] = 0.0;
    Lbz = 0.0;
    dzbMIN = dzbMAX = dzMIN;
  } else {			// we have bedrock
    switch (bed_vertical_spacing) {
    case EQUAL:
      {
	// equally-spaced in bedrock with spacing determined by Lbz and Mbz
	dzbMIN = dzbMAX = Lbz / (Mbz - 1);

	zblevels[0] = -Lbz;
	for (PetscInt kb=1; kb < Mbz - 1; kb++)
	  zblevels[kb] = -Lbz + dzbMIN * ((PetscScalar) kb);
	zblevels[Mbz - 1] = 0.0;
	break;
      }
    case QUADRATIC:
      {
	// this is the lambda we have in bedrock; "2" below means that bedrock
	// level thickness should increase with depth at about twice the rate of
	// increase of ice levels.
	PetscScalar lam = 2*DEFAULT_QUADZ_LAMBDA;

	zblevels[0] = -Lbz;	// make sure it's right on
	for (PetscInt k = 1; k < Mbz - 1; k++) {
	  const PetscScalar zeta = ((PetscScalar) Mbz - 1 - k) / ((PetscScalar) Mbz - 1);
	  zblevels[k] = -Lbz * (zeta / lam) * (1.0 + (lam - 1.0) * zeta);
	}
	zblevels[Mbz - 1] = 0;	// make sure this one is right on, too

	dzbMAX = zblevels[1] - zblevels[0];
	dzbMIN = zblevels[Mbz-1] - zblevels[Mbz-2];
	break;
      }
    default:
      SETERRQ(1, "PISM ERROR: Only equal and quadratic bedrock vertical spacing types are supported.");
    }
  } // end of } else {...
  
  return 0;
}


//! Print to stdout information on computational domain and grid (other than the vertical levels themselves).
PetscErrorCode IceGrid::printInfo(const int verbosity) {
  PetscErrorCode ierr;
  ierr = verbPrintf(verbosity,com,
         "  IceGrid parameters:\n"); CHKERRQ(ierr);
  ierr = verbPrintf(verbosity,com,
         "            Lx = %6.2f km, Ly = %6.2f km, Lz = %6.2f m, Lbz = %6.2f m,\n",
         Lx/1000.0,Ly/1000.0,Lz,Lbz); CHKERRQ(ierr);
  ierr = verbPrintf(verbosity,com,
         "            x0 = %6.2f km, y0 = %6.2f km,   (coordinates of center)\n",
		    x0/1000.0, y0/1000.0); CHKERRQ(ierr);
  ierr = verbPrintf(verbosity,com,
         "            Mx = %d, My = %d, Mz = %d, Mbz = %d,\n",
         Mx,My,Mz,Mbz); CHKERRQ(ierr);
  ierr = verbPrintf(verbosity,com,
         "            dx = %6.3f km, dy = %6.3f km, year = %8.4f]\n",
         dx/1000.0,dy/1000.0,year); CHKERRQ(ierr);
  return 0;
}


//! Print the vertical levels in \c zlevels[] and \c zblevels[] to stdout.
PetscErrorCode IceGrid::printVertLevels(const int verbosity) {
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

//! \brief From given vertical grid zlevels[], determine \c dzMIN, \c dzMAX, \c dzbMIN, \c dzbMAX and
//! determine whether ice and bedrock vertical spacings are equal.
/*! The standard for equal vertical spacing in the ice is \f$10^{-8}\f$ m max
  difference between \c dzMIN and \c dzMAX. Similar for the bedrock.
 */
PetscErrorCode IceGrid::get_dzMIN_dzMAX_spacingtype() {

  // ice:
  dzMIN = Lz; 
  dzMAX = 0.0;
  for (PetscInt k = 0; k < Mz - 1; k++) {
    const PetscScalar mydz = zlevels[k+1] - zlevels[k];
    dzMIN = PetscMin(mydz,dzMIN);
    dzMAX = PetscMax(mydz,dzMAX);
  }
  if (PetscAbs(dzMAX - dzMIN) <= 1.0e-8) {
    ice_vertical_spacing = EQUAL;
  } else {
    ice_vertical_spacing = UNKNOWN;
  }

  // bedrock:
  if (Mbz == 1) {
    dzbMIN = dzbMAX = dzMIN;
  } else {
    dzbMIN = Lbz; 
    dzbMAX = 0.0;
    for (PetscInt k = 0; k < Mbz - 1; k++) {
      const PetscScalar mydz = zblevels[k+1] - zblevels[k];
      dzbMIN = PetscMin(mydz,dzMIN);
      dzbMAX = PetscMax(mydz,dzMAX);
    }
  }

  if (PetscAbs(dzMAX - dzMIN) <= 1.0e-8) {
    bed_vertical_spacing = EQUAL;
  } else {
    bed_vertical_spacing = UNKNOWN;
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

  This method contains the "fundamental" transpose: "My,Mx" instead of "Mx,My"
  in the DACreate2d call; this transpose allows us to index arrays by "[i][j]"
  (where 'i' corresponds to 'x' and 'j' to 'y') and be consistent about
  meanings of 'x', 'y', 'u' and 'v'.

  Unfortunately this means that PETSc viewers appear transposed.

  This choice should be virtually invisible, unless you're using DALocalInfo
  structures.

  See IceModelVec3 for creation of three-dimensional DAs.

  \note PETSc order: x in columns, y in rows, indexing as array[y][x]. PISM
  order: x in rows, y in columns, indexing as array[x][y].
 */
PetscErrorCode IceGrid::createDA() {
  PetscErrorCode ierr;
  if (da2 != PETSC_NULL)
    SETERRQ(1, "IceGrid::createDA(): da2 != PETSC_NULL");

  // Transpose:
  ierr = DACreate2d(com, DA_XYPERIODIC, DA_STENCIL_BOX,
                    My, Mx, PETSC_DECIDE, PETSC_DECIDE, 1, 1,
                    PETSC_NULL, PETSC_NULL, &da2); CHKERRQ(ierr);

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

  if ( (!is_increasing(new_Mz, new_zlevels)) || (PetscAbs(new_zlevels[0]) > 1.0e-10) ) {
    SETERRQ(3, "IceGrid::set_vertical_levels(): invalid zlevels; must be strictly increasing and start with z=0.");
  }

  if ( (!is_increasing(new_Mbz, new_zblevels)) || (PetscAbs(new_zblevels[new_Mbz-1]) > 1.0e-10) ) {
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

//! Computes fine vertical spacing in the ice and bedrock.
/*! The computations in IceModel::temperatureStep() and IceModel::ageStep() use
  a fine equally-spaced grid.

  This method computes the number of levels and the levels themselves so that
  fine equally-spaced grids in ice and bedrock have the same spacing (if
  bedrock is present).

  Arrays fzlev and fzblev are allocated here and have to be freed by the
  caller.

  Mapping to and from the storage grid occurs in IceModelVec3 and
  IceModelVec3Bedrock methods.

  Note that the computational grid in the ice is allowed to exceed Lz; we need
  this to match spacings (see above). The temperature field is extrapolated to
  the extra level using the value from the topmost level.
 */
PetscErrorCode IceGrid::get_fine_vertical_grid(PetscInt &fMz, PetscInt &fMbz,
					       PetscScalar &fdz, PetscScalar &fdzb,
					       PetscScalar* &fzlev, PetscScalar* &fzblev) {

  // the smallest of the spacings used in ice and bedrock:
  PetscScalar dz_fine = PetscMin(dzMIN, dzbMIN);

  if (Lbz > 1e-9) {		// we have bedrock

    // the number of levels of the computational grid in bedrock:
    fMbz = static_cast<PetscInt>(ceil(Lbz / dz_fine) + 1);

    // recompute fine spacing so that Lbz is an integer multiple of dz_fine:
    dz_fine = Lbz / (fMbz - 1);

    // number of levels in the ice; this is one level too many, but spacing
    // matches the one used in the bedrock and the extra level is in the air,
    // anyway
    fMz = static_cast<PetscInt>(ceil(Lz / dz_fine) + 2);
  } else {			// we don't have bedrock

    fMbz = 1;
    fMz = static_cast<PetscInt>(ceil(Lz / dz_fine) + 1);
    dz_fine = Lz / (fMz - 1);
  }

  // both ice and bedrock will have this spacing
  fdz = fdzb = dz_fine;

  // allocate arrays:
  fzlev  = new PetscScalar[fMz];
  fzblev = new PetscScalar[fMbz];

  // compute levels in the ice:
  for (PetscInt k = 0; k < fMz; k++)
    fzlev[k] = ((PetscScalar) k) * dz_fine;
  // Note that it's allowed to go over Lz.

  // and bedrock:
  for (PetscInt kb = 0; kb < fMbz-1; kb++)
    fzblev[kb] = - Lbz + dz_fine * ((PetscScalar) kb);
  fzblev[fMbz-1] = 0.0;  // make sure it is right on
  
  return 0;
}

//! Computes fine vertical spacing in the ice.
/*! Similar to IceGrid::get_fine_vertical_grid(), but used in the age
  computation. The grid created by this method does not extend into the
  bedrock, so neither matching the spacing nor extending past Lz is necessary.
 */
PetscErrorCode IceGrid::get_fine_vertical_grid_ice(PetscInt &fMz, PetscScalar &fdz,
						   PetscScalar* &fzlev) {
  if (ice_vertical_spacing == EQUAL) {
    // just use the storage grid
    fdz = dzMIN;
    fMz = Mz;
    fzlev = new PetscScalar[fMz];
    for (PetscInt k = 0; k < fMz; k++)
      fzlev[k] = zlevels[k];
  } else {
    PetscScalar dz_fine = dzMIN;
    fMz = static_cast<PetscInt>(ceil(Lz / dz_fine) + 1);
    dz_fine = Lz / (fMz - 1);

    fdz = dz_fine;
    fzlev  = new PetscScalar[fMz];

    // compute levels in the ice:
    for (PetscInt k = 0; k < fMz - 1; k++)
      fzlev[k] = ((PetscScalar) k) * dz_fine;
    fzlev[fMz - 1] = Lz;
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
