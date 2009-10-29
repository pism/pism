// Copyright (C) 2007--2009 Ed Bueler
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

#ifndef __deformation_hh
#define __deformation_hh

#include "NCVariable.hh"
#include <petscvec.h>
#if (PISM_HAVE_FFTW)
#include <fftw3.h>
#endif

//! Class implementing the bed deformation model described in \ref BLKfastearth.
/*!
  This class implements the \ref LingleClark bed deformation model by a Fourier 
  spectral collocation method, as described in \ref BLKfastearth.  (The former
  reference is where the continuum model arose, and a flow-line application is given.
  The latter reference describes a new, fast method and gives verification results.
  See also \ref BLK2006earth if more technical detail and/or Matlab programs are desired.)
  
  Both a viscous half-space model (with elastic
  lithosphere) and a spherical elastic model are computed.  They are superposed
  because the underlying earth model is linear.

  The class assumes that the supplied Petsc Vecs are *sequential*.  It is expected to be 
  run only on processor zero (or possibly by each processor once each processor 
  owns the entire 2D gridded ice thicknesses and bed elevations.)
  
  This class SHOULD!
  include the scatter structures necessary to make this work in parallel.

  A test program for this class is pism/src/verif/tryLCbd.cc.
  
  Called by IceModel in pism/src/iMbeddef.cc.  
 */
class BedDeformLC {
public:
  BedDeformLC();
  ~BedDeformLC();
  PetscErrorCode settings(const NCConfigVariable &config,
			  PetscTruth  myinclude_elastic,
                          PetscInt myMx, PetscInt myMy, PetscScalar mydx, PetscScalar mydy,
                          PetscInt myZ, PetscScalar myicerho,
                          PetscScalar myrho, PetscScalar myeta, PetscScalar myD,
                          Vec* myHstart, Vec* mybedstart, Vec* myuplift,  // initial state
                          Vec* myH,     // generally gets changed by calling program
                                        // before each call to step
                          Vec* mybed);  // mybed gets modified by step()
  PetscErrorCode alloc();
  PetscErrorCode uplift_init();
  PetscErrorCode step(const PetscScalar dtyear, const PetscScalar yearFromStart);

protected:
  PetscTruth    include_elastic;
  PetscInt      Mx, My;
  PetscScalar   dx, dy;
  PetscInt      Z;       // factor by which fat FFT domain is larger than 
                         // region of physical interest
  PetscScalar   icerho,  // ice density (for computing load from volume)
                rho,     // earth density
                eta,     // mantle viscosity
                D;       // lithosphere flexural rigidity
  
private:
  PetscScalar   standard_gravity;
  PetscTruth    settingsDone, allocDone;
  PetscInt      Nx, Ny,      // fat sizes
                Nxge, Nyge;  // fat with boundary sizes
  PetscInt      i0_plate,  j0_plate; // indices into fat array for corner of thin
  PetscScalar   Lx, Ly;      // half-lengths of (thin) physical domain
  PetscScalar   fatLx, fatLy;      // half-lengths of fat computational domain
  PetscScalar   *cx, *cy;      // coeffs of derivs in Fourier space
  Vec           *H, *bed, *Hstart, *bedstart, *uplift; // pointers to sequential
  Vec           Hdiff, dbedElastic,    // sequential; working space
                platefat, plateoffset, // seq and fat
                vleft, vright,  // coeffs; sequential and fat
                lrmE;           // load response matrix (elastic); sequential and fat *with* boundary
#if (PISM_HAVE_FFTW)
  fftw_complex  *bdin, *bdout;  // 2D sequential
  fftw_plan     bdplanfor,bdplanback;
#endif
};

#endif	/* __deformation_hh */

