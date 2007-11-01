// Copyright (C) 2007 Ed Bueler
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

#ifndef __beddefLC_hh
#define __beddefLC_hh

#include <petscvec.h>
#if (WITH_FFTW)
#include <fftw3.h>
#endif

//--------
// see Bueler et al (2007) and Lingle & Clark (1985) for G^E(r), the Green's function
// of spherical layered elastic earth model

struct ge_params {
   double dx, dy;
   int p, q; 
};

double ge_integrand(unsigned ndimMUSTBETWO, const double* xiANDeta, void* paramsIN);
//---------

//--------
// for the solution of the disc load case of the viscous half-space model, see 
// appendix B of  E. Bueler, C. S. Lingle, and J. A. Kallen-Brown (2006)
// "Computation of a combined spherical-elastic and viscous-half-space Earth model
// for ice sheet simulation", arXiv:physics/0606074

struct vd_params {
   double t, R0, rk, rho, grav, D, eta;
};

double viscDiscIntegrand (double kap, void * paramsIN);

double viscDisc(double t, double H0, double R0, double r, 
                double rho, double grav, double D, double eta);
//--------


//--------
// this works on sequential Vecs; it is directly analogous to the 
// Matlab command "conv2(A,B,'same')"
PetscErrorCode conv2_same(Vec vA, const PetscInt mA, const PetscInt nA, 
                          Vec vB, const PetscInt mB, const PetscInt nB,
                          Vec &vresult);
//--------


//--------
/*
  This class implements the Lingle & Clark bed deformation model by a Fourier 
  spectral collocation method.  Both a viscous half-space model (with elastic
  lithosphere) and a spherical elastic model are computed.  They are superposed
  because the underlying earth model is linear.

  See
       E. Bueler, C. S. Lingle, and J. Brown (2007) "Fast computation of a viscoelastic
       deformable Earth model for ice sheet simulations", Ann. Glaciol. 46, 97--105
  for the method and verification results.  The continuum model arose in
       C. S. Lingle and J. A. Clark (1985) "A numerical model of interactions between
       a marine ice sheet and the solid earth: {A}pplication to a {W}est {A}ntarctic
       ice stream", J. Geophys. Res. 90 (C1), 1100--1114
  to model the deformation of the earth under an ice stream (flow line model), though
  different numerics are used there.
  
  See also the preprint Bueler, Lingle, Kallen-Brown (2006) if more technical detail 
  or Matlab programs are desired.

  The model can initialize from an uplift map as described in Bueler et al (2007).

  The class assumes that the supplied Petsc Vecs are *sequential*.  It is expected to be 
  run only on processor zero (or possibly by each processor once each processor 
  owns the entire 2D gridded ice thicknesses and bed elevations.)  This class does not
  include the scatter structures necessary to make this work in parallel;
  see pism/src/iMbeddef.cc.
  
  A test program for this class is pism/src/exact/tryLCbd.cc.
*/
//--------

class BedDeformLC {
public:
  BedDeformLC();
  ~BedDeformLC();
  PetscErrorCode settings(PetscTruth  myinclude_elastic,
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
#if (WITH_FFTW)
  fftw_complex  *bdin, *bdout;  // 2D sequential
  fftw_plan     bdplanfor,bdplanback;
#endif
};

#endif	/* __beddefLC_hh */

