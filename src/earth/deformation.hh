// Copyright (C) 2007--2009, 2011, 2012 Ed Bueler
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
#if (PISM_USE_FFTW==1)
#include <fftw3.h>
#endif

//! Class implementing the bed deformation model described in [\ref BLKfastearth].
/*!
  This class implements the [\ref LingleClark] bed deformation model by a Fourier 
  spectral collocation method, as described in [\ref BLKfastearth].  (The former
  reference is where the continuum model arose, and a flow-line application is given.
  The latter reference describes a new, fast method and gives verification results.
  See also [\ref BLK2006earth] if more technical detail and/or Matlab programs are desired.)

  Both a viscous half-space model (with elastic
  lithosphere) and a spherical elastic model are computed.  They are superposed
  because the underlying earth model is linear.

  The class assumes that the supplied Petsc Vecs are *sequential*.  It is expected to be 
  run only on processor zero (or possibly by each processor once each processor 
  owns the entire 2D gridded ice thicknesses and bed elevations.)

  This class SHOULD!
  include the scatter structures necessary to make this work in parallel.

  A test program for this class is pism/src/verif/tryLCbd.cc.
 */
class BedDeformLC {
public:
  BedDeformLC();
  ~BedDeformLC();
  PetscErrorCode settings(const NCConfigVariable &config,
			  PetscBool myinclude_elastic,
                          PetscInt myMx, PetscInt myMy, PetscScalar mydx, PetscScalar mydy,
                          PetscInt myZ,
                          Vec* myHstart, Vec* mybedstart, Vec* myuplift,  // initial state
                          Vec* myH,     // generally gets changed by calling program
                                        // before each call to step
                          Vec* mybed);  // mybed gets modified by step()
  PetscErrorCode alloc();
  PetscErrorCode uplift_init();
  PetscErrorCode step(const PetscScalar dtyear, const PetscScalar yearFromStart);

protected:
  PetscBool     include_elastic;
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
  PetscBool    settingsDone, allocDone;
  PetscInt      Nx, Ny,      // fat sizes
                Nxge, Nyge;  // fat with boundary sizes
  PetscInt      i0_plate,  j0_plate; // indices into fat array for corner of thin
  PetscScalar   Lx, Ly;      // half-lengths of the physical domain
  PetscScalar   Lx_fat, Ly_fat; // half-lengths of the FFT (spectral) computational domain
  PetscScalar   *cx, *cy;      // coeffs of derivatives in Fourier space
  Vec           *H, *bed, *H_start, *bed_start, *uplift; // pointers to sequential
  Vec           Hdiff, dbedElastic,    // sequential; working space
                U, U_start, // sequential and fat
                vleft, vright,  // coefficients; sequential and fat
                lrmE;           // load response matrix (elastic); sequential and fat *with* boundary
  fftw_complex  *fftw_input, *fftw_output, *loadhat;  // 2D sequential
  fftw_plan     dft_forward, dft_inverse;

  void tweak(PetscReal seconds_from_start);

  void clear_fftw_input();
  void copy_fftw_output(fftw_complex *buffer);
  void set_fftw_input(Vec input, PetscReal normalization, int M, int N, int i0, int j0);
  void get_fftw_output(Vec output, PetscReal normalization, int M, int N, int i0, int j0);
};

class PetscVecAccessor2D {
public:
  PetscVecAccessor2D(Vec vec, PetscInt my_Mx, PetscInt my_My)
    : Mx(my_Mx), My(my_My), i_offset(0), j_offset(0), v(vec)
  { VecGetArray2d(v, Mx, My, 0, 0, &array); }

  PetscVecAccessor2D(Vec vec, PetscInt my_Mx, PetscInt my_My, PetscInt i0, PetscInt j0)
    : Mx(my_Mx), My(my_My), i_offset(i0), j_offset(j0), v(vec)
  { VecGetArray2d(v, Mx, My, 0, 0, &array); }

  ~PetscVecAccessor2D()
  { VecRestoreArray2d(v, Mx, My, 0, 0, &array); }

  inline PetscScalar& operator()(int i, int j)
  { return array[i + i_offset][j + j_offset]; }
private:
  PetscInt Mx, My, i_offset, j_offset;
  Vec v;
  PetscScalar **array;
};

template <class T>
class VecAccessor2D {
public:
  VecAccessor2D(T* a, int my_Mx, int my_My, int my_i_offset, int my_j_offset)
    : Mx(my_Mx), My(my_My), i_offset(my_i_offset), j_offset(my_j_offset), array(a) {}

  VecAccessor2D(T* a, int my_Mx, int my_My)
    : Mx(my_Mx), My(my_My), i_offset(0), j_offset(0), array(a) {}

  inline T& operator()(int i, int j)
  { return array[(j_offset + j) + My * (i_offset + i)]; }

private:
  int Mx, My, i_offset, j_offset;
  T* array;
};

#endif	/* __deformation_hh */

