// Copyright (C) 2004-2009, 2011 Ed Bueler
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

#include <cmath>
#include <petscvec.h>
#include <fftw3.h>
#include "pism_const.hh"
#include "matlablike.hh"
#include "greens.hh"
#include "deformation.hh"

BedDeformLC::BedDeformLC() {
  settingsDone = PETSC_FALSE;
  allocDone = PETSC_FALSE;
}

BedDeformLC::~BedDeformLC() {
  if (allocDone == PETSC_TRUE) {

    fftw_destroy_plan(dft_forward);
    fftw_destroy_plan(dft_inverse);
    fftw_free(fftw_input);
    fftw_free(fftw_output);
    fftw_free(loadhat);

    VecDestroy(&Hdiff);
    VecDestroy(&dbedElastic);
    VecDestroy(&U);
    VecDestroy(&U_start);
    VecDestroy(&vleft);
    VecDestroy(&vright);
    VecDestroy(&lrmE);
    delete [] cx;  delete [] cy;
  }
  allocDone = PETSC_FALSE;
}

PetscErrorCode BedDeformLC::settings(const NCConfigVariable &config,
                                     PetscBool  myinclude_elastic,
                                     PetscInt myMx, PetscInt myMy,
                                     PetscScalar mydx, PetscScalar mydy,
                                     PetscInt myZ,
                                     Vec* myHstart, Vec* mybedstart, Vec* myuplift,
                                     Vec* myH, Vec* mybed) {

  // set parameters
  include_elastic = myinclude_elastic;
  Mx = myMx;
  My = myMy;
  dx = mydx;
  dy = mydy;
  Z = myZ;
  icerho = config.get("ice_density");
  rho = config.get("lithosphere_density");
  eta = config.get("mantle_viscosity");
  D = config.get("lithosphere_flexural_rigidity");

  standard_gravity = config.get("standard_gravity");

  // derive more parameters
  Lx = ((Mx - 1) / 2) * dx;
  Ly = ((My - 1) / 2) * dy;
  Nx = Z*(Mx - 1);
  Ny = Z*(My - 1);
  Lx_fat = (Nx / 2) * dx;
  Ly_fat = (Ny / 2) * dy;
  Nxge = Nx + 1;
  Nyge = Ny + 1;
  i0_plate = (Z - 1)*(Mx - 1) / 2;
  j0_plate = (Z - 1)*(My - 1) / 2;

  // attach to existing (must be allocated!) sequential Vecs
  H = myH;
  bed = mybed;
  H_start = myHstart;
  bed_start = mybedstart;
  uplift = myuplift;

  settingsDone = PETSC_TRUE;
  return 0;
}


PetscErrorCode BedDeformLC::alloc() {
  PetscErrorCode  ierr;
  if (settingsDone == PETSC_FALSE) {
    SETERRQ(PETSC_COMM_SELF, 1, "BedDeformLC must be set with settings() before alloc()\n");
  }
  if (allocDone == PETSC_TRUE) {
    SETERRQ(PETSC_COMM_SELF, 2, "BedDeformLC already allocated\n");
  }

  ierr = VecDuplicate(*H, &Hdiff); CHKERRQ(ierr);  // allocate working space
  ierr = VecDuplicate(*H, &dbedElastic); CHKERRQ(ierr);  // allocate working space

  // allocate plate displacement
  ierr = VecCreateSeq(PETSC_COMM_SELF, Nx * Ny, &U); CHKERRQ(ierr);
  ierr = VecDuplicate(U, &U_start); CHKERRQ(ierr);
  // FFT - side coefficient fields (i.e. multiplication form of operators)
  ierr = VecDuplicate(U, &vleft); CHKERRQ(ierr);
  ierr = VecDuplicate(U, &vright); CHKERRQ(ierr);
  ierr = VecCreateSeq(PETSC_COMM_SELF, Nxge * Nyge, &lrmE); CHKERRQ(ierr);

  // setup fftw stuff: FFTW builds "plans" based on observed performance

  fftw_input  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx * Ny);
  fftw_output = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx * Ny);
  loadhat     = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx * Ny);

  // fftw manipulates the data in setting up a plan, so fill with nonconstant junk
  {
    VecAccessor2D<fftw_complex> tmp(fftw_input, Nx, Ny);
    for (PetscInt i = 0; i < Nx; i++) {
      for (PetscInt j = 0; j < Ny; j++) {
        tmp(i, j)[0] = i - 3;
        tmp(i, j)[1] = j*j + 2;
      }
    }
  }
  dft_forward = fftw_plan_dft_2d(Nx, Ny, fftw_input, fftw_output, FFTW_FORWARD, FFTW_MEASURE);
  dft_inverse = fftw_plan_dft_2d(Nx, Ny, fftw_input, fftw_output, FFTW_BACKWARD, FFTW_MEASURE);

  // coeffs for Fourier spectral method Laplacian
  // Matlab version:  cx=(pi/Lx)*[0:Nx/2 Nx/2-1:-1:1]
  cx = new PetscScalar[Nx];
  cy = new PetscScalar[Ny];

  for (PetscInt i = 0; i <= Nx / 2; i++)
    cx[i] = (pi / Lx_fat) * i;

  for (PetscInt i = Nx / 2 + 1; i < Nx; i++)
    cx[i] = (pi / Lx_fat) * (Nx - i);

  for (PetscInt j = 0; j <= Ny / 2; j++)
    cy[j] = (pi / Ly_fat) * j;

  for (PetscInt j = Ny / 2 + 1; j < Ny; j++)
    cy[j] = (pi / Ly_fat) * (Ny - j);

  // compare geforconv.m
  if (include_elastic == PETSC_TRUE) {
    ierr = PetscPrintf(PETSC_COMM_SELF,
           "     computing spherical elastic load response matrix ..."); CHKERRQ(ierr);
    PetscVecAccessor2D II(lrmE, Nxge, Nyge);
    ge_params ge_data;
    ge_data.dx = dx;
    ge_data.dy = dy;
    //const PetscInt imid_ge = Nx/2, jmid_ge = Ny/2;
    for (PetscInt i = 0; i < Nxge; i++) {
      for (PetscInt j = 0; j < Nyge; j++) {
        ge_data.p = i;
        ge_data.q = j;
        //ge_data.p = i - imid_ge;
        //ge_data.q = j - jmid_ge;
        II(i, j) = dblquad_cubature(ge_integrand, -dx/2, dx/2, -dy/2, dy/2,
                                    1.0e-8, &ge_data);
      }
    }

    ierr = PetscPrintf(PETSC_COMM_SELF, " done\n"); CHKERRQ(ierr);
  }

  allocDone = PETSC_TRUE;
  return 0;
}


PetscErrorCode BedDeformLC::uplift_init() {
  // to initialize we solve:
  //   rho_r g U + D grad^4 U = 0 - 2 eta |grad| uplift
  // where U=U_start; yes it really should have "0" on right side because
  // load for future times will be "-rho g (H-H_start)", which is zero if no geometry
  // change.

  // Compare equation (16) in Bueler, Lingle, Brown (2007) "Fast computation of
  // a viscoelastic deformable Earth model for ice sheet simulations", Ann.
  // Glaciol. 46, 97--105
  // [NOTE PROBABLE SIGN ERROR in eqn (16)?:  load "rho g H" should be "- rho g H"]

  // spectral/FFT quantities are on fat computational grid but uplift is on thin
  PetscErrorCode ierr;
  PetscVecAccessor2D left(vleft, Nx, Ny), right(vright, Nx, Ny);

  // fft2(uplift)
  clear_fftw_input();
  set_fftw_input(*uplift, 1.0, Mx, My, i0_plate, j0_plate);
  fftw_execute(dft_forward);

  // compute left and right coefficients
  for (PetscInt i = 0; i < Nx; i++) {
    for (PetscInt j = 0; j < Ny; j++) {
      const PetscScalar cclap = cx[i]*cx[i] + cy[j]*cy[j];
      left(i, j) = rho * standard_gravity + D * cclap * cclap;
      right(i, j) = -2.0 * eta * sqrt(cclap);
    }
  }

  // Matlab version:
  //        frhs = right.*fft2(uplift);
  //        u = real(ifft2( frhs. / left ));
  {
    VecAccessor2D<fftw_complex> u0_hat(fftw_input, Nx, Ny),
      uplift_hat(fftw_output, Nx, Ny);

    for (PetscInt i = 0; i < Nx; i++) {
      for (PetscInt j = 0; j < Ny; j++) {
        u0_hat(i, j)[0] = (right(i, j) * uplift_hat(i, j)[0]) / left(i, j);
        u0_hat(i, j)[1] = (right(i, j) * uplift_hat(i, j)[1]) / left(i, j);
      }
    }
  }

  fftw_execute(dft_inverse);
  get_fftw_output(U_start, 1.0 / (Nx * Ny), Nx, Ny, 0, 0);

  {
    PetscVecAccessor2D u_start(U_start, Nx, Ny);

    PetscScalar av = 0.0;
    for (PetscInt i = 0; i < Nx; i++)
      av += u_start(i, 0);

    for (PetscInt j = 0; j < Ny; j++)
      av += u_start(0, j);

    av = av / ((PetscScalar) (Nx + Ny));

    ierr = VecShift(U_start, -av); CHKERRQ(ierr);
  }

  ierr = VecCopy(U_start, U); CHKERRQ(ierr);

  return 0;
}


PetscErrorCode BedDeformLC::step(const PetscScalar dt_seconds, const PetscScalar seconds_from_start) {
  // solves:
  //     (2 eta |grad| U^{n+1}) + (dt/2) * ( rho_r g U^{n+1} + D grad^4 U^{n+1} )
  //   = (2 eta |grad| U^n) - (dt/2) * ( rho_r g U^n + D grad^4 U^n ) - dt * rho g H_start
  // where U=plate; see equation (7) in
  // Bueler, Lingle, Brown (2007) "Fast computation of a viscoelastic
  // deformable Earth model for ice sheet simulations", Ann. Glaciol. 46, 97--105

  // note ice thicknesses and bed elevations only on physical ("thin") grid
  //   while spectral/FFT quantities are on fat computational grid

  PetscVecAccessor2D left(vleft, Nx, Ny), right(vright, Nx, Ny);

  // Compute Hdiff
  PetscErrorCode ierr = VecWAXPY(Hdiff, -1, *H_start, *H); CHKERRQ(ierr);

  // Compute fft2(-ice_rho * g * dH * dt), where H = H - H_start.
  clear_fftw_input();
  set_fftw_input(Hdiff, - icerho * standard_gravity * dt_seconds,
                 Mx, My, i0_plate, j0_plate);
  fftw_execute(dft_forward);

  // Save fft2(-ice_rho * g * dH * dt) in loadhat.
  copy_fftw_output(loadhat);

  // Compute fft2(u).
  // no need to clear fftw_input: all values are overwritten
  set_fftw_input(U, 1.0, Nx, Ny, 0, 0);
  fftw_execute(dft_forward);

  // Compute left and right coefficients; note they depend on the length of a
  // time-step and thus cannot be precomputed
  for (PetscInt i = 0; i < Nx; i++) {
    for (PetscInt j = 0; j < Ny; j++) {
      const PetscScalar cclap = cx[i]*cx[i] + cy[j]*cy[j],
        part1 = 2.0 * eta * sqrt(cclap),
        part2 = (dt_seconds / 2.0) * (rho * standard_gravity + D * cclap * cclap);
      left(i, j)  = part1 + part2;
      right(i, j) = part1 - part2;
    }
  }

  //         frhs = right.*fft2(uun) + fft2(dt*sszz);
  //         uun1 = real(ifft2( frhs./left ));
  {
    VecAccessor2D<fftw_complex> input(fftw_input, Nx, Ny),
      u_hat(fftw_output, Nx, Ny), load_hat(loadhat, Nx, Ny);
    for (PetscInt i = 0; i < Nx; i++) {
      for (PetscInt j = 0; j < Ny; j++) {
        input(i, j)[0] = (right(i, j) * u_hat(i, j)[0] + load_hat(i, j)[0]) / left(i, j);
        input(i, j)[1] = (right(i, j) * u_hat(i, j)[1] + load_hat(i, j)[1]) / left(i, j);
      }
    }
  }

  fftw_execute(dft_inverse);
  get_fftw_output(U, 1.0 / (Nx * Ny), Nx, Ny, 0, 0);

  // now tweak
  tweak(seconds_from_start);

  // now compute elastic response if desired; bed = ue at end of this block
  if (include_elastic == PETSC_TRUE) {
    // Matlab:     ue=rhoi*conv2(H-H_start, II, 'same')
    ierr = conv2_same(Hdiff, Mx, My, lrmE, Nxge, Nyge, dbedElastic);  CHKERRQ(ierr);
    ierr = VecScale(dbedElastic, icerho);  CHKERRQ(ierr);
  } else {
    ierr = VecSet(dbedElastic, 0.0); CHKERRQ(ierr);
  }

  // now sum contributions to get new bed elevation:
  //    (new bed) = ue + (bed start) + plate
  // (but use only central part of plate if Z>1)
  {
    PetscVecAccessor2D b(*bed, Mx, My), b_start(*bed_start, Mx, My), db_elastic(dbedElastic, Mx, My),
      u(U, Nx, Ny, i0_plate, j0_plate), u_start(U_start, Nx, Ny, i0_plate, j0_plate);

    for (PetscInt i = 0; i < Mx; i++) {
      for (PetscInt j = 0; j < My; j++) {
        b(i, j) = b_start(i, j) + db_elastic(i, j) + (u(i, j) - u_start(i, j));
      }
    }
  }

  return 0;
}

void BedDeformLC::tweak(PetscReal seconds_from_start) {
  PetscVecAccessor2D u(U, Nx, Ny);

  // find average value along "distant" boundary of [-Lx_fat, Lx_fat]X[-Ly_fat, Ly_fat]
  // note domain is periodic, so think of cut locus of torus (!)
  // (will remove it:   uun1=uun1-( sum(uun1(1, :))+sum(uun1(:, 1)) )/(2*N);)
  PetscScalar av = 0.0;
  for (PetscInt i = 0; i < Nx; i++)
    av += u(i, 0);

  for (PetscInt j = 0; j < Ny; j++)
    av += u(0, j);

  av = av / ((PetscScalar) (Nx + Ny));

  // tweak continued: replace far field with value for an equivalent disc load which has R0=Lx*(2/3)=L/3
  // (instead of 1000km in Matlab code: H0 = dx*dx*sum(sum(H))/(pi*1e6^2);  % trapezoid rule)
  const PetscScalar Lav = (Lx_fat + Ly_fat) / 2.0;
  const PetscScalar Requiv = Lav * (2.0 / 3.0);
  PetscScalar delvolume;
  VecSum(Hdiff, &delvolume);
  delvolume = delvolume * dx * dy;  // make into a volume
  const PetscScalar Hequiv = delvolume / (pi * Requiv * Requiv);

  const PetscScalar discshift = viscDisc(seconds_from_start,
                                         Hequiv, Requiv, Lav, rho, standard_gravity, D, eta) - av;

  VecShift(U, discshift);

}

//! \brief Fill fftw_input with zeros.
void BedDeformLC::clear_fftw_input() {
  VecAccessor2D<fftw_complex> fftw_in(fftw_input, Nx, Ny);
  for (int i = 0; i < Nx; ++i) {
    for (int j = 0; j < Ny; ++j) {
      fftw_in(i, j)[0] = 0;
      fftw_in(i, j)[1] = 0;
    }
  }
}

//! \brief Copy fftw_output to \c output.
void BedDeformLC::copy_fftw_output(fftw_complex *output) {
  VecAccessor2D<fftw_complex> fftw_out(fftw_output, Nx, Ny), out(output, Nx, Ny);
  for (int i = 0; i < Nx; ++i) {
    for (int j = 0; j < Ny; ++j) {
      out(i, j)[0] = fftw_out(i, j)[0];
      out(i, j)[1] = fftw_out(i, j)[1];
    }
  }
}

//! \brief Set the real part of fftw_input to vec_input.
/*!
 * Sets the imaginary part to zero.
 */
void BedDeformLC::set_fftw_input(Vec vec_input, PetscReal normalization, int M, int N, int i0, int j0) {
  PetscVecAccessor2D in(vec_input, M, N);
  VecAccessor2D<fftw_complex> input(fftw_input, Nx, Ny, i0, j0);
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < M; ++j) {
      input(i, j)[0] = in(i, j) * normalization;
      input(i, j)[1] = 0.0;
    }
  }
}

//! \brief Get the real part of fftw_output and put it in output.
void BedDeformLC::get_fftw_output(Vec output, PetscReal normalization, int M, int N, int i0, int j0) {
  PetscVecAccessor2D out(output, M, N);
  VecAccessor2D<fftw_complex> fftw_out(fftw_output, Nx, Ny, i0, j0);
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < M; ++j) {
      out(i, j) = fftw_out(i, j)[0] * normalization;
    }
  }
}
