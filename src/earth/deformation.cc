// Copyright (C) 2004-2009, 2011, 2013, 2014, 2015, 2016 Ed Bueler and Constantine Khroulev
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
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
#include <fftw3.h>
#include <cassert>
#include <gsl/gsl_math.h>       // M_PI

#include "base/util/pism_const.hh"
#include "matlablike.hh"
#include "greens.hh"
#include "deformation.hh"
#include "base/util/PISMConfigInterface.hh"
#include "base/util/error_handling.hh"
#include "base/util/petscwrappers/Vec.hh"

namespace pism {
namespace bed {

template <class T>
class VecAccessor2D {
public:
  VecAccessor2D(T* a, int my_Mx, int my_My, int my_i_offset, int my_j_offset)
    : m_Mx(my_Mx), m_My(my_My), m_i_offset(my_i_offset), m_j_offset(my_j_offset), m_array(a) {}

  VecAccessor2D(T* a, int my_Mx, int my_My)
    : m_Mx(my_Mx), m_My(my_My), m_i_offset(0), m_j_offset(0), m_array(a) {}

  inline T& operator()(int i, int j) {
    return m_array[(m_j_offset + j) + m_My * (m_i_offset + i)];
  }

private:
  int m_Mx, m_My, m_i_offset, m_j_offset;
  T* m_array;
};

BedDeformLC::BedDeformLC(const Config &config,
                         bool myinclude_elastic,
                         int myMx, int myMy,
                         double mydx, double mydy,
                         int myZ,
                         Vec myHstart, Vec mybedstart, Vec myuplift,
                         Vec myH, Vec mybed) {

  // set parameters
  m_include_elastic = myinclude_elastic;

  m_Mx     = myMx;
  m_My     = myMy;
  m_dx     = mydx;
  m_dy     = mydy;
  m_Z      = myZ;
  m_icerho = config.get_double("ice_density");
  m_rho    = config.get_double("lithosphere_density");
  m_eta    = config.get_double("mantle_viscosity");
  m_D      = config.get_double("lithosphere_flexural_rigidity");

  m_standard_gravity = config.get_double("standard_gravity");

  // derive more parameters
  m_Lx       = ((m_Mx - 1) / 2) * m_dx;
  m_Ly       = ((m_My - 1) / 2) * m_dy;
  m_Nx       = m_Z*(m_Mx - 1);
  m_Ny       = m_Z*(m_My - 1);
  m_Lx_fat   = (m_Nx / 2) *   m_dx;
  m_Ly_fat   = (m_Ny / 2) *   m_dy;
  m_Nxge     = m_Nx + 1;
  m_Nyge     = m_Ny + 1;
  m_i0_plate = (m_Z - 1)*(m_Mx - 1) / 2;
  m_j0_plate = (m_Z - 1)*(m_My - 1) / 2;

  // attach to existing (must be allocated!) sequential Vecs
  m_H         = myH;
  m_bed       = mybed;
  m_H_start   = myHstart;
  m_bed_start = mybedstart;
  m_uplift    = myuplift;

  // memory allocation
  PetscErrorCode  ierr;

  ierr = VecDuplicate(m_H, m_Hdiff.rawptr());
  PISM_CHK(ierr, "VecDuplicate");

  ierr = VecDuplicate(m_H, m_dbedElastic.rawptr());
  PISM_CHK(ierr, "VecDuplicate");

  // allocate plate displacement
  ierr = VecCreateSeq(PETSC_COMM_SELF, m_Nx * m_Ny, m_U.rawptr());
  PISM_CHK(ierr, "VecCreateSeq");

  ierr = VecDuplicate(m_U, m_U_start.rawptr());
  PISM_CHK(ierr, "VecDuplicate");

  // FFT - side coefficient fields (i.e. multiplication form of operators)
  ierr = VecDuplicate(m_U, m_vleft.rawptr());
  PISM_CHK(ierr, "VecDuplicate");

  ierr = VecDuplicate(m_U, m_vright.rawptr());
  PISM_CHK(ierr, "VecDuplicate");

  ierr = VecCreateSeq(PETSC_COMM_SELF, m_Nxge * m_Nyge, m_lrmE.rawptr());
  PISM_CHK(ierr, "VecCreateSeq");

  // setup fftw stuff: FFTW builds "plans" based on observed performance

  m_fftw_input  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * m_Nx * m_Ny);
  m_fftw_output = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * m_Nx * m_Ny);
  m_loadhat     = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * m_Nx * m_Ny);

  // fftw manipulates the data in setting up a plan, so fill with nonconstant junk
  {
    VecAccessor2D<fftw_complex> tmp(m_fftw_input, m_Nx, m_Ny);
    for (int j = 0; j < m_Ny; j++) {
      for (int i = 0; i < m_Nx; i++) {
        tmp(i, j)[0] = i - 3;
        tmp(i, j)[1] = j*j + 2;
      }
    }
  }

  // Limit the amount of time FFTW is allowed to spend choosing algorithms.
  fftw_set_timelimit(60.0);

  m_dft_forward = fftw_plan_dft_2d(m_Nx, m_Ny, m_fftw_input, m_fftw_output,
                                   FFTW_FORWARD, FFTW_MEASURE);
  m_dft_inverse = fftw_plan_dft_2d(m_Nx, m_Ny, m_fftw_input, m_fftw_output,
                                   FFTW_BACKWARD, FFTW_MEASURE);

  // Note: FFTW is weird. If a malloc() call fails it will just call
  // abort() on you without giving you a chance to recover or tell the
  // user what happened. This is why we don't check return values of
  // fftw_malloc() and fftw_plan_dft_2d() calls here...
  //
  // (Constantine Khroulev, February 1, 2015)

  m_cx.resize(m_Nx);
  m_cy.resize(m_Ny);

  precompute_coefficients();
}

BedDeformLC::~BedDeformLC() {
  fftw_destroy_plan(m_dft_forward);
  fftw_destroy_plan(m_dft_inverse);
  fftw_free(m_fftw_input);
  fftw_free(m_fftw_output);
  fftw_free(m_loadhat);
}

/**
 * Pre-compute coefficients used by the model.
 */
void BedDeformLC::precompute_coefficients() {
  PetscErrorCode ierr;

  // coeffs for Fourier spectral method Laplacian
  // Matlab version:  cx=(pi/Lx)*[0:Nx/2 Nx/2-1:-1:1]
  for (int i = 0; i <= m_Nx / 2; i++) {
    m_cx[i] = (M_PI / m_Lx_fat) * i;
  }

  for (int i = m_Nx / 2 + 1; i < m_Nx; i++) {
    m_cx[i] = (M_PI / m_Lx_fat) * (m_Nx - i);
  }

  for (int j = 0; j <= m_Ny / 2; j++) {
    m_cy[j] = (M_PI / m_Ly_fat) * j;
  }

  for (int j = m_Ny / 2 + 1; j < m_Ny; j++) {
    m_cy[j] = (M_PI / m_Ly_fat) * (m_Ny - j);
  }

  // compare geforconv.m
  if (m_include_elastic == true) {
    ierr = PetscPrintf(PETSC_COMM_SELF,
                       "     computing spherical elastic load response matrix ...");
    PISM_CHK(ierr, "PetscPrintf");

    petsc::VecArray2D II(m_lrmE, m_Nxge, m_Nyge);
    ge_params ge_data;
    ge_data.dx = m_dx;
    ge_data.dy = m_dy;
    for (int j = 0; j < m_Nyge; j++) {
      for (int i = 0; i < m_Nxge; i++) {
        ge_data.p = i;
        ge_data.q = j;
        II(i, j) = dblquad_cubature(ge_integrand, -m_dx/2, m_dx/2, -m_dy/2, m_dy/2,
                                    1.0e-8, &ge_data);
      }
    }

    ierr = PetscPrintf(PETSC_COMM_SELF, " done\n");
    PISM_CHK(ierr, "PetscPrintf");
  }
}

void BedDeformLC::uplift_init() {
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
  petsc::VecArray2D left(m_vleft, m_Nx, m_Ny), right(m_vright, m_Nx, m_Ny);

  // fft2(uplift)
  clear_fftw_input();
  set_fftw_input(m_uplift, 1.0, m_Mx, m_My, m_i0_plate, m_j0_plate);
  fftw_execute(m_dft_forward);

  // compute left and right coefficients
  for (int j = 0; j < m_Ny; j++) {
    for (int i = 0; i < m_Nx; i++) {
      const double cclap = m_cx[i]*m_cx[i] + m_cy[j]*m_cy[j];
      left(i, j) = m_rho * m_standard_gravity + m_D * cclap * cclap;
      right(i, j) = -2.0 * m_eta * sqrt(cclap);
    }
  }

  // Matlab version:
  //        frhs = right.*fft2(uplift);
  //        u = real(ifft2(frhs. / left));
  {
    VecAccessor2D<fftw_complex> u0_hat(m_fftw_input, m_Nx, m_Ny),
      uplift_hat(m_fftw_output, m_Nx, m_Ny);

    for (int j = 0; j < m_Ny; j++) {
      for (int i = 0; i < m_Nx; i++) {
        u0_hat(i, j)[0] = (right(i, j) * uplift_hat(i, j)[0]) / left(i, j);
        u0_hat(i, j)[1] = (right(i, j) * uplift_hat(i, j)[1]) / left(i, j);
      }
    }
  }

  fftw_execute(m_dft_inverse);
  get_fftw_output(m_U_start, 1.0 / (m_Nx * m_Ny), m_Nx, m_Ny, 0, 0);

  {
    petsc::VecArray2D u_start(m_U_start, m_Nx, m_Ny);

    double av = 0.0;
    for (int i = 0; i < m_Nx; i++) {
      av += u_start(i, 0);
    }

    for (int j = 0; j < m_Ny; j++) {
      av += u_start(0, j);
    }

    av = av / ((double) (m_Nx + m_Ny));

    ierr = VecShift(m_U_start, -av);
    PISM_CHK(ierr, "VecShift");
  }

  ierr = VecCopy(m_U_start, m_U);
  PISM_CHK(ierr, "VecCopy");
}


void BedDeformLC::step(double dt_seconds, double seconds_from_start) {
  // solves:
  //     (2 eta |grad| U^{n+1}) + (dt/2) * (rho_r g U^{n+1} + D grad^4 U^{n+1})
  //   = (2 eta |grad| U^n) - (dt/2) * (rho_r g U^n + D grad^4 U^n) - dt * rho g H_start
  // where U=plate; see equation (7) in
  // Bueler, Lingle, Brown (2007) "Fast computation of a viscoelastic
  // deformable Earth model for ice sheet simulations", Ann. Glaciol. 46, 97--105

  // note ice thicknesses and bed elevations only on physical ("thin") grid
  //   while spectral/FFT quantities are on fat computational grid

  petsc::VecArray2D left(m_vleft, m_Nx, m_Ny), right(m_vright, m_Nx, m_Ny);

  // Compute Hdiff
  PetscErrorCode ierr = VecWAXPY(m_Hdiff, -1, m_H_start, m_H);
  PISM_CHK(ierr, "VecWAXPY");

  // Compute fft2(-ice_rho * g * dH * dt), where H = H - H_start.
  clear_fftw_input();
  set_fftw_input(m_Hdiff, - m_icerho * m_standard_gravity * dt_seconds,
                 m_Mx, m_My, m_i0_plate, m_j0_plate);
  fftw_execute(m_dft_forward);

  // Save fft2(-ice_rho * g * dH * dt) in loadhat.
  copy_fftw_output(m_loadhat);

  // Compute fft2(u).
  // no need to clear fftw_input: all values are overwritten
  set_fftw_input(m_U, 1.0, m_Nx, m_Ny, 0, 0);
  fftw_execute(m_dft_forward);

  // Compute left and right coefficients; note they depend on the length of a
  // time-step and thus cannot be precomputed
  for (int j = 0; j < m_Ny; j++) {
    for (int i = 0; i < m_Nx; i++) {
      const double cclap = m_cx[i]*m_cx[i] + m_cy[j]*m_cy[j],
        part1 = 2.0 * m_eta * sqrt(cclap),
        part2 = (dt_seconds / 2.0) * (m_rho * m_standard_gravity + m_D * cclap * cclap);
      left(i, j)  = part1 + part2;
      right(i, j) = part1 - part2;
    }
  }

  //         frhs = right.*fft2(uun) + fft2(dt*sszz);
  //         uun1 = real(ifft2(frhs./left));
  {
    VecAccessor2D<fftw_complex> input(m_fftw_input, m_Nx, m_Ny),
      u_hat(m_fftw_output, m_Nx, m_Ny), load_hat(m_loadhat, m_Nx, m_Ny);
    for (int j = 0; j < m_Ny; j++) {
      for (int i = 0; i < m_Nx; i++) {
        input(i, j)[0] = (right(i, j) * u_hat(i, j)[0] + load_hat(i, j)[0]) / left(i, j);
        input(i, j)[1] = (right(i, j) * u_hat(i, j)[1] + load_hat(i, j)[1]) / left(i, j);
      }
    }
  }

  fftw_execute(m_dft_inverse);
  get_fftw_output(m_U, 1.0 / (m_Nx * m_Ny), m_Nx, m_Ny, 0, 0);

  // now tweak
  tweak(seconds_from_start);

  // now compute elastic response if desired; bed = ue at end of this block
  if (m_include_elastic == true) {
    // Matlab:     ue=rhoi*conv2(H-H_start, II, 'same')
    conv2_same(m_Hdiff, m_Mx, m_My, m_lrmE, m_Nxge, m_Nyge, m_dbedElastic);

    ierr = VecScale(m_dbedElastic, m_icerho);
    PISM_CHK(ierr, "VecScale");
  } else {
    ierr = VecSet(m_dbedElastic, 0.0);
    PISM_CHK(ierr, "VecSet");
  }

  // now sum contributions to get new bed elevation:
  //    (new bed) = ue + (bed start) + plate
  // (but use only central part of plate if Z>1)
  {
    petsc::VecArray2D b(m_bed, m_Mx, m_My), b_start(m_bed_start, m_Mx, m_My), db_elastic(m_dbedElastic, m_Mx, m_My),
      u(m_U, m_Nx, m_Ny, m_i0_plate, m_j0_plate), u_start(m_U_start, m_Nx, m_Ny, m_i0_plate, m_j0_plate);

    for (int j = 0; j < m_My; j++) {
      for (int i = 0; i < m_Mx; i++) {
        b(i, j) = b_start(i, j) + db_elastic(i, j) + (u(i, j) - u_start(i, j));
      }
    }
  }
}

void BedDeformLC::tweak(double seconds_from_start) {
  petsc::VecArray2D u(m_U, m_Nx, m_Ny);

  // find average value along "distant" boundary of [-Lx_fat, Lx_fat]X[-Ly_fat, Ly_fat]
  // note domain is periodic, so think of cut locus of torus (!)
  // (will remove it:   uun1=uun1-(sum(uun1(1, :))+sum(uun1(:, 1)))/(2*N);)
  double av = 0.0;
  for (int i = 0; i < m_Nx; i++) {
    av += u(i, 0);
  }

  for (int j = 0; j < m_Ny; j++) {
    av += u(0, j);
  }

  av = av / ((double) (m_Nx + m_Ny));

  // tweak continued: replace far field with value for an equivalent disc load which has R0=Lx*(2/3)=L/3
  // (instead of 1000km in Matlab code: H0 = dx*dx*sum(sum(H))/(pi*1e6^2);  % trapezoid rule)
  const double Lav = (m_Lx_fat + m_Ly_fat) / 2.0;
  const double Requiv = Lav * (2.0 / 3.0);

  double delvolume;
  PetscErrorCode ierr = VecSum(m_Hdiff, &delvolume);
  PISM_CHK(ierr, "VecSum");

  delvolume = delvolume * m_dx * m_dy;  // make into a volume
  const double Hequiv = delvolume / (M_PI * Requiv * Requiv);

  const double discshift = viscDisc(seconds_from_start,
                                    Hequiv, Requiv, Lav, m_rho, m_standard_gravity, m_D, m_eta) - av;

  ierr = VecShift(m_U, discshift);
  PISM_CHK(ierr, "VecShift");
}

//! \brief Fill fftw_input with zeros.
void BedDeformLC::clear_fftw_input() {
  VecAccessor2D<fftw_complex> fftw_in(m_fftw_input, m_Nx, m_Ny);
  for (int j = 0; j < m_Ny; ++j) {
    for (int i = 0; i < m_Nx; ++i) {
      fftw_in(i, j)[0] = 0;
      fftw_in(i, j)[1] = 0;
    }
  }
}

//! \brief Copy fftw_output to `output`.
void BedDeformLC::copy_fftw_output(fftw_complex *output) {
  VecAccessor2D<fftw_complex> fftw_out(m_fftw_output, m_Nx, m_Ny), out(output, m_Nx, m_Ny);
  for (int j = 0; j < m_Ny; ++j) {
    for (int i = 0; i < m_Nx; ++i) {
      out(i, j)[0] = fftw_out(i, j)[0];
      out(i, j)[1] = fftw_out(i, j)[1];
    }
  }
}

//! \brief Set the real part of fftw_input to vec_input.
/*!
 * Sets the imaginary part to zero.
 */
void BedDeformLC::set_fftw_input(Vec vec_input, double normalization, int M, int N, int i0, int j0) {
  petsc::VecArray2D in(vec_input, M, N);
  VecAccessor2D<fftw_complex> input(m_fftw_input, m_Nx, m_Ny, i0, j0);
  for (int j = 0; j < N; ++j) {
    for (int i = 0; i < M; ++i) {
      input(i, j)[0] = in(i, j) * normalization;
      input(i, j)[1] = 0.0;
    }
  }
}

//! \brief Get the real part of fftw_output and put it in output.
void BedDeformLC::get_fftw_output(Vec output, double normalization, int M, int N, int i0, int j0) {
  petsc::VecArray2D out(output, M, N);
  VecAccessor2D<fftw_complex> fftw_out(m_fftw_output, m_Nx, m_Ny, i0, j0);
  for (int j = 0; j < N; ++j) {
    for (int i = 0; i < M; ++i) {
      out(i, j) = fftw_out(i, j)[0] * normalization;
    }
  }
}

} // end of namespace bed
} // end of namespace pism
