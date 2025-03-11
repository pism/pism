// Copyright (C) 2004-2009, 2011, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2023, 2025 Ed Bueler and Constantine Khroulev
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

#include <cassert>
#include <cmath>                // sqrt
#include <fftw3.h>
#include <gsl/gsl_math.h>       // M_PI

#include "pism/earth/matlablike.hh"
#include "pism/earth/greens.hh"
#include "pism/earth/LingleClarkSerial.hh"

#include "pism/util/Config.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/petscwrappers/Vec.hh"
#include "pism/util/fftw_utilities.hh"

namespace pism {
namespace bed {

/*!
 * @param[in] config configuration database
 * @param[in] include_elastic include elastic deformation component
 * @param[in] Mx grid size in the X direction
 * @param[in] My grid size in the Y direction
 * @param[in] dx grid spacing in the X direction
 * @param[in] dy grid spacing in the Y direction
 * @param[in] Nx extended grid size in the X direction
 * @param[in] Ny extended grid size in the Y direction
 */
LingleClarkSerial::LingleClarkSerial(Logger::ConstPtr log,
                                     const Config &config,
                                     bool include_elastic,
                                     int Mx, int My,
                                     double dx, double dy,
                                     int Nx, int Ny)
  : m_t_infty(1e16),            // around 317 million years
    m_log(log) {

  // set parameters
  m_include_elastic = include_elastic;

  if (include_elastic) {
    // check if the extended grid is large enough (it has to be at least twice the size of
    // the physical grid so that the load in one corner of the domain affects the grid
    // point in the opposite corner).

    if (config.get_number("bed_deformation.lc.grid_size_factor") < 2) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                    "bed_deformation.lc.elastic_model"
                                    " requires bed_deformation.lc.grid_size_factor > 1");
    }
  }

  // grid parameters
  m_Mx = Mx;
  m_My = My;
  m_dx = dx;
  m_dy = dy;
  m_Nx = Nx;
  m_Ny = Ny;

  m_load_density   = config.get_number("constants.ice.density");
  m_mantle_density = config.get_number("bed_deformation.mantle_density");
  m_eta            = config.get_number("bed_deformation.mantle_viscosity");
  m_D              = config.get_number("bed_deformation.lithosphere_flexural_rigidity");

  m_standard_gravity = config.get_number("constants.standard_gravity");

  // derive more parameters
  m_Lx        = 0.5 * (m_Nx - 1.0) * m_dx;
  m_Ly        = 0.5 * (m_Ny - 1.0) * m_dy;
  m_i0_offset = (Nx - Mx) / 2;
  m_j0_offset = (Ny - My) / 2;

  // memory allocation
  PetscErrorCode ierr = 0;

  // total displacement
  ierr = VecCreateSeq(PETSC_COMM_SELF, m_Mx * m_My, m_U.rawptr());;
  PISM_CHK(ierr, "VecCreateSeq");

  // elastic displacement
  ierr = VecCreateSeq(PETSC_COMM_SELF, m_Mx * m_My, m_Ue.rawptr());;
  PISM_CHK(ierr, "VecCreateSeq");

  // viscous displacement
  ierr = VecCreateSeq(PETSC_COMM_SELF, m_Nx * m_Ny, m_Uv.rawptr());
  PISM_CHK(ierr, "VecCreateSeq");

  // setup fftw stuff: FFTW builds "plans" based on observed performance
  m_fftw_input  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * m_Nx * m_Ny);
  m_fftw_output = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * m_Nx * m_Ny);
  m_loadhat     = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * m_Nx * m_Ny);
  m_lrm_hat = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * m_Nx * m_Ny);

  clear_fftw_array(m_fftw_input, m_Nx, m_Ny);
  m_dft_forward = fftw_plan_dft_2d(m_Nx, m_Ny, m_fftw_input, m_fftw_output,
                                   FFTW_FORWARD, FFTW_ESTIMATE);
  m_dft_inverse = fftw_plan_dft_2d(m_Nx, m_Ny, m_fftw_input, m_fftw_output,
                                   FFTW_BACKWARD, FFTW_ESTIMATE);

  // Note: FFTW is weird. If a malloc() call fails it will just call
  // abort() on you without giving you a chance to recover or tell the
  // user what happened. This is why we don't check return values of
  // fftw_malloc() and fftw_plan_dft_2d() calls here...
  //
  // (Constantine Khroulev, February 1, 2015)

  precompute_coefficients();
}

LingleClarkSerial::~LingleClarkSerial() {
  fftw_destroy_plan(m_dft_forward);
  fftw_destroy_plan(m_dft_inverse);
  fftw_free(m_fftw_input);
  fftw_free(m_fftw_output);
  fftw_free(m_loadhat);
  fftw_free(m_lrm_hat);
}

/*!
 * Return total displacement.
 */
const petsc::Vec& LingleClarkSerial::total_displacement() const {
  return m_U;
}

/*!
 * Return viscous plate displacement.
 */
const petsc::Vec& LingleClarkSerial::viscous_displacement() const {
  return m_Uv;
}

/*!
 * Return elastic plate displacement.
 */
const petsc::Vec& LingleClarkSerial::elastic_displacement() const {
  return m_Ue;
}

void LingleClarkSerial::compute_load_response_matrix(fftw_complex *output) {

  FFTWArray LRM(output, m_Nx, m_Ny);

  greens_elastic G;
  ge_data ge_data {m_dx, m_dy, 0, 0, &G};

  int Nx2 = m_Nx / 2;
  int Ny2 = m_Ny / 2;

  // Top half
  for (int j = 0; j <= Ny2; ++j) {
    // Top left quarter
    for (int i = 0; i <= Nx2; ++i) {
      ge_data.p = Nx2 - i;
      ge_data.q = Ny2 - j;

      LRM(i, j) = dblquad_cubature(ge_integrand,
                                   -m_dx / 2, m_dx / 2,
                                   -m_dy / 2, m_dy / 2,
                                   1.0e-8, &ge_data);
    }

    // Top right quarter
    //
    // Note: Nx2 = m_Nx / 2 (using integer division!), so
    //
    // - If m_Nx is even then 2 * Nx2 == m_Nx. So i < m_Nx implies i < 2 * Nx2 and
    //   2 * Nx2 - i > 0.
    //
    // - If m_Nx is odd then 2 * Nx2 == m_Nx - 1 or m_Nx == 2 * Nx2 + 1. So i < m_Nx
    //   implies i < 2 * Nx2 + 1, which is the same as 2 * Nx2 - i > -1 or
    //   2 * Nx2 - i >= 0.
    //
    // Also, i == Nx2 + 1 gives 2 * Nx2 - i == Nx2 - 1
    //
    // So, in both cases (even and odd) 0 <= 2 * Nx2 - i <= Nx2 - 1.
    //
    // This means that LRM(2 * Nx2 - i, j) will not use indexes that are out of bounds
    // *and* will only use values computed in the for loop above.
    for (int i = Nx2 + 1; i < m_Nx; ++i) {
      assert(2 * Nx2 - i >= 0);
      LRM(i, j) = LRM(2 * Nx2 - i, j);
    }
  } // End of the loop over the top half

    // Bottom half
    //
    // See the comment above the "top right quarter" loop.
  for (int j = Ny2 + 1; j < m_Ny; ++j) {
    for (int i = 0; i < m_Nx; ++i) {
      assert(2 * Ny2 - j >= 0);
      LRM(i, j) = LRM(i, 2 * Ny2 - j);
    }
  }
}

/**
 * Pre-compute coefficients used by the model.
 */
void LingleClarkSerial::precompute_coefficients() {

  // Coefficients for Fourier spectral method Laplacian
  // MATLAB version:  cx=(pi/Lx)*[0:Nx/2 Nx/2-1:-1:1]
  m_cx = fftfreq(m_Nx, m_Lx / (m_Nx * M_PI));
  m_cy = fftfreq(m_Ny, m_Ly / (m_Ny * M_PI));

  // compare geforconv.m
  if (m_include_elastic) {
    m_log->message(2, "     computing spherical elastic load response matrix ...");
    {
      compute_load_response_matrix(m_fftw_input);
      // Compute fft2(LRM) and save it in m_lrm_hat
      fftw_execute(m_dft_forward);
      copy_fftw_array(m_fftw_output, m_lrm_hat, m_Nx, m_Ny);
    }
    m_log->message(2, " done\n");
  }
}

/*!
 * Solve
 *
 * @f$ 2 \nu |\nabla| \diff{u}{t} + \rho_r g U + D\nabla^4 U = \sigma_{zz}@f$
 *
 * for @f$ U @f$, treating @f$ \diff{u}{t} @f$ and @f$ \sigma_{zz} @f$ as known.
 *
 * @param[in] load_thickness load thickness, meters
 * @param[in] bed_uplift bed uplift, m/second
 *
 * Here `load_thickness` is used to compute the load @f$ \sigma_{zz} @f$ and `bed_uplift` is
 * @f$ \diff{u}{t} @f$ itself.
 *
 */
void LingleClarkSerial::uplift_problem(petsc::Vec &load_thickness,
                                       petsc::Vec &bed_uplift,
                                       petsc::Vec &output) {

  // Compute fft2(-load_density * g * load_thickness)
  {
    clear_fftw_array(m_fftw_input, m_Nx, m_Ny);
    set_real_part(load_thickness, - m_load_density * m_standard_gravity,
                  m_Mx, m_My, m_Nx, m_Ny, m_i0_offset, m_j0_offset,
                  m_fftw_input);
    fftw_execute(m_dft_forward);
    // Save fft2(-load_density * g * load_thickness) in loadhat.
    copy_fftw_array(m_fftw_output, m_loadhat, m_Nx, m_Ny);
  }

  // fft2(uplift)
  {
    clear_fftw_array(m_fftw_input, m_Nx, m_Ny);
    set_real_part(bed_uplift, 1.0, m_Mx, m_My, m_Nx, m_Ny, m_i0_offset, m_j0_offset,
                  m_fftw_input);
    fftw_execute(m_dft_forward);
  }

  {
    FFTWArray
      u0_hat(m_fftw_input, m_Nx, m_Ny),
      load_hat(m_loadhat, m_Nx, m_Ny),
      uplift_hat(m_fftw_output, m_Nx, m_Ny);

    for (int i = 0; i < m_Nx; i++) {
      for (int j = 0; j < m_Ny; j++) {
        const double
          C = m_cx[i]*m_cx[i] + m_cy[j]*m_cy[j],
          A = - 2.0 * m_eta * sqrt(C),
          B = m_mantle_density * m_standard_gravity + m_D * C * C;

        u0_hat(i, j) = (load_hat(i, j) + A * uplift_hat(i, j)) / B;
      }
    }
  }

  fftw_execute(m_dft_inverse);
  get_real_part(m_fftw_output, 1.0 / (m_Nx * m_Ny), m_Nx, m_Ny, m_Nx, m_Ny, 0, 0, output);

  tweak(load_thickness, output, m_Nx, m_Ny);
}

/*! Initialize using provided load thickness and the bed uplift rate.
 *
 * Here we solve:
 *
 *   rho_r g U + D grad^4 U = -rho g H - 2 eta |grad| uplift
 *
 * Compare equation (16) in Bueler, Lingle, Brown (2007) "Fast computation of a viscoelastic
 * deformable Earth model for ice sheet simulations", Ann. Glaciol. 46, 97--105
 *
 * @note Probable sign error in eqn (16)?:  load "rho g H" should be "- rho g H"]
 *
 * This initialization method is used to "bootstrap" the model. It should not be used to re-start a
 * stopped modeling run.
 *
 * @param[in] thickness load thickness, meters
 * @param[in] uplift initial bed uplift on the PISM grid
 *
 * Sets m_Uv, m_Ue, m_U.
 */
void LingleClarkSerial::bootstrap(petsc::Vec &thickness, petsc::Vec &uplift) {

  // compute viscous displacement
  uplift_problem(thickness, uplift, m_Uv);

  if (m_include_elastic) {
    compute_elastic_response(thickness, m_Ue);
  } else {
    PetscErrorCode ierr = VecSet(m_Ue, 0.0); PISM_CHK(ierr, "VecSet");
  }

  update_displacement(m_Uv, m_Ue, m_U);
}

/*!
 * Initialize using provided plate displacement.
 *
 * @param[in] viscous_displacement initial viscous plate displacement (meters) on the extended grid
 * @param[in] elastic_displacement initial viscous plate displacement (meters) on the regular grid
 *
 * Sets m_Uv, m_Ue, m_U.
 */
void LingleClarkSerial::init(petsc::Vec &viscous_displacement,
                             petsc::Vec &elastic_displacement) {
  PetscErrorCode ierr = 0;

  ierr = VecCopy(viscous_displacement, m_Uv);  PISM_CHK(ierr, "VecCopy");

  if (m_include_elastic) {
    ierr = VecCopy(elastic_displacement, m_Ue);  PISM_CHK(ierr, "VecCopy");
  } else {
    ierr = VecSet(m_Ue, 0.0); PISM_CHK(ierr, "VecSet");
  }

  update_displacement(m_Uv, m_Ue, m_U);
}

/*!
 * Perform a time step.
 *
 * @param[in] dt time step length
 * @param[in] H load thickness on the physical (Mx*My) grid
 */
void LingleClarkSerial::step(double dt, petsc::Vec &H) {
  // solves:
  //     (2 eta |grad| U^{n+1}) + (dt/2) * (rho_r g U^{n+1} + D grad^4 U^{n+1})
  //   = (2 eta |grad| U^n) - (dt/2) * (rho_r g U^n + D grad^4 U^n) - dt * rho g H_start
  // where U=plate displacement; see equation (7) in
  // Bueler, Lingle, Brown (2007) "Fast computation of a viscoelastic
  // deformable Earth model for ice sheet simulations", Ann. Glaciol. 46, 97--105

  // Compute viscous displacement if dt > 0 and bypass this computation if dt == 0.
  //
  // This makes it easier to test the elastic part of the model.
  if (dt > 0.0) {
    // Non-zero time step: include the viscous part of the model.

    // Compute fft2(-load_density * g * dt * H)
    {
      clear_fftw_array(m_fftw_input, m_Nx, m_Ny);
      set_real_part(H,
                    - m_load_density * m_standard_gravity * dt,
                    m_Mx, m_My, m_Nx, m_Ny, m_i0_offset, m_j0_offset,
                    m_fftw_input);
      fftw_execute(m_dft_forward);

      // Save fft2(-load_density * g * H * dt) in loadhat.
      copy_fftw_array(m_fftw_output, m_loadhat, m_Nx, m_Ny);
    }

    // Compute fft2(u).
    // no need to clear fftw_input: all values are overwritten
    {
      set_real_part(m_Uv, 1.0, m_Nx, m_Ny, m_Nx, m_Ny, 0, 0, m_fftw_input);
      fftw_execute(m_dft_forward);
    }

    // frhs = right.*fft2(uun) + fft2(dt*sszz);
    // uun1 = real(ifft2(frhs./left));
    {
      FFTWArray input(m_fftw_input, m_Nx, m_Ny),
        u_hat(m_fftw_output, m_Nx, m_Ny), load_hat(m_loadhat, m_Nx, m_Ny);
      for (int i = 0; i < m_Nx; i++) {
        for (int j = 0; j < m_Ny; j++) {
          const double
            C     = m_cx[i]*m_cx[i] + m_cy[j]*m_cy[j],
            part1 = 2.0 * m_eta * sqrt(C),
            part2 = (dt / 2.0) * (m_mantle_density * m_standard_gravity + m_D * C * C),
            A = part1 - part2,
            B = part1 + part2;

          input(i, j) = (load_hat(i, j) + A * u_hat(i, j)) / B;
        }
      }
    }

    fftw_execute(m_dft_inverse);
    get_real_part(m_fftw_output, 1.0 / (m_Nx * m_Ny), m_Nx, m_Ny, m_Nx, m_Ny, 0, 0, m_Uv);

    // Now tweak. (See the "correction" in section 5 of BuelerLingleBrown.)
    tweak(H, m_Uv, m_Nx, m_Ny);
  } else {
    // zero time step: viscous displacement is zero
    PetscErrorCode ierr = VecSet(m_Uv, 0.0); PISM_CHK(ierr, "VecSet");
  }

  // now compute elastic response if desired
  if (m_include_elastic) {
    compute_elastic_response(H, m_Ue);
  }

  update_displacement(m_Uv, m_Ue, m_U);
}

/*!
 * Compute elastic response to the load H
 *
 * @param[in] H load thickness (ice equivalent meters)
 * @param[out] dE elastic plate displacement
 */
void LingleClarkSerial::compute_elastic_response(petsc::Vec &H, petsc::Vec &dE) {

  // Compute fft2(load_density * H)
  //
  // Note that here the load is placed in the corner of the array on the extended grid
  // (offsets i0 and j0 are zero).
  {
    clear_fftw_array(m_fftw_input, m_Nx, m_Ny);
    set_real_part(H, m_load_density, m_Mx, m_My, m_Nx, m_Ny, 0, 0, m_fftw_input);
    fftw_execute(m_dft_forward);
  }

  // fft2(m_response_matrix) * fft2(load_density*H)
  //
  // Compute the product of Fourier transforms of the LRM and the load. This uses C++'s
  // native support for complex arithmetic.
  {
    FFTWArray
      input(m_fftw_input, m_Nx, m_Ny),
      LRM_hat(m_lrm_hat, m_Nx, m_Ny),
      load_hat(m_fftw_output, m_Nx, m_Ny);
    for (int i = 0; i < m_Nx; i++) {
      for (int j = 0; j < m_Ny; j++) {
        input(i, j) = LRM_hat(i, j) * load_hat(i, j);
      }
    }
  }

  // Compute the inverse transform and extract the elastic response.
  //
  // Here the offsets are:
  // i0 = m_Nx / 2,
  // j0 = m_Ny / 2.
  fftw_execute(m_dft_inverse);
  get_real_part(m_fftw_output, 1.0 / (m_Nx * m_Ny), m_Mx, m_My, m_Nx, m_Ny,
                m_Nx/2, m_Ny/2, dE);
}

/*! Compute total displacement by combining viscous and elastic contributions.
 *
 * @param[in] V viscous displacement
 * @param[in] dE elastic displacement
 * @param[out] dU total displacement
 */
void LingleClarkSerial::update_displacement(petsc::Vec &Uv, petsc::Vec &Ue, petsc::Vec &U) {
  // PISM grid
  petsc::VecArray2D
    u(U, m_Mx, m_My),
    u_elastic(Ue, m_Mx, m_My);
  // extended grid
  petsc::VecArray2D
    u_viscous(Uv, m_Nx, m_Ny, m_i0_offset, m_j0_offset);

  for (int i = 0; i < m_Mx; i++) {
    for (int j = 0; j < m_My; j++) {
      u(i, j) = u_viscous(i, j) + u_elastic(i, j);
    }
  }
}


/*!
 * Modify the plate displacement to correct for the effect of imposing periodic boundary conditions
 * at a finite distance.
 *
 * See Section 5 in [@ref BuelerLingleBrown].
 *
 * @param[in] load_thickness thickness of the load (used to compute the corresponding disc volume)
 * @param[in,out] U viscous plate displacement
 * @param[in] Nx grid size
 * @param[in] Ny grid size
 */
void LingleClarkSerial::tweak(petsc::Vec &load_thickness, petsc::Vec &U, int Nx, int Ny) {
  PetscErrorCode ierr = 0;
  petsc::VecArray2D u(U, Nx, Ny);

  // find average value along "distant" boundary of [-Lx, Lx]X[-Ly, Ly]
  // note domain is periodic, so think of cut locus of torus (!)
  // (will remove it:   uun1=uun1-(sum(uun1(1, :))+sum(uun1(:, 1)))/(2*N);)
  double average = 0.0;
  for (int i = 0; i < Nx; i++) {
    average += u(i, 0);
  }

  for (int j = 0; j < Ny; j++) {
    average += u(0, j);
  }

  average /= (double) (Nx + Ny);

  double shift = 0.0;

  {
    // tweak continued: replace far field with value for an equivalent disc load which has
    // R0=Lx*(2/3)=L/3 (instead of 1000km in MATLAB code: H0 = dx*dx*sum(sum(H))/(pi*1e6^2); %
    // trapezoid rule)
    const double L_average = (m_Lx + m_Ly) / 2.0;
    const double R         = L_average * (2.0 / 3.0);

    double H_sum = 0.0;
    ierr = VecSum(load_thickness, &H_sum); PISM_CHK(ierr, "VecSum");

    // compute disc thickness by dividing its volume by the area
    const double H = (H_sum * m_dx * m_dy) / (M_PI * R * R);

    shift = viscDisc(m_t_infty,                        // time in seconds
                     H,                                // disc thickness
                     R,                                // disc radius
                     L_average,                        // compute deflection at this radius
                     m_mantle_density, m_load_density, // mantle and load densities
                     m_standard_gravity,               //
                     m_D,                              // flexural rigidity
                     m_eta);                           // mantle viscosity
  }

  ierr = VecShift(U, shift - average); PISM_CHK(ierr, "VecShift");
}

} // end of namespace bed
} // end of namespace pism
