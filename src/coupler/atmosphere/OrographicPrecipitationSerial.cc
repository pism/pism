// Copyright (C) 2018 Andy Aschwanden and Constantine Khroulev
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

#include <cmath>                // sqrt
#include <fftw3.h>
#include <gsl/gsl_math.h>       // M_PI

#include "matlablike.hh"
#include "greens.hh"
#include "OrographicPrecipitationSerial.hh"

#include "pism/util/pism_utilities.hh"
#include "pism/util/ConfigInterface.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/petscwrappers/Vec.hh"

namespace pism {
namespace atmosphere {

template <class T>
class VecAccessor2D {
public:
  VecAccessor2D(T* a, int Mx, int My, int i_offset, int j_offset)
    : m_Mx(Mx), m_My(My), m_i_offset(i_offset), m_j_offset(j_offset), m_array(a) {}

  VecAccessor2D(T* a, int Mx, int My)
    : m_Mx(Mx), m_My(My), m_i_offset(0), m_j_offset(0), m_array(a) {}

  inline T& operator()(int i, int j) {
    return m_array[(m_j_offset + j) + m_My * (m_i_offset + i)];
  }

private:
  int m_Mx, m_My, m_i_offset, m_j_offset;
  T* m_array;
};

//! \brief Fill `input` with zeros.
static void clear_fftw_input(fftw_complex *input, int Nx, int Ny) {
  VecAccessor2D<fftw_complex> fftw_in(input, Nx, Ny);
  for (int i = 0; i < Nx; ++i) {
    for (int j = 0; j < Ny; ++j) {
      fftw_in(i, j)[0] = 0;
      fftw_in(i, j)[1] = 0;
    }
  }
}

//! @brief Copy `source` to `destination`.
static void copy_fftw_output(fftw_complex *source, fftw_complex *destination,
                             int Nx, int Ny) {
  VecAccessor2D<fftw_complex> S(source, Nx, Ny), D(destination, Nx, Ny);
  for (int i = 0; i < Nx; ++i) {
    for (int j = 0; j < Ny; ++j) {
      D(i, j)[0] = S(i, j)[0];
      D(i, j)[1] = S(i, j)[1];
    }
  }
}

/*!
 * @param[in] config configuration database
 * @param[in] Mx grid size in the X direction
 * @param[in] My grid size in the Y direction
 * @param[in] dx grid spacing in the X direction
 * @param[in] dy grid spacing in the Y direction
 * @param[in] Nx extended grid size in the X direction
 * @param[in] Ny extended grid size in the Y direction
 */
OrographicPrecipitationSerial::OrographicPrecipitationSerial(const Config &config,
                                     int Mx, int My,
                                     double dx, double dy,
                                     int Nx, int Ny) {

  // set parameters
  bool m_include_elastic = 0;

  // grid parameters
  m_Mx = Mx;
  m_My = My;
  m_dx = dx;
  m_dy = dy;
  m_Nx = Nx;
  m_Ny = Ny;

  m_tau_c   = config.get_double("atmosphere.orographic_precipitation.conversion_time");
  m_tau_f   = config.get_double("atmosphere.orographic_precipitation.fallout_time");
  m_Hw   = config.get_double("atmosphere.orographic_precipitation.water_vapor_scale_height");
  m_Nm   = config.get_double("atmosphere.orographic_precipitation.moist_stability");
  m_wind_speed   = config.get_double("atmosphere.orographic_precipitation.wind_speed");
  m_wind_direction   = config.get_double("atmosphere.orographic_precipitation.wind_direction");
  
  m_load_density   = config.get_double("constants.ice.density");
  m_mantle_density = config.get_double("atmosphere_deformation.mantle_density");
  m_eta            = config.get_double("atmosphere_deformation.mantle_viscosity");
  m_D              = config.get_double("atmosphere_deformation.lithosphere_flexural_rigidity");

  m_standard_gravity = config.get_double("constants.standard_gravity");

  // derive more parameters
  m_Lx        = 0.5 * (m_Nx - 1.0) * m_dx;
  m_Ly        = 0.5 * (m_Ny - 1.0) * m_dy;
  m_Nxge      = m_Nx + 1;
  m_Nyge      = m_Ny + 1;
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

  // elastic load response matrix
  ierr = VecCreateSeq(PETSC_COMM_SELF, m_Nxge * m_Nyge, m_load_response_matrix.rawptr());
  PISM_CHK(ierr, "VecCreateSeq");

  // setup fftw stuff: FFTW builds "plans" based on observed performance
  m_fftw_input  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * m_Nx * m_Ny);
  m_fftw_output = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * m_Nx * m_Ny);
  m_loadhat     = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * m_Nx * m_Ny);

  // fill m_fftw_input with zeros
  {
    VecAccessor2D<fftw_complex> tmp(m_fftw_input, m_Nx, m_Ny);
    for (int j = 0; j < m_Ny; j++) {
      for (int i = 0; i < m_Nx; i++) {
        tmp(i, j)[0] = 0.0;
        tmp(i, j)[1] = 0.0;
      }
    }
  }

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

OrographicPrecipitationSerial::~OrographicPrecipitationSerial() {
  fftw_destroy_plan(m_dft_forward);
  fftw_destroy_plan(m_dft_inverse);
  fftw_free(m_fftw_input);
  fftw_free(m_fftw_output);
  fftw_free(m_loadhat);
}

/*!
 * Return total displacement.
 */
Vec OrographicPrecipitationSerial::total_displacement() const {
  return m_U;
}

/*!
 * Return viscous plate displacement.
 */
Vec OrographicPrecipitationSerial::viscous_displacement() const {
  return m_Uv;
}

/**
 * Pre-compute coefficients used by the model.
 */
void OrographicPrecipitationSerial::precompute_coefficients() {
  PetscErrorCode ierr = 0;

  m_cx.resize(m_Nx);
  m_cy.resize(m_Ny);

  // Coefficients for Fourier spectral method Laplacian
  // MATLAB version:  cx=(pi/Lx)*[0:Nx/2 Nx/2-1:-1:1]
  for (int i = 0; i <= m_Nx / 2; i++) {
    m_cx[i] = (M_PI / m_Lx) * i;
  }

  for (int i = m_Nx / 2 + 1; i < m_Nx; i++) {
    m_cx[i] = (M_PI / m_Lx) * (m_Nx - i);
  }

  for (int j = 0; j <= m_Ny / 2; j++) {
    m_cy[j] = (M_PI / m_Ly) * j;
  }

  for (int j = m_Ny / 2 + 1; j < m_Ny; j++) {
    m_cy[j] = (M_PI / m_Ly) * (m_Ny - j);
  }

  // compare geforconv.m
  if (m_include_elastic) {
    ierr = PetscPrintf(PETSC_COMM_SELF,
                       "     computing spherical elastic load response matrix ...");
    PISM_CHK(ierr, "PetscPrintf");

    petsc::VecArray2D II(m_load_response_matrix, m_Nxge, m_Nyge);
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

/*!
 * Solve
 *
 * @f$ 2 \nu |\nabla| \diff{u}{t} + \rho_r g U + D\nabla^4 U = \sigma_{zz}@f$
 *
 * for @f$ U @f$, treating @f$ \diff{u}{t} @f$ and @f$ \sigma_{zz} @f$ as known.
 *
 * @param[in] load_thickness load thickness, meters
 * @param[in] atmosphere_uplift atmosphere uplift, m/second
 *
 * Here `load_thickness` is used to compute the load @f$ \sigma_{zz} @f$ and `atmosphere_uplift` is
 * @f$ \diff{u}{t} @f$ itself.
 *
 */
void OrographicPrecipitationSerial::uplift_problem(Vec load_thickness, Vec atmosphere_uplift,
                                       Vec output) {

  // Compute fft2(-load_density * g * load_thickness)
  {
    clear_fftw_input(m_fftw_input, m_Nx, m_Ny);
    set_fftw_input(load_thickness, - m_load_density * m_standard_gravity, m_Mx, m_My, m_i0_offset, m_j0_offset);
    fftw_execute(m_dft_forward);
    // Save fft2(-load_density * g * load_thickness) in loadhat.
    copy_fftw_output(m_fftw_output, m_loadhat, m_Nx, m_Ny);
  }

  // fft2(uplift)
  {
    clear_fftw_input(m_fftw_input, m_Nx, m_Ny);
    set_fftw_input(atmosphere_uplift, 1.0, m_Mx, m_My, m_i0_offset, m_j0_offset);
    fftw_execute(m_dft_forward);
  }

  {
    VecAccessor2D<fftw_complex>
      u0_hat(m_fftw_input, m_Nx, m_Ny),
      load_hat(m_loadhat, m_Nx, m_Ny),
      uplift_hat(m_fftw_output, m_Nx, m_Ny);

    for (int i = 0; i < m_Nx; i++) {
      for (int j = 0; j < m_Ny; j++) {
        const double
          C = m_cx[i]*m_cx[i] + m_cy[j]*m_cy[j],
          A = - 2.0 * m_eta * sqrt(C),
          B = m_mantle_density * m_standard_gravity + m_D * C * C;

        u0_hat(i, j)[0] = (load_hat(i, j)[0] + A * uplift_hat(i, j)[0]) / B;
        u0_hat(i, j)[1] = (load_hat(i, j)[1] + A * uplift_hat(i, j)[1]) / B;
      }
    }
  }

  fftw_execute(m_dft_inverse);
  get_fftw_output(output, 1.0 / (m_Nx * m_Ny), m_Nx, m_Ny, 0, 0);

  tweak(load_thickness, output, m_Nx, m_Ny, 0.0);
}

/*! Initialize using provided load thickness and the atmosphere uplift rate.
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
 * @param[in] uplift initial atmosphere uplift on the PISM grid
 *
 * Sets m_Uv, m_Ue, m_U.
 */
void OrographicPrecipitationSerial::bootstrap(Vec thickness, Vec uplift) {

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
 * @param[in] V initial viscous plate displacement on the extended grid
 *
 * Sets m_Uv, m_Ue, m_U.
 */
void OrographicPrecipitationSerial::init(Vec thickness, Vec viscous_displacement) {
  PetscErrorCode ierr = 0;

  ierr = VecCopy(viscous_displacement, m_Uv);  PISM_CHK(ierr, "VecCopy");

  if (m_include_elastic) {
    compute_elastic_response(thickness, m_Ue);
  } else {
    ierr = VecSet(m_Ue, 0.0); PISM_CHK(ierr, "VecSet");
  }

  update_displacement(m_Uv, m_Ue, m_U);
}

/*!
 * Perform a time step.
 *
 * @param[in] dt_seconds time step length
 * @param[in] H load thickness on the physical (Mx*My) grid
 */
void OrographicPrecipitationSerial::step(double dt_seconds, Vec H) {
  // solves:
  //     (2 eta |grad| U^{n+1}) + (dt/2) * (rho_r g U^{n+1} + D grad^4 U^{n+1})
  //   = (2 eta |grad| U^n) - (dt/2) * (rho_r g U^n + D grad^4 U^n) - dt * rho g H_start
  // where U=plate displacement; see equation (7) in
  // Bueler, Lingle, Brown (2007) "Fast computation of a viscoelastic
  // deformable Earth model for ice sheet simulations", Ann. Glaciol. 46, 97--105

  // Compute fft2(-load_density * g * dt * H)
  {
    clear_fftw_input(m_fftw_input, m_Nx, m_Ny);
    set_fftw_input(H,
                   - m_load_density * m_standard_gravity * dt_seconds,
                   m_Mx, m_My, m_i0_offset, m_j0_offset);
    fftw_execute(m_dft_forward);

    // Save fft2(-load_density * g * H * dt) in loadhat.
    copy_fftw_output(m_fftw_output, m_loadhat, m_Nx, m_Ny);
  }

  // Compute fft2(u).
  // no need to clear fftw_input: all values are overwritten
  {
    set_fftw_input(m_Uv, 1.0, m_Nx, m_Ny, 0, 0);
    fftw_execute(m_dft_forward);
  }

  // frhs = right.*fft2(uun) + fft2(dt*sszz);
  // uun1 = real(ifft2(frhs./left));
  {
    VecAccessor2D<fftw_complex> input(m_fftw_input, m_Nx, m_Ny),
      u_hat(m_fftw_output, m_Nx, m_Ny), load_hat(m_loadhat, m_Nx, m_Ny);
    for (int i = 0; i < m_Nx; i++) {
      for (int j = 0; j < m_Ny; j++) {
        const double
          C     = m_cx[i]*m_cx[i] + m_cy[j]*m_cy[j],
          part1 = 2.0 * m_eta * sqrt(C),
          part2 = (dt_seconds / 2.0) * (m_mantle_density * m_standard_gravity + m_D * C * C),
          A = part1 - part2,
          B = part1 + part2;

        input(i, j)[0] = (load_hat(i, j)[0] + A * u_hat(i, j)[0]) / B;
        input(i, j)[1] = (load_hat(i, j)[1] + A * u_hat(i, j)[1]) / B;
      }
    }
  }

  fftw_execute(m_dft_inverse);
  get_fftw_output(m_Uv, 1.0 / (m_Nx * m_Ny), m_Nx, m_Ny, 0, 0);

  // Now tweak. (See the "correction" in section 5 of BuelerLingleBrown.)
  //
  // Here 1e16 approximates t = \infty.
  tweak(H, m_Uv, m_Nx, m_Ny, 1e16);

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
void OrographicPrecipitationSerial::compute_elastic_response(Vec H, Vec dE) {

  conv2_same(H, m_Mx, m_My, m_load_response_matrix, m_Nxge, m_Nyge, dE);

  PetscErrorCode ierr = VecScale(m_Ue, m_load_density); PISM_CHK(ierr, "VecScale");
}

/*! Compute total displacement by combining viscous and elastic contributions.
 *
 * @param[in] V viscous displacement
 * @param[in] dE elastic displacement
 * @param[out] dU total displacement
 */
void OrographicPrecipitationSerial::update_displacement(Vec Uv, Vec Ue, Vec U) {
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
 * @param[in] time time, seconds (usually 0 or a large number approximating \infty)
 */
void OrographicPrecipitationSerial::tweak(Vec load_thickness, Vec U, int Nx, int Ny, double time) {
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

  if (time > 0.0) {
    // tweak continued: replace far field with value for an equivalent disc load which has
    // R0=Lx*(2/3)=L/3 (instead of 1000km in MATLAB code: H0 = dx*dx*sum(sum(H))/(pi*1e6^2); %
    // trapezoid rule)
    const double L_average = (m_Lx + m_Ly) / 2.0;
    const double R         = L_average * (2.0 / 3.0);

    double H_sum = 0.0;
    ierr = VecSum(load_thickness, &H_sum); PISM_CHK(ierr, "VecSum");

    // compute disc thickness by dividing its volume by the area
    const double H = (H_sum * m_dx * m_dy) / (M_PI * R * R);

    shift = viscDisc(time,               // time in seconds
                     H,                  // disc thickness
                     R,                  // disc radius
                     L_average,          // compute deflection at this radius
                     m_mantle_density, m_load_density,    // mantle and load densities
                     m_standard_gravity, //
                     m_D,                // flexural rigidity
                     m_eta);             // mantle viscosity
  }

  ierr = VecShift(U, shift - average); PISM_CHK(ierr, "VecShift");
}

//! \brief Set the real part of fftw_input to vec_input.
/*!
 * Sets the imaginary part to zero.
 */
void OrographicPrecipitationSerial::set_fftw_input(Vec vec_input, double normalization,
                                       int Mx, int My, int i0, int j0) {
  petsc::VecArray2D in(vec_input, Mx, My);
  VecAccessor2D<fftw_complex> input(m_fftw_input, m_Nx, m_Ny, i0, j0);
  for (int j = 0; j < My; ++j) {
    for (int i = 0; i < Mx; ++i) {
      input(i, j)[0] = in(i, j) * normalization;
      input(i, j)[1] = 0.0;
    }
  }
}

//! \brief Get the real part of fftw_output and put it in output.
void OrographicPrecipitationSerial::get_fftw_output(Vec output, double normalization,
                                        int Mx, int My, int i0, int j0) {
  petsc::VecArray2D out(output, Mx, My);
  VecAccessor2D<fftw_complex> fftw_out(m_fftw_output, m_Nx, m_Ny, i0, j0);
  for (int j = 0; j < My; ++j) {
    for (int i = 0; i < Mx; ++i) {
      out(i, j) = fftw_out(i, j)[0] * normalization;
    }
  }
}

} // end of namespace atmosphere
} // end of namespace pism


    // eps = 1e-18
    // pad = 250

    // ny, nx = orography.shape
    // logger.debug('Raster shape before padding ({},{})'.format(nx, ny))

    // padded_orography = np.pad(orography, pad, 'constant')
    // nx, ny = padded_orography.shape
    // logger.debug('Raster shape after padding ({},{})'.format(ny, nx))

    // logger.info('Fourier transform orography')
    // padded_orography_fft = np.fft.fft2(padded_orography)

    // x_n_value = np.fft.fftfreq(ny, (1.0 / ny))
    // y_n_value = np.fft.fftfreq(nx, (1.0 / nx))

    // x_len = nx * dx
    // y_len = ny * dy
    // kx_line = np.divide(np.multiply(2.0 * np.pi, x_n_value), x_len)
    // ky_line = np.divide(
    //     np.multiply(
    //         2.0 * np.pi,
    //         y_n_value),
    //     y_len)[
    //     np.newaxis].T

    // kx = np.tile(kx_line, (nx, 1))
    // ky = np.tile(ky_line, (1, ny))

    // # Intrinsic frequency sigma = U*k + V*l
    // u0 = constants.u
    // v0 = constants.v

    // logger.info('Calculate sigma')
    // sigma = np.add(np.multiply(kx, u0), np.multiply(ky, v0))
    // sigma_sqr_reg = sigma ** 2
    // m_denom = np.power(sigma, 2.) - constants.f**2

    // sigma_sqr_reg[
    //     np.logical_and(
    //         np.fabs(sigma_sqr_reg) < eps,
    //         np.fabs(
    //             sigma_sqr_reg >= 0))] = eps
    // sigma_sqr_reg[
    //     np.logical_and(
    //         np.fabs(sigma_sqr_reg) < eps,
    //         np.fabs(
    //             sigma_sqr_reg < 0))] = -eps

    // # The vertical wave number
    // # Eqn. 12
    // # Regularization
    // m_denom[
    //     np.logical_and(
    //         np.fabs(m_denom) < eps,
    //         np.fabs(m_denom) >= 0)] = eps
    // m_denom[
    //     np.logical_and(
    //         np.fabs(m_denom) < eps,
    //         np.fabs(m_denom) < 0)] = -eps

    // m1 = np.divide(
    //     np.subtract(
    //         constants.Nm**2,
    //         np.power(
    //             sigma,
    //             2.)),
    //     m_denom)
    // m2 = np.add(np.power(kx, 2.), np.power(ky, 2.))
    // m_sqr = np.multiply(m1, m2)
    // logger.info('Calculating m')
    // m = np.sqrt(-1 * m_sqr)
    // # Regularization
    // m[np.logical_and(m_sqr >= 0, sigma == 0)] = np.sqrt(
    //     m_sqr[np.logical_and(m_sqr >= 0, sigma == 0)])
    // m[np.logical_and(m_sqr >= 0, sigma != 0)] = np.sqrt(m_sqr[np.logical_and(
    //     m_sqr >= 0, sigma != 0)]) * np.sign(sigma[np.logical_and(m_sqr >= 0, sigma != 0)])
    // # Numerator in Eqn. 49
    // P_karot_num = np.multiply(np.multiply(np.multiply(
    //     constants.Cw, 1j), sigma), padded_orography_fft)
    // P_karot_denom_Hw = np.subtract(1, np.multiply(
    //     np.multiply(constants.Hw, m), 1j))
    // P_karot_denom_tauc = np.add(1, np.multiply(np.multiply(
    //     sigma, constants.tau_c), 1j))
    // P_karot_denom_tauf = np.add(1, np.multiply(np.multiply(
    //     sigma, constants.tau_f), 1j))
    // # Denominator in Eqn. 49
    // P_karot_denom = np.multiply(
    //     P_karot_denom_Hw, np.multiply(
    //         P_karot_denom_tauc, P_karot_denom_tauf))
    // P_karot = np.divide(P_karot_num, P_karot_denom)

    // # Converting from wave domain back to space domain
    // logger.info('Performing inverse Fourier transform')
    // P = np.fft.ifft2(P_karot)
    // spy = 31556925.9747
    // logger.info('De-pad array')
    // P = P[pad:-pad, pad:-pad]
    // P = np.multiply(np.real(P), 3600)   # mm hr-1
    // # Add background precip
    // P0 = constants.P0
    // logger.info('Adding background precpipitation {} mm hr-1'.format(P0))
    // P += P0
    // # Truncation

    // if truncate:
    //     logger.info('Truncate precipitation')
    //     P[P < 0] = 0
    // P_scale = constants.P_scale
    // logger.info('Scale precipitation P = P * {}'.format(P_scale))
    // P *= P_scale
