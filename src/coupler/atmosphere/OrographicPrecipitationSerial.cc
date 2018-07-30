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

#include "../../earth/matlablike.hh"
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
  m_gamma   = config.get_double("atmosphere.orographic_precipitation.lapse_rate");
  m_Theta_m   = config.get_double("atmosphere.orographic_precipitation.moist_adiabatic_lapse_rate");
  m_rho_Sref   = config.get_double("atmosphere.orographic_precipitation.reference_density");
  m_latitude   = config.get_double("atmosphere.orographic_precipitation.latitude");

  // derive more parameters
  m_Lx        = 0.5 * (m_Nx - 1.0) * m_dx;
  m_Ly        = 0.5 * (m_Ny - 1.0) * m_dy;
  m_Nxge      = m_Nx + 1;
  m_Nyge      = m_Ny + 1;
  m_i0_offset = (Nx - Mx) / 2;
  m_j0_offset = (Ny - My) / 2;

  // setup fftw stuff: FFTW builds "plans" based on observed performance
  m_fftw_input  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * m_Nx * m_Ny);
  m_fftw_output = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * m_Nx * m_Ny);
  m_Hhat     = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * m_Nx * m_Ny);
  m_Phat     = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * m_Nx * m_Ny);
  m_sigma     = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * m_Nx * m_Ny);

  m_eps = 1.0e-16;
  
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

  precompute_derived_constants();
}

OrographicPrecipitationSerial::~OrographicPrecipitationSerial() {
  fftw_destroy_plan(m_dft_forward);
  fftw_destroy_plan(m_dft_inverse);
  fftw_free(m_fftw_input);
  fftw_free(m_fftw_output);
  fftw_free(m_Hhat);
  fftw_free(m_Phat);
  fftw_free(m_sigma);
}

/**
 * Pre-compute coefficients used by the model.
 */
void OrographicPrecipitationSerial::precompute_coefficients() {

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

}

/**
 * Pre-compute derived constants.
 */
void OrographicPrecipitationSerial::precompute_derived_constants() {

  m_f = 2.0 * 7.2921e-5 * sin(m_latitude * M_PI / 180.0);

  m_u = -sin(m_wind_direction * 2.0 * M_PI / 360.0) * m_wind_speed;
  m_v = cos(m_wind_direction * 2.0 * M_PI / 360.0) * m_wind_speed;

  m_Cw = m_rho_Sref * m_Theta_m / m_gamma;
  
}

/**
 * Compute intrinsic frequency.
 */
void OrographicPrecipitationSerial::compute_intrinsic_frequency() {

  
}
  
/*!
 * Update precipitation.
 *
 * @param[in] H surface on the physical (Mx*My) grid
 */
void OrographicPrecipitationSerial::step(Vec H) {
  // solves:
  // Phat(k,l) = (Cw * i * sigma * Hhat(k,l)) /
  //             (1 - i * m * Hw) * (1 + i * sigma * tauc) * (1 + i * sigma * tauc);
  // see equation (49) in
  // R. B. Smith and I. Barstad, 2004:
  // A Linear Theory of Orographic Precipitation. J. Atmos. Sci. 61, 1377-1391.
  
  // Compute fft2(orography)
  {
    clear_fftw_input(m_fftw_input, m_Nx, m_Ny);
    set_fftw_input(H,
                   1.0,
                   m_Mx, m_My, m_i0_offset, m_j0_offset);
    fftw_execute(m_dft_forward);

    // Save fft2(orography) in Hhat.
    copy_fftw_output(m_fftw_output, m_Hhat, m_Nx, m_Ny);

  {
    VecAccessor2D<fftw_complex>
      sigma(m_sigma, m_Nx, m_Ny), sigma2(m_sigma2, m_Nx, m_Ny);
    for (int i = 0; i < m_Nx; i++) {
      for (int j = 0; j < m_Ny; j++) {
        sigma(i, j)[0] = m_u * m_cx[i] + m_v * m_cy[j];
        sigma(i, j)[1] = m_u * m_cx[i] + m_v * m_cy[j];
        // Regularization
        sigma2(i, j)[0] = pow(sigma(i,j)[0], 2.0);
        if ((fabs(sigma2(i, j)[0]) < m_eps) and fabs(sigma2(i, j)[0] >= 0.0)) {
            sigma2(i, j)[0] = m_eps;
          }
        if ((fabs(sigma2(i, j)[0]) < m_eps) and fabs(sigma2(i, j)[0] < 0.0)) {
            sigma2(i, j)[0] = -m_eps;
          }
        sigma2(i, j)[1] = pow(sigma(i,j)[1], 2.0);
        if ((fabs(sigma2(i, j)[1]) < m_eps) and fabs(sigma2(i, j)[1] >= 0.0)) {
            sigma2(i, j)[1] = m_eps;
          }
        if ((fabs(sigma2(i, j)[1]) < m_eps) and fabs(sigma2(i, j)[1] < 0.0)) {
            sigma2(i, j)[1] = -m_eps;
          }
      }
    }
  }

  // fftw_execute(m_dft_inverse);
  // get_fftw_output(m_Phat, 1.0 / (m_Nx * m_Ny), m_Nx, m_Ny, 0, 0);

  }

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
