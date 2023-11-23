// Copyright (C) 2018, 2019, 2020, 2021, 2023 Andy Aschwanden and Constantine Khroulev
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

#include "pism/coupler/atmosphere/OrographicPrecipitationSerial.hh"

#include <complex> // std::complex<double>, std::sqrt()
#include <gsl/gsl_math.h> // M_PI
#include <cassert>        // assert()
#include <cmath>          // std::exp()

#include "pism/util/ConfigInterface.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/fftw_utilities.hh"

namespace pism {
namespace atmosphere {

/*!
 * @param[in] config configuration database
 * @param[in] log logger
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
                                                             int Nx, int Ny)
  : m_Mx(Mx), m_My(My), m_Nx(Nx), m_Ny(Ny) {

  m_eps = 1.0e-18;

  // derive more parameters
  {
    m_i0_offset = (Nx - Mx) / 2;
    m_j0_offset = (Ny - My) / 2;

    m_kx = fftfreq(m_Nx, dx / (2.0 * M_PI));
    m_ky = fftfreq(m_Ny, dy / (2.0 * M_PI));
  }

  {
    m_background_precip_pre  = config.get_number("atmosphere.orographic_precipitation.background_precip_pre", "mm/s");
    m_background_precip_post = config.get_number("atmosphere.orographic_precipitation.background_precip_post", "mm/s");

    m_precip_scale_factor = config.get_number("atmosphere.orographic_precipitation.scale_factor");
    m_tau_c               = config.get_number("atmosphere.orographic_precipitation.conversion_time");
    m_tau_f               = config.get_number("atmosphere.orographic_precipitation.fallout_time");
    m_Hw                  = config.get_number("atmosphere.orographic_precipitation.water_vapor_scale_height");
    m_Nm                  = config.get_number("atmosphere.orographic_precipitation.moist_stability_frequency");
    m_wind_speed          = config.get_number("atmosphere.orographic_precipitation.wind_speed");
    m_wind_direction      = config.get_number("atmosphere.orographic_precipitation.wind_direction");
    m_gamma               = config.get_number("atmosphere.orographic_precipitation.lapse_rate");
    m_Theta_m             = config.get_number("atmosphere.orographic_precipitation.moist_adiabatic_lapse_rate");
    m_rho_Sref            = config.get_number("atmosphere.orographic_precipitation.reference_density");
    m_latitude            = config.get_number("atmosphere.orographic_precipitation.coriolis_latitude");
    m_truncate            = config.get_flag("atmosphere.orographic_precipitation.truncate");


    // derived constants
    m_f = 2.0 * 7.2921e-5 * sin(m_latitude * M_PI / 180.0);

    m_u = -sin(m_wind_direction * 2.0 * M_PI / 360.0) * m_wind_speed;
    m_v = -cos(m_wind_direction * 2.0 * M_PI / 360.0) * m_wind_speed;

    m_Cw = m_rho_Sref * m_Theta_m / m_gamma;
  }

  // memory allocation
  {
    PetscErrorCode ierr = 0;

    // precipitation
    ierr = VecCreateSeq(PETSC_COMM_SELF, m_Mx * m_My, m_precipitation.rawptr());
    PISM_CHK(ierr, "VecCreateSeq");

    // FFTW arrays
    m_fftw_input  = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * m_Nx * m_Ny);
    m_fftw_output = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * m_Nx * m_Ny);
    m_G_hat       = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * m_Nx *    m_Ny);

    // FFTW plans
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
  }

  // initialize the Gaussian filter
  {
    FFTWArray G_hat(m_G_hat, m_Nx, m_Ny);
    double sigma = config.get_number("atmosphere.orographic_precipitation.smoothing_standard_deviation");

    if (sigma > 0.0) {
      FFTWArray
        fftw_output(m_fftw_output, m_Nx, m_Ny),
        fftw_input(m_fftw_input, m_Nx, m_Ny);

      int
        Nx2 = Nx / 2,
        Ny2 = Ny / 2;

      double sum = 0.0;
      for (int i = 0; i < m_Nx; i++) {
        for (int j = 0; j < m_Ny; j++) {
          int
            p = i <= Nx2 ? i : m_Nx - i,
            q = j <= Ny2 ? j : m_Ny - j;
          double
            x = p * dx,
            y = q * dy;

          double G = std::exp(-0.5 * (x * x + y * y) / (sigma * sigma));
          sum += G;

          fftw_input(i, j) = G;
        }
      }

      // normalize:
      assert(sum > 0.0);
      for (int i = 0; i < m_Nx; i++) {
        for (int j = 0; j < m_Ny; j++) {
          fftw_input(i, j) /= sum;
        }
      }

      // compute FFT of the Gaussian
      fftw_execute(m_dft_forward);

      // copy to m_G_hat
      for (int i = 0; i < m_Nx; i++) {
        for (int j = 0; j < m_Ny; j++) {
          G_hat(i, j) = fftw_output(i, j);
        }
      }

    } else {
      // fill m_G_hat with ones to disable smoothing
      for (int i = 0; i < m_Nx; i++) {
        for (int j = 0; j < m_Ny; j++) {
          G_hat(i, j) = 1.0;
        }
      }
    }
  }
}

OrographicPrecipitationSerial::~OrographicPrecipitationSerial() {
  fftw_destroy_plan(m_dft_forward);
  fftw_destroy_plan(m_dft_inverse);
  fftw_free(m_fftw_input);
  fftw_free(m_fftw_output);
  fftw_free(m_G_hat);
}

/*!
 * Return precipitation (FIXME: units?)
 */
Vec OrographicPrecipitationSerial::precipitation() const {
  return m_precipitation;
}

/*!
 * Update precipitation.
 *
 * @param[in] surface_elevation surface on the physical (Mx*My) grid
 */
void OrographicPrecipitationSerial::update(petsc::Vec &surface_elevation) {
  // solves:
  // Phat(k,l) = (Cw * i * sigma * Hhat(k,l)) /
  //             (1 - i * m * Hw) * (1 + i * sigma * tauc) * (1 + i * sigma * tauf);
  // see equation (49) in
  // R. B. Smith and I. Barstad, 2004:
  // A Linear Theory of Orographic Precipitation. J. Atmos. Sci. 61, 1377-1391.

  std::complex<double> I(0.0, 1.0);

  // Compute fft2(surface_elevation)
  {
    clear_fftw_array(m_fftw_input, m_Nx, m_Ny);
    set_real_part(surface_elevation,
                  1.0,
                  m_Mx, m_My,
                  m_Nx, m_Ny,
                  m_i0_offset, m_j0_offset,
                  m_fftw_input);
    fftw_execute(m_dft_forward);
  }

  {
    FFTWArray
      fftw_output(m_fftw_output, m_Nx, m_Ny),
      fftw_input(m_fftw_input, m_Nx, m_Ny),
      G_hat(m_G_hat, m_Nx, m_Ny);

    for (int i = 0; i < m_Nx; i++) {
      const double kx = m_kx[i];
      for (int j = 0; j < m_Ny; j++) {
        const double ky = m_ky[j];

        // FFT(h) * FFT(Gaussian), i.e. FFT(smoothed ice surface elevation)
        const auto &h_hat = fftw_output(i, j) * G_hat(i, j);

        double sigma = m_u * kx + m_v * ky;

        // See equation (6) in [@ref SmithBarstadBonneau2005]
        std::complex<double> m;
        {
          double denominator = sigma * sigma - m_f * m_f;

          // avoid dividing by zero:
          if (fabs(denominator) < m_eps) {
            denominator = denominator >= 0 ? m_eps : -m_eps;
          }

          double m_squared = (m_Nm * m_Nm - sigma * sigma) * (kx * kx + ky * ky) / denominator;

          // Note: this is a *complex* square root.
          m = std::sqrt(std::complex<double>(m_squared));

          if (m_squared >= 0.0 and sigma != 0.0) {
            m *= sigma > 0.0 ? 1.0 : -1.0;
          }
        }

        // avoid dividing by zero:
        double delta = 0.0;
        if (std::abs(1.0 - I * m * m_Hw) < m_eps) {
          delta = m_eps;
        }

        // See equation (49) in [@ref SmithBarstad2004] or equation (3) in [@ref
        // SmithBarstadBonneau2005].
        auto P_hat = h_hat * (m_Cw * I * sigma /
                              ((1.0 - I * m * m_Hw + delta) *
                               (1.0 + I * sigma * m_tau_c) *
                               (1.0 + I * sigma * m_tau_f)));
        // Note: sigma, m_tau_c, and m_tau_f are purely real, so the second and the third
        // factors in the denominator are never zero.
        //
        // The first factor (1 - i m H_w) *could* be zero. Here we check if it is and
        // "regularize" if necessary.

        fftw_input(i, j) = P_hat;
      }
    }
  }

  fftw_execute(m_dft_inverse);

  // get m_fftw_output and put it into m_precipitation
  get_real_part(m_fftw_output,
                1.0 / (m_Nx * m_Ny),
                m_Mx, m_My,
                m_Nx, m_Ny,
                m_i0_offset, m_j0_offset,
                m_precipitation);

  petsc::VecArray2D p(m_precipitation, m_Mx, m_My);
  for (int i = 0; i < m_Mx; i++) {
    for (int j = 0; j < m_My; j++) {
      p(i, j) += m_background_precip_pre;
      if (m_truncate) {
        p(i, j) = std::max(p(i, j), 0.0);
      }
      p(i, j) *= m_precip_scale_factor;
      p(i, j) += m_background_precip_post;
    }
  }
}

} // end of namespace atmosphere
} // end of namespace pism
