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

#include "OrographicPrecipitationSerial.hh"
#include "../../earth/matlablike.hh"
#include <cmath>   // sqrt
#include <complex> // I
#include <fftw3.h>
#include <gsl/gsl_math.h> // M_PI
#include <iostream>

#include "pism/util/ConfigInterface.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/petscwrappers/Vec.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/util/fftw_utilities.hh"

namespace pism {
namespace atmosphere {

template <class T>
class VecAccessor2D {
public:
  VecAccessor2D(T *a, int Mx, int My, int i_offset, int j_offset)
    : m_Mx(Mx), m_My(My), m_i_offset(i_offset), m_j_offset(j_offset), m_array(a) {
  }

  VecAccessor2D(T *a, int Mx, int My)
    : m_Mx(Mx), m_My(My), m_i_offset(0), m_j_offset(0), m_array(a) {
  }

  inline T &operator()(int i, int j) {
    return m_array[(m_j_offset + j) + m_My * (m_i_offset + i)];
  }

private:
  int m_Mx, m_My, m_i_offset, m_j_offset;
  T *m_array;
};
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
                                                             const Logger::ConstPtr log,
                                                             int Mx, int My,
                                                             double dx, double dy,
                                                             int Nx, int Ny)
  : m_Mx(Mx), m_My(My), m_dx(dx), m_dy(dy), m_Nx(Nx), m_Ny(Ny), m_log(log) {

  // derive more parameters
  m_Lx        = 0.5 * (m_Nx - 1.0) * m_dx;
  m_Ly        = 0.5 * (m_Ny - 1.0) * m_dy;
  m_Nxge      = m_Nx + 1;
  m_Nyge      = m_Ny + 1;
  m_i0_offset = (Nx - Mx) / 2;
  m_j0_offset = (Ny - My) / 2;

  // FIXME: use these!
  // constants
  // double precip_pre = config.get_double("atmosphere.orographic_precipitation.background_precip_pre", "m/s");
  // double precip_post = config.get_double("atmosphere.orographic_precipitation.background_precip_post", "m/s");

  m_precip_scale_factor = config.get_double("atmosphere.orographic_precipitation.scale_factor");
  m_tau_c               = config.get_double("atmosphere.orographic_precipitation.conversion_time");
  m_tau_f               = config.get_double("atmosphere.orographic_precipitation.fallout_time");
  m_Hw                  = config.get_double("atmosphere.orographic_precipitation.water_vapor_scale_height");
  m_Nm                  = config.get_double("atmosphere.orographic_precipitation.moist_stability_frequency");
  m_wind_speed          = config.get_double("atmosphere.orographic_precipitation.wind_speed");
  m_wind_direction      = config.get_double("atmosphere.orographic_precipitation.wind_direction");
  m_gamma               = config.get_double("atmosphere.orographic_precipitation.lapse_rate");
  m_Theta_m             = config.get_double("atmosphere.orographic_precipitation.moist_adiabatic_lapse_rate");
  m_rho_Sref            = config.get_double("atmosphere.orographic_precipitation.reference_density");
  m_latitude            = config.get_double("atmosphere.orographic_precipitation.coriolis_latitude");
  m_truncate            = config.get_boolean("atmosphere.orographic_precipitation.truncate");

  // memory allocation
  PetscErrorCode ierr = 0;

  // precipitation
  ierr = VecCreateSeq(PETSC_COMM_SELF, m_Mx * m_My, m_p.rawptr());
  PISM_CHK(ierr, "VecCreateSeq");
  ierr = VecCopy(orographic_precipitation(), m_p);
  PISM_CHK(ierr, "VecCopy");

  // setup fftw stuff: FFTW builds "plans" based on observed performance
  m_fftw_input  = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * m_Nx * m_Ny);
  m_fftw_output = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * m_Nx * m_Ny);
  m_Hhat        = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * m_Nx * m_Ny);
  m_Phat        = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * m_Nx * m_Ny);
  m_m           = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * m_Nx * m_Ny);
  m_sigma       = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * m_Nx * m_Ny);

  m_eps = 1.0e-18;

  // fill m_fftw_input with zeros
  clear_fftw_array(m_fftw_output, m_Nx, m_Ny);

  // fill m_sigma with zeros
  clear_fftw_array(m_sigma, m_Nx, m_Ny);

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

  // Coefficients for Fourier spectral method Laplacian
  // MATLAB version:  cx=(pi/Lx)*[0:Nx/2 Nx/2-1:-1:1]
  m_cx = fftfreq(m_Nx, M_PI / m_Lx);
  m_cy = fftfreq(m_Ny, M_PI / m_Ly);

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
  fftw_free(m_m);
}

/*!
 * Return viscous plate displacement.
 */
Vec OrographicPrecipitationSerial::orographic_precipitation() const {
  return m_p;
}

/**
 * Pre-compute derived constants.
 */
void OrographicPrecipitationSerial::precompute_derived_constants() {

  m_log->message(4, "* Precompute derived constants\n");

  m_f = 2.0 * 7.2921e-5 * sin(m_latitude * M_PI / 180.0);

  m_u = -sin(m_wind_direction * 2.0 * M_PI / 360.0) * m_wind_speed;
  m_v = cos(m_wind_direction * 2.0 * M_PI / 360.0) * m_wind_speed;

  m_Cw = m_rho_Sref * m_Theta_m / m_gamma;
}

void OrographicPrecipitationSerial::compute_intrinsic_frequency() {

  m_log->message(4, "* Compute intrinsic frequency\n");

  {
    VecAccessor2D<fftw_complex> sigma(m_sigma, m_Nx, m_Ny);

    for (int i = 0; i < m_Nx; i++) {
      for (int j = 0; j < m_Ny; j++) {
        sigma(i, j)[0] = m_u * m_cx[i] + m_v * m_cy[j];
        sigma(i, j)[1] = 0.0;
      }
    }
  }
}

void OrographicPrecipitationSerial::compute_vertical_wave_number() {
  // Computes:
  // m = [ ((Nm^2 - sigma^2) / sigma^2) * (k^2 + l^2) ]^(1/2)

  m_log->message(4, "* Compute vertical wave number\n");

  {
    VecAccessor2D<fftw_complex> m(m_m, m_Nx, m_Ny), sigma(m_sigma, m_Nx, m_Ny);

    for (int i = 0; i < m_Nx; i++) {
      for (int j = 0; j < m_Ny; j++) {
        double sigma2_0 = pow(sigma(i, j)[0], 2.0);
        double sigma2_1 = pow(sigma(i, j)[1], 2.0);

        double nom_0 = pow(m_Nm, 2.0) - sigma2_0;
        double nom_1 = pow(m_Nm, 2.0) - sigma2_1;

        double denom_0 = sigma2_0;
        double denom_1 = sigma2_1;

        // regularization
        if (fabs(sigma2_0) < m_eps and fabs(sigma2_0 >= 0)) {
          denom_0 = m_eps;
        }
        if (fabs(sigma2_0) < m_eps and fabs(sigma2_0 < 0)) {
          denom_0 = -m_eps;
        }
        if (fabs(sigma2_1) < m_eps and fabs(sigma2_1 >= 0)) {
          denom_1 = m_eps;
        }
        if (fabs(sigma2_1) < m_eps and fabs(sigma2_1 < 0)) {
          denom_1 = -m_eps;
        }

        m(i, j)[0] = pow(nom_0 / denom_0 * (m_cx[i] * m_cx[i] + m_cy[j] * m_cy[j]), 0.5);
        m(i, j)[1] = pow(nom_1 / denom_1 * (m_cx[i] * m_cx[i] + m_cy[j] * m_cy[j]), 0.5);
      }
    }
  }
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
    clear_fftw_array(m_fftw_input, m_Nx, m_Ny);
    set_real_part(H,
                  1.0,
                  m_Mx, m_My, m_Nx, m_Ny,
                  m_i0_offset, m_j0_offset,
                  m_fftw_input);
    fftw_execute(m_dft_forward);

    // Save fft2(orography) in Hhat.
    copy_fftw_array(m_fftw_output, m_Hhat, m_Nx, m_Ny);
  }

  compute_intrinsic_frequency();
  compute_vertical_wave_number();

  std::complex<double> I(0.0, 1.0);

  {
    VecAccessor2D<fftw_complex>
      m(m_m, m_Nx, m_Ny),
      Hhat(m_Hhat, m_Nx, m_Ny),
      Phat(m_Phat, m_Nx, m_Ny),
      sigma(m_sigma, m_Nx, m_Ny);

    for (int i = 0; i < m_Nx; i++) {
      for (int j = 0; j < m_Ny; j++) {
        double nom_0 = m_Cw * I.real() * sigma(i, j)[0] * Hhat(i, j)[0];
        double nom_1 = m_Cw * I.imag() * sigma(i, j)[1] * Hhat(i, j)[1];

        double denom_0 = (1.0 - I.real() * m(i, j)[0] * m_Hw) * (1.0 + I.real() * sigma(i, j)[0] * m_tau_f) *
          (1.0 + I.real() * sigma(i, j)[0] * m_tau_c);
        double denom_1 = (1.0 - I.imag() * m(i, j)[1] * m_Hw) * (1.0 + I.imag() * sigma(i, j)[1] * m_tau_f) *
          (1.0 + I.imag() * sigma(i, j)[1] * m_tau_c);

        Phat(i, j)[0] = nom_0 / denom_0;
        Phat(i, j)[1] = nom_1 / denom_1;
      }
    }
  }

  // Save Phat in m_fftw_output.
  copy_fftw_array(m_Phat, m_fftw_output, m_Nx, m_Ny);
  fftw_execute(m_dft_inverse);

  // get m_fftw_output and put it into m_p
  get_real_part(m_fftw_output, 1.0 / (m_Nx * m_Ny), m_Mx, m_My, m_Nx, m_Ny, 0, 0,
                m_p);

  petsc::VecArray2D p(m_p, m_Mx, m_My);
  for (int i = 0; i < m_Mx; i++) {
    for (int j = 0; j < m_My; j++) {
      p(i, j) += m_background_precip_pre;
      if (m_truncate) {
        p(i, j) = std::min(p(i, j), 0.0);
      }
      p(i, j) *= m_precip_scale_factor;
      p(i, j) += m_background_precip_post;
    }
  }
}

} // end of namespace atmosphere
} // end of namespace pism
