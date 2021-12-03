// Copyright (C) 2018, 2020 Constantine Khroulev and Andy Aschwanden
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

#ifndef OROGRAPHICPRECIPITATIONSERIAL_H
#define OROGRAPHICPRECIPITATIONSERIAL_H

#include <vector>
#include <fftw3.h>

#include "pism/util/petscwrappers/Vec.hh"

namespace pism {

class Config;

namespace atmosphere {

//! Class implementing the linear model of orographic precipitation [@ref
//! SmithBarstad2004], [@ref SmithBarstadBonneau2005].
class OrographicPrecipitationSerial {
public:
  OrographicPrecipitationSerial(const Config &config,
                                int Mx, int My,
                                double dx, double dy,
                                int Nx, int Ny);
  ~OrographicPrecipitationSerial();

  Vec precipitation() const;

  void update(petsc::Vec &surface_elevation);

private:
  // regularization
  double m_eps;

  // grid size
  int m_Mx;
  int m_My;

  //! truncate
  bool m_truncate;
  //! precipitation scale factor
  double m_precip_scale_factor;
  //! background precipitation
  double m_background_precip_pre, m_background_precip_post;
  //! cloud conversion time
  double m_tau_c;
  //! cloud fallout time
  double m_tau_f;
  //! water vapor scale height
  double m_Hw;
  //! moist stability frequency
  double m_Nm;
  //! wind direction
  double m_wind_direction;
  //! wind speed
  double m_wind_speed;
  //! moist adiabatic lapse rate
  double m_Theta_m;
  //! moist lapse rate
  double m_gamma;
  //! reference density
  double m_rho_Sref;
  //! Coriolis force
  double m_f;
  //! uplift sensitivity factor
  double m_Cw;
  //! latitude for Coriolis force
  double m_latitude;
  //! horizontal wind component
  double m_u;
  //! vertical wind component
  double m_v;

  // extended grid size
  int m_Nx;
  int m_Ny;

  // indices into extended grid for the corner of the physical grid
  int m_i0_offset;
  int m_j0_offset;

  std::vector<double> m_kx, m_ky;

  // resulting orographic precipitation
  petsc::Vec m_precipitation;

  fftw_complex *m_fftw_input;
  fftw_complex *m_fftw_output;

  // FFT(Gaussian) used to smooth surface elevation
  fftw_complex *m_G_hat;

  fftw_plan m_dft_forward;
  fftw_plan m_dft_inverse;
};

} // end of namespace atmosphere
} // end of namespace pism

#endif /* OROGRAPHICPRECIPITATIONSERIAL_H */
