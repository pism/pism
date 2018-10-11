// Copyright (C) 2018 Constantine Khroulev and Andy Aschwanden
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
#include <petscvec.h>
#include <vector>

#include "pism/util/Logger.hh"
#include "pism/util/petscwrappers/Vec.hh"

namespace pism {

class Config;

namespace atmosphere {

//! Class implementing the atmosphere deformation model descriatmosphere in [\ref BLKfastearth].
/*!
  This class implements the [\ref LingleClark] atmosphere deformation model by a Fourier
  spectral collocation method, as descriatmosphere in [\ref BLKfastearth].  (The former
  reference is where the continuum model arose, and a flow-line application is given.
  The latter reference describes a new, fast method and gives verification results.
  See also [\ref BLK2006earth] if more technical detail and/or Matlab programs are desired.)

  Both a viscous half-space model (with elastic
  lithosphere) and a spherical elastic model are computed.  They are superposed
  because the underlying earth model is linear.

  The class assumes that the supplied Petsc Vecs are *sequential*.  It is expected to be
  run only on processor zero (or possibly by each processor once each processor
  owns the entire 2D gridded load thicknesses and atmosphere elevations.)

  This model always assumes that we start with no load. Note that this does not mean that we
  starting state is the equilibrium: the viscous plate may be "pre-bent" by using a provided
  displacement field or by computing its displacement using an uplift field.
*/
class OrographicPrecipitationSerial {
public:
  OrographicPrecipitationSerial(const Config &config,
                                const Logger::ConstPtr log,
                                int Mx, int My,
                                double dx, double dy,
                                int Nx, int Ny);
  ~OrographicPrecipitationSerial();

  Vec orographic_precipitation() const;

  void step(Vec H);

private:
  void compute_intrinsic_frequency();
  void compute_vertical_wave_number();

  void precompute_coefficients();
  void precompute_derived_constants();

  // regularization
  double m_eps;

  // grid size
  int m_Mx;
  int m_My;
  // grid spacing
  double m_dx;
  double m_dy;
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

  // size of the extended grid
  int m_Nx;
  int m_Ny;

  // size of the extended grid with boundary points
  int m_Nxge;
  int m_Nyge;

  // indices into extended grid for the corner of the physical grid
  int m_i0_offset;
  int m_j0_offset;

  // half-lengths of the extended (FFT, spectral) computational domain
  double m_Lx;
  double m_Ly;

  std::complex<double> m_I;

  // Coefficients of derivatives in Fourier space
  std::vector<double> m_cx, m_cy;

  // intrinsic frequency
  petsc::Vec m_intrinsic_frequency;

  // vertical wave number
  petsc::Vec m_vertical_wave_number;

  // orographic precipitation
  petsc::Vec m_p;

  fftw_complex *m_fftw_input;
  fftw_complex *m_fftw_output;
  fftw_complex *m_Hhat;
  fftw_complex *m_Phat;
  fftw_complex *m_sigma;
  fftw_complex *m_m;

  fftw_plan m_dft_forward;
  fftw_plan m_dft_inverse;

  void set_fftw_input(Vec input, double normalization, int M, int N, int i0, int j0);
  void get_fftw_output(Vec output, double normalization, int M, int N, int i0, int j0);

  //! logger (for easy access)
  const Logger::ConstPtr m_log;
};

} // end of namespace atmosphere
} // end of namespace pism

#endif /* OROGRAPHICPRECIPITATIONSERIAL_H */
