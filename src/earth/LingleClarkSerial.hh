// Copyright (C) 2007--2009, 2011, 2012, 2013, 2014, 2015, 2017, 2018 Ed Bueler and Constantine Khroulev
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

#ifndef LINGLECLARKSERIAL_H
#define LINGLECLARKSERIAL_H

#include <vector>

#include <petscvec.h>
#include <fftw3.h>
#include <vector>

#include "pism/util/petscwrappers/Vec.hh"

namespace pism {

class Config;

namespace bed {

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
  owns the entire 2D gridded load thicknesses and bed elevations.)

  This model always assumes that we start with no load. Note that this does not mean that we
  starting state is the equilibrium: the viscous plate may be "pre-bent" by using a provided
  displacement field or by computing its displacement using an uplift field.
*/
class LingleClarkSerial {
public:
  LingleClarkSerial(const Config &config,
                    bool include_elastic,
                    int Mx, int My,
                    double dx, double dy,
                    int Nx, int Ny);
  ~LingleClarkSerial();

  void init(Vec thickness, Vec viscous_displacement);

  void bootstrap(Vec thickness, Vec uplift);

  void step(double dt_seconds, Vec H);

  Vec total_displacement() const;

  Vec viscous_displacement() const;
private:
  void compute_elastic_response(Vec H, Vec dE);

  void uplift_problem(Vec load_thickness, Vec bed_uplift, Vec output);

  void precompute_coefficients();

  void update_displacement(Vec V, Vec dE, Vec dU);

  bool m_include_elastic;
  // grid size
  int m_Mx;
  int m_My;
  // grid spacing
  double m_dx;
  double m_dy;
  //! load density (for computing load from its thickness)
  double m_load_density;
  //! mantle density
  double m_mantle_density;
  //! mantle viscosity
  double m_eta;
  //! lithosphere flexural rigidity
  double m_D;

  // acceleration due to gravity
  double m_standard_gravity;

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

  // Coefficients of derivatives in Fourier space
  std::vector<double> m_cx, m_cy;

  // viscous displacement on the extended grid
  petsc::Vec m_Uv;

  // load response matrix (elastic); sequential and fat *with* boundary
  petsc::Vec m_load_response_matrix;
  // elastic plate displacement
  petsc::Vec m_Ue;

  // total (viscous and elastic) plate displacement
  petsc::Vec m_U;

  fftw_complex *m_fftw_input;
  fftw_complex *m_fftw_output;
  fftw_complex *m_loadhat;

  fftw_plan m_dft_forward;
  fftw_plan m_dft_inverse;

  void tweak(Vec load_thickness, Vec U, int Nx, int Ny, double time);

  void set_fftw_input(Vec input, double normalization, int M, int N, int i0, int j0);
  void get_fftw_output(Vec output, double normalization, int M, int N, int i0, int j0);
};

} // end of namespace bed
} // end of namespace pism

#endif /* LINGLECLARKSERIAL_H */
