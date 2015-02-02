// Copyright (C) 2007--2009, 2011, 2012, 2013, 2014, 2015 Ed Bueler
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

#ifndef __deformation_hh
#define __deformation_hh

#include <petscvec.h>
#include <fftw3.h>

#include "Vec.hh"

namespace pism {

class Config;

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
  owns the entire 2D gridded ice thicknesses and bed elevations.)

  A test program for this class is pism/src/verif/tryLCbd.cc.
*/
class BedDeformLC {
public:
  BedDeformLC();
  ~BedDeformLC();
  void settings(const Config &config,
                bool myinclude_elastic,
                int myMx, int myMy, double mydx, double mydy,
                int myZ,
                Vec myHstart, Vec mybedstart, Vec myuplift,  // initial state
                Vec myH,     // generally gets changed by calling program
                // before each call to step
                Vec mybed);  // mybed gets modified by step()
  void alloc();
  void init();
  void uplift_init();
  void step(const double dtyear, const double yearFromStart);

protected:
  bool        m_include_elastic;
  int         m_Mx, m_My;
  double   m_dx, m_dy;
  //! Factor by which fat FFT domain is larger than region of physical
  //! interest.
  int m_Z;
  double m_icerho,           // ice density (for computing load from volume)
    m_rho,                        // earth density
    m_eta,                        // mantle viscosity
    m_D;                          // lithosphere flexural rigidity

private:
  double m_standard_gravity;
  bool m_settingsDone, m_allocDone;
  int m_Nx, m_Ny,         // fat sizes
    m_Nxge, m_Nyge;     // fat with boundary sizes
  int      m_i0_plate,  m_j0_plate; // indices into fat array for corner of thin
  double   m_Lx, m_Ly;         // half-lengths of the physical domain
  double   m_Lx_fat, m_Ly_fat; // half-lengths of the FFT (spectral) computational domain
  std::vector<double>  m_cx, m_cy;        // coeffs of derivatives in Fourier space

  // point to storage owned elsewhere
  Vec m_H, m_bed, m_H_start, m_bed_start, m_uplift;

  petsc::Vec m_Hdiff, m_dbedElastic, // sequential; working space
    m_U, m_U_start,     // sequential and fat
    m_vleft, m_vright,  // coefficients; sequential and fat
    m_lrmE;           // load response matrix (elastic); sequential and fat *with* boundary

  fftw_complex *m_fftw_input, *m_fftw_output, *m_loadhat; // 2D sequential
  fftw_plan m_dft_forward, m_dft_inverse;

  void tweak(double seconds_from_start);

  void clear_fftw_input();
  void copy_fftw_output(fftw_complex *buffer);
  void set_fftw_input(Vec input, double normalization, int M, int N, int i0, int j0);
  void get_fftw_output(Vec output, double normalization, int M, int N, int i0, int j0);
};

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

} // end of namespace pism

#endif  /* __deformation_hh */

