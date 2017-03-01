// Copyright (C) 2004-2009, 2011, 2013, 2014, 2015, 2016, 2017 Ed Bueler and Constantine Khroulev
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
#include "deformation.hh"

#include "base/util/pism_const.hh"
#include "base/util/PISMConfigInterface.hh"
#include "base/util/error_handling.hh"
#include "base/util/petscwrappers/Vec.hh"

namespace pism {
namespace bed {

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
 * @param[in] include_elastic include elastic deformation component
 * @param[in] Mx grid size in the X direction
 * @param[in] My grid size in the Y direction
 * @param[in] dx grid spacing in the X direction
 * @param[in] dy grid spacing in the Y direction
 * @param[in] Z integer factor (the spectral grid is Z times bigger than the PISM grid)
 */
BedDeformLC::BedDeformLC(const Config &config,
                         bool include_elastic,
                         int Mx, int My,
                         double dx, double dy,
                         int Z) {

  // set parameters
  m_include_elastic = include_elastic;

  m_Mx     = Mx;
  m_My     = My;
  m_dx     = dx;
  m_dy     = dy;
  m_Z      = Z;
  m_load_density = config.get_double("constants.ice.density");
  m_mantle_density    = config.get_double("bed_deformation.mantle_density");
  m_eta    = config.get_double("bed_deformation.mantle_viscosity");
  m_D      = config.get_double("bed_deformation.lithosphere_flexural_rigidity");

  m_standard_gravity = config.get_double("constants.standard_gravity");

  // derive more parameters
  m_Nx        = m_Z*(m_Mx - 1) + 1;
  m_Ny        = m_Z*(m_My - 1) + 1;
  m_Lx        = 0.5 * (m_Nx - 1.0) * m_dx;
  m_Ly        = 0.5 * (m_Ny - 1.0) * m_dy;
  m_Nxge      = m_Nx + 1;
  m_Nyge      = m_Ny + 1;
  m_i0_offset = (m_Z - 1)*(m_Mx - 1) / 2;
  m_j0_offset = (m_Z - 1)*(m_My - 1) / 2;

  // memory allocation
  PetscErrorCode ierr = 0;

  // total plate displacement
  ierr = VecCreateSeq(PETSC_COMM_SELF, m_Mx * m_My, m_db.rawptr());;
  PISM_CHK(ierr, "VecCreateSeq");

  // ice thickness difference (the load)
  ierr = VecCreateSeq(PETSC_COMM_SELF, m_Mx * m_My, m_load_thickness.rawptr());;
  PISM_CHK(ierr, "VecCreateSeq");

  // elastic plate displacement
  ierr = VecCreateSeq(PETSC_COMM_SELF, m_Mx * m_My, m_db_elastic.rawptr());;
  PISM_CHK(ierr, "VecCreateSeq");

  // viscous plate displacement
  ierr = VecCreateSeq(PETSC_COMM_SELF, m_Nx * m_Ny, m_U.rawptr());
  PISM_CHK(ierr, "VecCreateSeq");

  // initial viscous plate displacement
  ierr = VecDuplicate(m_U, m_U_start.rawptr());
  PISM_CHK(ierr, "VecDuplicate");

  // elastic load response matrix
  ierr = VecCreateSeq(PETSC_COMM_SELF, m_Nxge * m_Nyge, m_load_response_matrix.rawptr());
  PISM_CHK(ierr, "VecCreateSeq");

  // setup fftw stuff: FFTW builds "plans" based on observed performance
  m_fftw_input  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * m_Nx * m_Ny);
  m_fftw_output = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * m_Nx * m_Ny);
  m_loadhat     = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * m_Nx * m_Ny);

  // fftw manipulates the data in setting up a plan, so fill with non-constant junk
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

  precompute_coefficients();
}

BedDeformLC::~BedDeformLC() {
  fftw_destroy_plan(m_dft_forward);
  fftw_destroy_plan(m_dft_inverse);
  fftw_free(m_fftw_input);
  fftw_free(m_fftw_output);
  fftw_free(m_loadhat);
}

Vec BedDeformLC::plate_displacement() const {
  return m_db;
}

/**
 * Pre-compute coefficients used by the model.
 */
void BedDeformLC::precompute_coefficients() {
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
 * @param[in] ice_thickness ice thickness, meters
 * @param[in] bed_uplift bed uplift, m/second
 *
 * Here `ice_thickness` is used to compute the load @f$ \sigma_{zz} @f$ and `bed_uplift` is
 * @f$ \diff{u}{t} @f$ itself.
 *
 * Sets m_U and m_U_start.
 */
void BedDeformLC::uplift_problem(Vec ice_thickness, Vec bed_uplift) {
  PetscErrorCode ierr = 0;

  // Compute fft2(-ice_rho * g * ice_thickness)
  {
    clear_fftw_input(m_fftw_input, m_Nx, m_Ny);
    set_fftw_input(ice_thickness, - m_load_density * m_standard_gravity, m_Mx, m_My, m_i0_offset, m_j0_offset);
    fftw_execute(m_dft_forward);
    // Save fft2(-ice_rho * g * ice_thickness) in loadhat.
    copy_fftw_output(m_fftw_output, m_loadhat, m_Nx, m_Ny);
  }

  // fft2(uplift)
  {
    clear_fftw_input(m_fftw_input, m_Nx, m_Ny);
    set_fftw_input(bed_uplift, 1.0, m_Mx, m_My, m_i0_offset, m_j0_offset);
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
  get_fftw_output(m_U_start, 1.0 / (m_Nx * m_Ny), m_Nx, m_Ny, 0, 0);

  tweak(m_U_start, m_Nx, m_Ny, 0.0);

  ierr = VecCopy(m_U_start, m_U); PISM_CHK(ierr, "VecCopy");

  // put m_U into "plate displacement" m_db.
  {
    // PISM grid
    petsc::VecArray2D db(m_db, m_Mx, m_My);
    // extended grid
    petsc::VecArray2D u(m_U, m_Nx, m_Ny, m_i0_offset, m_j0_offset);

    for (int i = 0; i < m_Mx; i++) {
      for (int j = 0; j < m_My; j++) {
        db(i, j) = u(i, j);
      }
    }
  }
}

/*!
 * to initialize we solve:
 *   rho_r g U + D grad^4 U = 0 - 2 eta |grad| uplift
 * where U=U_start; yes it really should have "0" on right side because
 * load for future times will be "-rho g (H-H_start)", which is zero if no geometry
 * change.
 *
 * Compare equation (16) in Bueler, Lingle, Brown (2007) "Fast computation of
 * a viscoelastic deformable Earth model for ice sheet simulations", Ann.
 * Glaciol. 46, 97--105
 * [NOTE PROBABLE SIGN ERROR in eqn (16)?:  load "rho g H" should be "- rho g H"]
 *
 */
void BedDeformLC::init(Vec uplift) {
  // zero load
  PetscErrorCode ierr = VecSet(m_load_thickness, 0.0); PISM_CHK(ierr, "VecSet");

  uplift_problem(m_load_thickness, uplift);
}

void BedDeformLC::step(double dt_seconds, Vec H_start, Vec H) {
  // solves:
  //     (2 eta |grad| U^{n+1}) + (dt/2) * (rho_r g U^{n+1} + D grad^4 U^{n+1})
  //   = (2 eta |grad| U^n) - (dt/2) * (rho_r g U^n + D grad^4 U^n) - dt * rho g H_start
  // where U=plate; see equation (7) in
  // Bueler, Lingle, Brown (2007) "Fast computation of a viscoelastic
  // deformable Earth model for ice sheet simulations", Ann. Glaciol. 46, 97--105

  // note ice thicknesses and bed elevations only on physical ("thin") grid
  //   while spectral/FFT quantities are on fat computational grid

  // Compute Hdiff
  // Hdiff = H - H_start
  PetscErrorCode ierr = VecWAXPY(m_load_thickness, -1, H_start, H); PISM_CHK(ierr, "VecWAXPY");

  // Compute fft2(-ice_rho * g * dH * dt), where H = H - H_start.
  {
    clear_fftw_input(m_fftw_input, m_Nx, m_Ny);
    set_fftw_input(m_load_thickness, - m_load_density * m_standard_gravity * dt_seconds,
                   m_Mx, m_My, m_i0_offset, m_j0_offset);
    fftw_execute(m_dft_forward);

    // Save fft2(-ice_rho * g * dH * dt) in loadhat.
    copy_fftw_output(m_fftw_output, m_loadhat, m_Nx, m_Ny);
  }

  // Compute fft2(u).
  // no need to clear fftw_input: all values are overwritten
  {
    set_fftw_input(m_U, 1.0, m_Nx, m_Ny, 0, 0);
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
  get_fftw_output(m_U, 1.0 / (m_Nx * m_Ny), m_Nx, m_Ny, 0, 0);

  // Now tweak. (See the "correction" in section 5 of BuelerLingleBrown.)
  //
  // Here 1e16 approximates t = \infty.
  tweak(m_U, m_Nx, m_Ny, 1e16);

  // now compute elastic response if desired; bed = ue at end of this block
  if (m_include_elastic) {
    // MATLAB:     ue=rhoi*conv2(H-H_start, II, 'same')
    conv2_same(m_load_thickness, m_Mx, m_My, m_load_response_matrix, m_Nxge, m_Nyge, m_db_elastic);

    ierr = VecScale(m_db_elastic, m_load_density); PISM_CHK(ierr, "VecScale");
  } else {
    ierr = VecSet(m_db_elastic, 0.0); PISM_CHK(ierr, "VecSet");
  }

  // compute total plate displacement relative to the initial bed elevation
  {
    // PISM grid
    petsc::VecArray2D
      db(m_db, m_Mx, m_My),
      db_elastic(m_db_elastic, m_Mx, m_My);
    // extended grid
    petsc::VecArray2D
      u(m_U, m_Nx, m_Ny, m_i0_offset, m_j0_offset),
      u_start(m_U_start, m_Nx, m_Ny, m_i0_offset, m_j0_offset);

    for (int i = 0; i < m_Mx; i++) {
      for (int j = 0; j < m_My; j++) {
        db(i, j) = (u(i, j) - u_start(i, j)) + db_elastic(i, j);
      }
    }
  }
}

void BedDeformLC::tweak(Vec U, int Nx, int Ny, double time) {
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
    ierr = VecSum(m_load_thickness, &H_sum); PISM_CHK(ierr, "VecSum");

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
void BedDeformLC::set_fftw_input(Vec vec_input, double normalization, int Mx, int My, int i0, int j0) {
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
void BedDeformLC::get_fftw_output(Vec output, double normalization, int Mx, int My, int i0, int j0) {
  petsc::VecArray2D out(output, Mx, My);
  VecAccessor2D<fftw_complex> fftw_out(m_fftw_output, m_Nx, m_Ny, i0, j0);
  for (int j = 0; j < My; ++j) {
    for (int i = 0; i < Mx; ++i) {
      out(i, j) = fftw_out(i, j)[0] * normalization;
    }
  }
}

} // end of namespace bed
} // end of namespace pism
