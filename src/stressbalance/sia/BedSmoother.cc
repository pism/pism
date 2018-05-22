// Copyright (C) 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018 Ed Bueler and Constantine Khroulev
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

#include "BedSmoother.hh"
#include "pism/util/Mask.hh"
#include "pism/util/IceGrid.hh"
#include "pism/util/petscwrappers/Vec.hh"
#include "pism/util/IceModelVec2CellType.hh"

#include "pism/util/error_handling.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/util/Logger.hh"

namespace pism {
namespace stressbalance {

BedSmoother::BedSmoother(IceGrid::ConstPtr g, int MAX_GHOSTS)
    : m_grid(g), m_config(g->ctx()->config()) {

  const Logger &log = *m_grid->ctx()->log();

  {
    // allocate Vecs that live on all procs; all have to be as "wide" as any of
    //   their prospective uses
    m_topgsmooth.create(m_grid, "topgsmooth", WITH_GHOSTS, MAX_GHOSTS);
    m_topgsmooth.set_attrs("bed_smoother_tool",
                         "smoothed bed elevation, in bed roughness parameterization",
                         "m", "");
    m_maxtl.create(m_grid, "maxtl", WITH_GHOSTS, MAX_GHOSTS);
    m_maxtl.set_attrs("bed_smoother_tool",
                    "maximum elevation in local topography patch, in bed roughness parameterization",
                    "m", "");
    m_C2.create(m_grid, "C2bedsmooth", WITH_GHOSTS, MAX_GHOSTS);
    m_C2.set_attrs("bed_smoother_tool",
                 "polynomial coeff of H^-2, in bed roughness parameterization",
                 "m2", "");
    m_C3.create(m_grid, "C3bedsmooth", WITH_GHOSTS, MAX_GHOSTS);
    m_C3.set_attrs("bed_smoother_tool",
                 "polynomial coeff of H^-3, in bed roughness parameterization",
                 "m3", "");
    m_C4.create(m_grid, "C4bedsmooth", WITH_GHOSTS, MAX_GHOSTS);
    m_C4.set_attrs("bed_smoother_tool",
                 "polynomial coeff of H^-4, in bed roughness parameterization",
                 "m4", "");

    // allocate Vecs that live on processor 0:
    m_topgp0 = m_topgsmooth.allocate_proc0_copy();
    m_topgsmoothp0 = m_topgsmooth.allocate_proc0_copy();
    m_maxtlp0 = m_maxtl.allocate_proc0_copy();
    m_C2p0 = m_C2.allocate_proc0_copy();
    m_C3p0 = m_C3.allocate_proc0_copy();
    m_C4p0 = m_C4.allocate_proc0_copy();
  }

  m_Glen_exponent = m_config->get_double("stress_balance.sia.Glen_exponent"); // choice is SIA; see #285
  m_smoothing_range = m_config->get_double("stress_balance.sia.bed_smoother.range");

  if (m_smoothing_range > 0.0) {
    log.message(2,
                "* Initializing bed smoother object with %.3f km half-width ...\n",
                units::convert(m_grid->ctx()->unit_system(), m_smoothing_range, "m", "km"));
  }

  // Make sure that Nx and Ny are initialized. In most cases SIAFD::update() will call
  // preprocess_bed() and set appropriate values, but in a zero-length (-y 0) run IceModel does not
  // call SIAFD::update()... We may need to re-structure this class so that everything is
  // initialized right after construction and users don't have to call preprocess_bed() manually.
  m_Nx = -1;
  m_Ny = -1;
}


BedSmoother::~BedSmoother() {
  // empty
}

/*!
Input lambda gives physical half-width (in m) of square over which to do the
average.  Only square smoothing domains are allowed with this call, which is the
default case.
 */
void BedSmoother::preprocess_bed(const IceModelVec2S &topg) {

  if (m_smoothing_range <= 0.0) {
    // smoothing completely inactive.  we transfer the original bed topg,
    //   including ghosts, to public member topgsmooth ...
    topg.update_ghosts(m_topgsmooth);
    // and we tell theta() to return theta=1
    m_Nx = -1;
    m_Ny = -1;
    return;
  }

  // determine Nx, Ny, which are always at least one if m_smoothing_range > 0
  m_Nx = static_cast<int>(ceil(m_smoothing_range / m_grid->dx()));
  m_Ny = static_cast<int>(ceil(m_smoothing_range / m_grid->dy()));
  if (m_Nx < 1) {
    m_Nx = 1;
  }
  if (m_Ny < 1) {
    m_Ny = 1;
  }

  preprocess_bed(topg, m_Nx, m_Ny);
}

const IceModelVec2S& BedSmoother::smoothed_bed() const {
  return m_topgsmooth;
}

/*!
Inputs Nx,Ny gives half-width in number of grid points, over which to do the
average.
 */
void BedSmoother::preprocess_bed(const IceModelVec2S &topg,
                                 unsigned int Nx, unsigned int Ny) {

  if ((Nx >= m_grid->Mx()) || (Ny >= m_grid->My())) {
    throw RuntimeError(PISM_ERROR_LOCATION, "input Nx, Ny in bed smoother is too large because\n"
                       "domain of smoothing exceeds IceGrid domain");
  }
  m_Nx = Nx;
  m_Ny = Ny;

  topg.put_on_proc0(*m_topgp0);
  smooth_the_bed_on_proc0();
  // next call *does indeed* fill ghosts in topgsmooth
  m_topgsmooth.get_from_proc0(*m_topgsmoothp0);

  compute_coefficients_on_proc0();
  // following calls *do* fill the ghosts
  m_maxtl.get_from_proc0(*m_maxtlp0);
  m_C2.get_from_proc0(*m_C2p0);
  m_C3.get_from_proc0(*m_C3p0);
  m_C4.get_from_proc0(*m_C4p0);
}


//! Computes the smoothed bed by a simple average over a rectangle of grid points.
void BedSmoother::smooth_the_bed_on_proc0() {

  ParallelSection rank0(m_grid->com);
  try {
    if (m_grid->rank() == 0) {
      const int Mx = (int)m_grid->Mx();
      const int My = (int)m_grid->My();

      petsc::VecArray2D
        b0(*m_topgp0,       Mx, My),
        bs(*m_topgsmoothp0, Mx, My);

      for (int j=0; j < My; j++) {
        for (int i=0; i < Mx; i++) {
          // average only over those points which are in the grid; do
          // not wrap periodically
          double sum = 0.0, count = 0.0;
          for (int r = -m_Nx; r <= m_Nx; r++) {
            for (int s = -m_Ny; s <= m_Ny; s++) {
              if ((i+r >= 0) and (i+r < Mx) and (j+s >= 0) and (j+s < My)) {
                sum   += b0(i+r, j+s);
                count += 1.0;
              }
            }
          }
          // unprotected division by count but r=0,s=0 case guarantees count>=1
          bs(i, j) = sum / count;
        }
      }
    }
  } catch (...) {
    rank0.failed();
  }
  rank0.check();
}


void BedSmoother::compute_coefficients_on_proc0() {

  const unsigned int Mx = m_grid->Mx(), My = m_grid->My();

  ParallelSection rank0(m_grid->com);
  try {
    if (m_grid->rank() == 0) {
      petsc::VecArray2D
        b0(*m_topgp0,       Mx, My),
        bs(*m_topgsmoothp0, Mx, My),
        mt(*m_maxtlp0,      Mx, My),
        c2(*m_C2p0,         Mx, My),
        c3(*m_C3p0,         Mx, My),
        c4(*m_C4p0,         Mx, My);

      for (int j=0; j < (int)My; j++) {
        for (int i=0; i < (int)Mx; i++) {
          // average only over those points which are in the grid
          // do not wrap periodically
          double
            topgs     = bs(i, j),
            maxtltemp = 0.0,
            sum2      = 0.0,
            sum3      = 0.0,
            sum4      = 0.0,
            count     = 0.0;

          for (int r = -m_Nx; r <= m_Nx; r++) {
            for (int s = -m_Ny; s <= m_Ny; s++) {
              if ((i+r >= 0) && (i+r < (int)Mx) && (j+s >= 0) && (j+s < (int)My)) {
                // tl is elevation of local topography at a pt in patch
                const double tl  = b0(i+r, j+s) - topgs;
                maxtltemp = std::max(maxtltemp, tl);
                // accumulate 2nd, 3rd, and 4th powers with only 3 multiplications
                const double tl2 = tl * tl;
                sum2 += tl2;
                sum3 += tl2 * tl;
                sum4 += tl2 * tl2;
                count += 1.0;
              }
            }
          }
          mt(i, j) = maxtltemp;

          // unprotected division by count but r=0,s=0 case guarantees count>=1
          c2(i, j) = sum2 / count;
          c3(i, j) = sum3 / count;
          c4(i, j) = sum4 / count;
        }
      }

      // scale the coeffs in Taylor series
      const double
        n = m_Glen_exponent,
        k  = (n + 2) / n,
        s2 = k * (2 * n + 2) / (2 * n),
        s3 = s2 * (3 * n + 2) / (3 * n),
        s4 = s3 * (4 * n + 2) / (4 * n);

      PetscErrorCode ierr;
      ierr = VecScale(*m_C2p0,s2);
      PISM_CHK(ierr, "VecScale");

      ierr = VecScale(*m_C3p0,s3);
      PISM_CHK(ierr, "VecScale");

      ierr = VecScale(*m_C4p0,s4);
      PISM_CHK(ierr, "VecScale");
    }
  } catch (...) {
    rank0.failed();
  }
  rank0.check();
}


//! Computes a smoothed thickness map.
/*!
The result `thksmooth` is the difference between the given upper surface
elevation (`usurf`) and the stored smoothed bed topography (`topgsmooth`),
except where the given original thickness (`thk`) is zero.  In places where
the original thickness is zero, the result `thksmooth` is also set to zero.

Ghosted values are updated directly and no communication occurs.  In fact,
we \e assume `usurf`, `thk`, and `thksmooth` all have stencil width at least
equal to GHOSTS.  We \e check whether `topgsmooth`, which has stencil width
maxGHOSTS, has at least GHOSTS stencil width, and throw an error if not.

Call preprocess_bed() first.
 */
void BedSmoother::smoothed_thk(const IceModelVec2S &usurf,
                               const IceModelVec2S &thk,
                               const IceModelVec2CellType &mask,
                               IceModelVec2S &result) const {

  IceModelVec::AccessList list{&mask, &m_maxtl, &result, &thk, &m_topgsmooth, &usurf};

  unsigned int GHOSTS = result.stencil_width();
  assert(mask.stencil_width()         >= GHOSTS);
  assert(m_maxtl.stencil_width()      >= GHOSTS);
  assert(thk.stencil_width()          >= GHOSTS);
  assert(m_topgsmooth.stencil_width() >= GHOSTS);
  assert(usurf.stencil_width()        >= GHOSTS);

  ParallelSection loop(m_grid->com);
  try {
    for (PointsWithGhosts p(*m_grid, GHOSTS); p; p.next()) {
      const int i = p.i(), j = p.j();

      if (thk(i, j) < 0.0) {
        throw RuntimeError::formatted(PISM_ERROR_LOCATION, "BedSmoother detects negative original thickness\n"
                                      "at location (i, j) = (%d, %d) ... ending", i, j);
      } else if (thk(i, j) == 0.0) {
        result(i, j) = 0.0;
      } else if (m_maxtl(i, j) >= thk(i, j)) {
        result(i, j) = thk(i, j);
      } else {
        if (mask.grounded(i, j)) {
          // if grounded, compute smoothed thickness as the difference of ice
          // surface elevation and smoothed bed elevation
          const double thks_try = usurf(i, j) - m_topgsmooth(i, j);
          result(i, j) = (thks_try > 0.0) ? thks_try : 0.0;
        } else {
          // if floating, use original thickness (note: surface elevation was
          // computed using this thickness and the sea level elevation)
          result(i, j) = thk(i, j);
        }
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();
}


/*!
Implements the strategy for computing \f$\theta(h,x,y)\f$ from previously-
stored coefficients, described on [Bed roughness parameterization](@ref bedrough) page and in [\ref
Schoofbasaltopg2003].

Specifically, \f$\theta = \omega^{-n}\f$ where \f$\omega\f$ is a local average
of a rational function of surface elevation, approximated here by a Taylor polynomial:
  \f[ \omega = \fint \left(1 - \frac{\tilde b(x_1,x_2,\xi_1,\xi_2)}{H}
                           \right)^{-(n+2)/n}\,d\xi_1\,d\xi_2
             \approx 1 + C_2 H^{-2} + C_3 H^{-3} + C_4 H^{-4} \f]
where \f$h =\f$ usurf, \f$H = h -\f$ topgsmooth and \f$\tilde b\f$ is the local
bed topography, a function with mean zero.  The coefficients \f$C_2,C_3,C_4\f$,
which depend on \f$x,y\f$, are precomputed by `preprocess_bed()`.

Ghosted values are updated directly and no communication occurs.  In fact,
we \e assume `usurf` and `theta` have stencil width at least
equal to GHOSTS.  We \e check whether `topgsmooth`, which has stencil width
maxGHOSTS, has at least GHOSTS stencil width, and throw an error if not.

Call preprocess_bed() first.
 */
void BedSmoother::theta(const IceModelVec2S &usurf, IceModelVec2S &result) const {

  if ((m_Nx < 0) || (m_Ny < 0)) {
    result.set(1.0);
    return;
  }

  IceModelVec::AccessList list{&m_C2, &m_C3, &m_C4, &m_maxtl, &result, &m_topgsmooth, &usurf};

  unsigned int GHOSTS = result.stencil_width();
  assert(m_C2.stencil_width()         >= GHOSTS);
  assert(m_C3.stencil_width()         >= GHOSTS);
  assert(m_C4.stencil_width()         >= GHOSTS);
  assert(m_maxtl.stencil_width()      >= GHOSTS);
  assert(m_topgsmooth.stencil_width() >= GHOSTS);
  assert(usurf.stencil_width()        >= GHOSTS);

  const double
    theta_min = m_config->get_double("stress_balance.sia.bed_smoother.theta_min"),
    theta_max = 1.0;

  ParallelSection loop(m_grid->com);
  try {
    for (PointsWithGhosts p(*m_grid, GHOSTS); p; p.next()) {
      const int i = p.i(), j = p.j();

      const double H = usurf(i, j) - m_topgsmooth(i, j);
      if (H > m_maxtl(i, j)) {
        // thickness exceeds maximum variation in patch of local topography,
        // so ice buries local topography; note maxtl >= 0 always
        const double Hinv = 1.0 / std::max(H, 1.0);
        double omega = 1.0 + Hinv*Hinv * (m_C2(i, j) + Hinv * (m_C3(i, j) + Hinv*m_C4(i, j)));
        if (omega <= 0) {  // this check *should not* be necessary: p4(s) > 0
          throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                        "omega is negative for i=%d, j=%d\n"
                                        "in BedSmoother.theta()", i, j);
        }

        if (omega < 0.001) {      // this check *should not* be necessary
          omega = 0.001;
        }

        result(i, j) = pow(omega, -m_Glen_exponent);
      } else {
        result(i, j) = 0.00;
      }

      result(i, j) = clip(result(i, j), theta_min, theta_max);
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();
}


} // end of namespace stressbalance
} // end of namespace pism
