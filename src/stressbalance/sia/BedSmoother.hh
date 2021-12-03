// Copyright (C) 2010, 2011, 2013, 2014, 2015, 2016, 2017, 2020, 2021 Ed Bueler and Constantine Khroulev
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

#ifndef __BedSmoother_hh
#define __BedSmoother_hh

#include <petsc.h>

#include "pism/util/iceModelVec.hh"
#include "pism/util/ConfigInterface.hh"

namespace pism {

class IceGrid;
class Config;
class IceModelVec2CellType;

namespace stressbalance {

//! PISM bed smoother, plus bed roughness parameterization, based on Schoof (2003).
/*!
  This class both smooths the bed and computes coefficients for an approximation
  to Schoof's \f$\theta\f$.  The factor \f$0\le \theta \le 1\f$ multiplies the diffusivity
  in the theory of the effect of bed roughness in the SIA by Christian Schoof
  (2003; *The effect of basal topography on ice sheet dynamics*) [\ref
  Schoofbasaltopg2003].

  For additional information on this class see page \ref bedrough.

  The user of this class hands BedSmoother an "original" topography, and it
  is preprocessed to fill the smoothed topography `topgsmooth`, and the
  coefficients in an approximation to \f$\theta\f$.  This is done by a call to 
  `preprocess_bed()`.  The call requires the half-width of the smoothing square
  (a distance in m), or the number of grid points in each direction in the
  smoothing rectangle, and the Glen exponent.

  The call to `preprocess_bed()` *must be repeated* any time the "original"
  topography changes, for instance at the start of an IceModel run, or at a bed
  deformation step in an IceModel run.

  BedSmoother then provides three major functionalities, all of which \e must
  \e follow the call to `preprocess_bed()`:
  -# User accesses public IceModelVec2S `topgsmooth`, the smoothed bed itself.
  -# User asks `get_smoothed_thk()` for gridded values of the consistent smoothed
  version of the ice thickness, which is the thickness corresponding to a given
  surface elevation and the pre-computed smoothed bed.
  -# User asks for gridded values of \f$\theta(h,x,y)\f$ using `get_theta()`.

  Here is a basic example of the creation and usage of a BedSmoother instance.
  Note that BedSmoother will update ghosted values in `topgsmooth`, and other
  internal fields, and update them in the return fields `thksmooth`, `theta`,
  if asked.  In IceModel::velocitySIAStaggered()
  \code
  BedSmoother smoother(grid, 2);
  const double n = 3.0,
  lambda = 50.0e3;
  ierr = smoother.preprocess_bed(topg, n, lambda); CHKERRQ(ierr);
  ierr = smoother.get_smoothed_thk(usurf, thk, 1, &thksmooth); CHKERRQ(ierr);
  ierr = smoother.get_theta(usurf, n, 1, &theta); CHKERRQ(ierr);
  \endcode
  See IceGrid documentation for initializing `grid`. Note we assume
  `topg`, `usurf`, `thk`, `thksmooth`, and `theta` are all created
  IceModelVec2S instances.
*/
class BedSmoother {
public:
  BedSmoother(IceGrid::ConstPtr g, int MAX_GHOSTS);
  virtual ~BedSmoother() = default;

  void preprocess_bed(const IceModelVec2S &topg);

  void smoothed_thk(const IceModelVec2S &usurf,
                    const IceModelVec2S &thk,
                    const IceModelVec2CellType &mask,
                    IceModelVec2S &thksmooth) const;

  void theta(const IceModelVec2S &usurf, IceModelVec2S &result) const;

  const IceModelVec2S& smoothed_bed() const;
protected:
  IceGrid::ConstPtr m_grid;
  const Config::ConstPtr m_config;

  //! smoothed bed elevation; set by calling preprocess_bed()
  IceModelVec2S m_topgsmooth;

  IceModelVec2S m_maxtl, m_C2, m_C3, m_C4;

  /* number of grid points to smooth over; e.g. i=-Nx,-Nx+1,...,-1,0,1,...,Nx-1,Nx; note
    Nx>=1 and Ny>=1 always, unless lambda<=0
   */
  int m_Nx, m_Ny;

  double m_Glen_exponent, m_smoothing_range;

  //! original bed elevation on processor 0
  std::shared_ptr<petsc::Vec> m_topgp0;
  //! smoothed bed elevation on processor 0
  std::shared_ptr<petsc::Vec> m_topgsmoothp0;
  //! maximum elevation at (i,j) of local topography (nearby patch)
  std::shared_ptr<petsc::Vec> m_maxtlp0, m_C2p0, m_C3p0, m_C4p0;

  void preprocess_bed(const IceModelVec2S &topg,
                      unsigned int Nx_in, unsigned int Ny_in);

  void smooth_the_bed_on_proc0();
  void compute_coefficients_on_proc0();
};

} // end of namespace stressbalance
} // end of namespace pism

#endif  // __BedSmoother_hh

