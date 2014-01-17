// Copyright (C) 2010, 2011, 2013, 2014 Ed Bueler and Constantine Khroulev
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

#ifndef __PISMBedSmoother_hh
#define __PISMBedSmoother_hh

#include <petsc.h>
#include "iceModelVec.hh"

class IceGrid;
class PISMConfig;

//! PISM bed smoother, plus bed roughness parameterization, based on Schoof (2003).
/*!
This class both smooths the bed and computes coefficients for an approximation
to Schoof's \f$\theta\f$.  The factor \f$0\le \theta \le 1\f$ multiplies the diffusivity
in the theory of the effect of bed roughness in the SIA by Christian Schoof
(2003; *The effect of basal topography on ice sheet dynamics*) [\ref
Schoofbasaltopg2003].

For additional information on this class see page \ref bedrough.

The user of this class hands PISMBedSmoother an "original" topography, and it
is preprocessed to fill the smoothed topography `topgsmooth`, and the
coefficients in an approximation to \f$\theta\f$.  This is done by a call to 
`preprocess_bed()`.  The call requires the half-width of the smoothing square
(a distance in m), or the number of grid points in each direction in the
smoothing rectangle, and the Glen exponent.

The call to `preprocess_bed()` <b>must be repeated</b> any time the "original"
topography changes, for instance at the start of an IceModel run, or at a bed
deformation step in an IceModel run.

PISMBedSmoother then provides three major functionalities, all of which \e must
\e follow the call to `preprocess_bed()`:
-# User accesses public IceModelVec2S `topgsmooth`, the smoothed bed itself.
-# User asks `get_smoothed_thk()` for gridded values of the consistent smoothed
   version of the ice thickness, which is the thickness corresponding to a given
   surface elevation and the pre-computed smoothed bed.
-# User asks for gridded values of \f$\theta(h,x,y)\f$ using `get_theta()`.

Here is a basic example of the creation and usage of a PISMBedSmoother instance.
Note that PISMBedSmoother will update ghosted values in `topgsmooth`, and other
internal fields, and update them in the return fields `thksmooth`, `theta`,
if asked.  In IceModel::velocitySIAStaggered()
\code
    PISMBedSmoother smoother(grid, config, 2);
    const PetscReal n = 3.0,
                    lambda = 50.0e3;
    ierr = smoother.preprocess_bed(topg, n, lambda); CHKERRQ(ierr);
    ierr = smoother.get_smoothed_thk(usurf, thk, 1, &thksmooth); CHKERRQ(ierr);
    ierr = smoother.get_theta(usurf, n, 1, &theta); CHKERRQ(ierr);
\endcode
See IceGrid and PISMConfig documentation for initializing `grid` and
`config`.  Note we assume `topg`, `usurf`, `thk`, `thksmooth`, and `theta`
are all created IceModelVec2S instances.
 */
class PISMBedSmoother {
public:
  PISMBedSmoother(IceGrid &g, const PISMConfig &conf, PetscInt MAX_GHOSTS);
  virtual ~PISMBedSmoother();

  virtual PetscErrorCode preprocess_bed(IceModelVec2S &topg);

  // FIXME: this method is used exactly once in bedrough_test.cc. Consider removing it.
  virtual PetscErrorCode get_smoothing_domain(PetscInt &Nx_out, PetscInt &Ny_out);

  virtual PetscErrorCode get_smoothed_thk(IceModelVec2S &usurf, IceModelVec2S &thk,
                                          IceModelVec2Int &mask, IceModelVec2S *thksmooth);
  virtual PetscErrorCode get_theta(IceModelVec2S &usurf, IceModelVec2S *theta);

  //! smoothed bed elevation; publicly-available; set by calling preprocess_bed()
  IceModelVec2S topgsmooth;

protected:
  IceGrid &grid;
  const PISMConfig &config;
  IceModelVec2S maxtl, C2, C3, C4;

  PetscInt Nx, Ny;  //!< number of grid points to smooth over; e.g.
                    //!i=-Nx,-Nx+1,...,-1,0,1,...,Nx-1,Nx; note Nx>=1 and Ny>=1
                    //!always, unless lambda<=0

  PetscReal m_Glen_exponent, m_smoothing_range;

  PetscErrorCode allocate(int MAX_GHOSTS);
  PetscErrorCode deallocate();

  Vec g2, g2natural;  //!< global Vecs used to transfer data to/from processor 0.
  VecScatter scatter; //!< VecScatter used to transfer data to/from processor 0.
  Vec topgp0,         //!< original bed elevation on processor 0
      topgsmoothp0,   //!< smoothed bed elevation on processor 0
      maxtlp0,        //!< maximum elevation at (i,j) of local topography (nearby patch)
      C2p0, C3p0, C4p0;

  virtual PetscErrorCode preprocess_bed(IceModelVec2S &topg,
                                        PetscInt Nx_in, PetscInt Ny_in);

  PetscErrorCode smooth_the_bed_on_proc0();
  PetscErrorCode compute_coefficients_on_proc0();
};

#endif	// __PISMBedSmoother_hh

