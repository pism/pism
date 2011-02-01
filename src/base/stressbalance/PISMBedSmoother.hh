// Copyright (C) 2010, 2011 Ed Bueler and Constantine Khroulev
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 2 of the License, or (at your option) any later
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
#include "NCVariable.hh"
#include "grid.hh"
#include "iceModelVec.hh"

//! PISM bed smoother, plus bed roughness parameterization, based on Schoof (2003).
/*!
This class both smooths the bed and computes coefficients for an approximation
to Schoof's \f$\theta\f$.  The factor \f$0\le \theta \le 1\f$ multiplies the diffusivity
in the theory of the effect of bed roughness in the SIA by Christian Schoof
(2003; <i>The effect of basal topography on ice sheet dynamics</i>) [\ref
Schoofbasaltopg2003].

For additional information on this class see page \ref bedrough.

The user of this class hands PISMBedSmoother an "original" topography, and it
is preprocessed to fill the smoothed topography \c topgsmooth, and the
coefficients in an approximation to \f$\theta\f$.  This is done by a call to 
\c preprocess_bed().  The call requires the half-width of the smoothing square
(a distance in m), or the number of grid points in each direction in the
smoothing rectangle, and the Glen exponent.

The call to \c preprocess_bed() <b>must be repeated</b> any time the "original"
topography changes, for instance at the start of an IceModel run, or at a bed
deformation step in an IceModel run.

PISMBedSmoother then provides three major functionalities, all of which \e must
\e follow the call to \c preprocess_bed():
-# User accesses public IceModelVec2S \c topgsmooth, the smoothed bed itself.
-# User asks \c get_smoothed_thk() for gridded values of the consistent smoothed
   version of the ice thickness, which is the thickness corresponding to a given
   surface elevation and the pre-computed smoothed bed.
-# User asks for gridded values of \f$\theta(h,x,y)\f$ using \c get_theta().

Here is a basic example of the creation and usage of a PISMBedSmoother instance.
Note that PISMBedSmoother will update ghosted values in \c topgsmooth, and other
internal fields, and update them in the return fields \c thksmooth, \c theta,
if asked.  In IceModel::velocitySIAStaggered(), MAX_GHOSTS=2
\code
    PISMBedSmoother smoother(grid, config, 2);
    const PetscReal n = 3.0,
                    lambda = 50.0e3;
    ierr = smoother.preprocess_bed(topg, n, lambda); CHKERRQ(ierr);
    ierr = smoother.get_smoothed_thk(usurf, thk, 1, &thksmooth); CHKERRQ(ierr);
    ierr = smoother.get_theta(usurf, n, 1, &theta); CHKERRQ(ierr);
\endcode
See IceGrid and NCConfigVariable documentation for initializing \c grid and
\c config.  Note we assume \c topg, \c usurf, \c thk, \c thksmooth, and \c theta
are all created IceModelVec2S instances.
 */
class PISMBedSmoother {
public:
  PISMBedSmoother(IceGrid &g, const NCConfigVariable &conf, PetscInt MAX_GHOSTS);
  virtual ~PISMBedSmoother();

  virtual PetscErrorCode preprocess_bed(IceModelVec2S topg, PetscReal n, 
                                        PetscReal lambda);
  virtual PetscErrorCode preprocess_bed(IceModelVec2S topg, PetscReal n,
                                        PetscInt Nx_in, PetscInt Ny_in);

  virtual PetscErrorCode get_smoothing_domain(PetscInt &Nx_out, PetscInt &Ny_out);
  virtual PetscInt       get_max_ghosts() { return maxGHOSTS; }

  virtual PetscErrorCode get_smoothed_thk(IceModelVec2S usurf, IceModelVec2S thk, IceModelVec2Mask mask,
                                          PetscInt GHOSTS, IceModelVec2S *thksmooth);
  virtual PetscErrorCode get_theta(IceModelVec2S usurf, PetscReal n,
                                   PetscInt GHOSTS, IceModelVec2S *theta);

  IceModelVec2S topgsmooth;  //!< smoothed bed elevation; publicly-available; set by calling preprocess_bed(); has ghosts with width get_max_ghosts()

protected:
  IceGrid &grid;
  const NCConfigVariable &config;
  IceModelVec2S maxtl,C2,C3,C4;

  PetscInt Nx,Ny;  //!< number of grid points to smooth over; e.g.
                   //!i=-Nx,-Nx+1,...,-1,0,1,...,Nx-1,Nx; note Nx>=1 and Ny>=1
                   //!always, unless lambda<=0
  PetscInt maxGHOSTS; //!< topg will be read, and topgsmooth will be created,
                      //!with this ghosting width; thksmooth and theta can be
                      //!filled at this ghosting level or less

  PetscErrorCode allocate();
  PetscErrorCode deallocate();

  Vec g2, g2natural;  //!< global Vecs used to transfer data to/from processor 0.
  VecScatter scatter; //!< VecScatter used to transfer data to/from processor 0.
  Vec topgp0,         //!< original bed elevation on processor 0
      topgsmoothp0,   //!< smoothed bed elevation on processor 0
      maxtlp0,        //!< maximum elevation at (i,j) of local topography (nearby patch)
      C2p0,C3p0,C4p0;

  PetscErrorCode smooth_the_bed_on_proc0();
  PetscErrorCode compute_coefficients_on_proc0(PetscReal n);
};

#endif	// __PISMBedSmoother_hh

