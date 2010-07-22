// Copyright (C) 2010 Ed Bueler
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
This class both smooths the bed and computes coefficients which allow the computation of a factor \f$0\le \theta \le 1\f$ which multiplies the diffusivity in the SIA.  The theory is from Christian Schoof's (2003) <i>The effect of basal topography on ice sheet dynamics</i> [\ref Schoofbasaltopg2003].

PISMBedSmoother provides three functionalities, besides construction and
destruction:
-# User hands PISMBedSmoother an "original" topography, and it is preprocessed.
   This call to preprocess_bed() requires the half-width of the smoothing square
   (a distance in m) and the Glen exponent.
-# User asks for gridded values of \f$\theta(H,x,y)\f$ using get_theta().  This
   call must follow a call to preprocess_bed().
-# User accesses public IceModelVec2S topgsmooth, the smoothed bed itself.
   This access must also follow a call to preprocess_bed().
 */
class PISMBedSmoother {
public:
  PISMBedSmoother(IceGrid &g, const NCConfigVariable &conf);
  virtual ~PISMBedSmoother();

  virtual PetscErrorCode preprocess_bed(IceModelVec2S topg, PetscReal lambda, PetscReal n);

  virtual PetscErrorCode get_theta(IceModelVec2S thk, PetscReal n, IceModelVec2S *theta);

  IceModelVec2S topgsmooth;  //!< smoothed bed elevation; publicly-available; set by calling preprocess_bed()

protected:
  IceGrid &grid;
  const NCConfigVariable &config;
  IceModelVec2S maxtl,C2,C3,C4,C5;
  PetscInt Nx,Ny;  //!< number of grid points to smooth over; e.g. i=-Nx,-Nx+1,...,-1,0,1,...,Nx-1,Nx; note Nx>=1 and Ny>=1 always, unless lambda<=0

  PetscErrorCode smooth_the_bed();
  PetscErrorCode compute_coefficients(PetscReal n);

  PetscErrorCode allocate();
  PetscErrorCode deallocate();
  PetscErrorCode transfer_to_proc0(IceModelVec2S source, Vec *result);
  PetscErrorCode transfer_from_proc0(Vec source, IceModelVec2S *result);
  Vec g2, g2natural;  //!< global Vecs used to transfer data to/from processor 0.
  VecScatter scatter; //!< VecScatter used to transfer data to/from processor 0.
  Vec topgp0,         //!< original bed elevation on processor 0
      topgsmoothp0,   //!< smoothed bed elevation on processor 0
      maxtlp0,        //!< maximum elevation at (i,j) of local topography (nearby patch)
      C2p0,C3p0,C4p0,C5p0;
};

#endif	// __PISMBedSmoother_hh

