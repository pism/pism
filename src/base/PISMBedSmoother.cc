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

#include "PISMBedSmoother.hh"

PISMBedSmoother::PISMBedSmoother(IceGrid &g, const NCConfigVariable &conf)
    : grid(g), config(conf) {

  if (allocate() != 0) {
    PetscPrintf(grid.com, "PISMBedSmoother constructor: allocate() failed\n");
    PetscEnd();
  }
  
}

PISMBedSmoother::~PISMBedSmoother() {

  if (deallocate() != 0) {
    PetscPrintf(grid.com, "PISMBedSmoother destructor: deallocate() failed\n");
    PetscEnd();
  }

}

/* the following is from the PETSc FAQ page:

How do I collect all the values from a parallel PETSc vector into a vector 
on the zeroth processor?

    * Create the scatter context that will do the communication
          o VecScatterCreateToZero(v,&ctx,&w);
    * Actually do the communication; this can be done repeatedly as needed
          o VecScatterBegin(ctx,v,w,INSERT_VALUES,SCATTER_FORWARD);
          o VecScatterEnd(ctx,v,w,INSERT_VALUES,SCATTER_FORWARD);
    * Remember to free the scatter context when no longer needed
          o VecScatterDestroy(ctx);

Note that this simply concatenates in the parallel ordering of the vector.
If you are using a vector from DACreateGlobalVector() you likely want to 
first call DAGlobalToNaturalBegin/End() to scatter the original vector into 
the natural ordering in a new global vector before calling 
VecScatterBegin/End() to scatter the natural vector onto process 0.
*/


PetscErrorCode PISMBedSmoother::allocate() {
  PetscErrorCode ierr;

  ierr = DACreateGlobalVector(grid.da2, &g2); CHKERRQ(ierr);

  // note we want a global Vec but reordered in the natural ordering so when it
  // is scattered to proc zero it is not all messed up; see above
  ierr = DACreateNaturalVector(grid.da2, &g2natural); CHKERRQ(ierr);
  // next get scatter context *and* allocate one of Vecs on proc zero
  ierr = VecScatterCreateToZero(g2natural, &scatter, &topgp0); CHKERRQ(ierr);
      
  // allocate Vecs that live on proc 0
  ierr = VecDuplicate(topgp0,&topgsmoothp0); CHKERRQ(ierr);
  ierr = VecDuplicate(topgp0,&C2p0); CHKERRQ(ierr);
  ierr = VecDuplicate(topgp0,&C3p0); CHKERRQ(ierr);
  ierr = VecDuplicate(topgp0,&C4p0); CHKERRQ(ierr);
  ierr = VecDuplicate(topgp0,&C5p0); CHKERRQ(ierr);

  // allocate Vecs that live on all procs
  ierr = topgsmooth.create(grid, "topgsmooth", false); CHKERRQ(ierr);
  ierr = topgsmooth.set_attrs(
     "bed_smoother_storage", 
     "smoothed bed elevation, in Schoof (2003) bed roughness parameterization",
		 "m", ""); CHKERRQ(ierr);
  ierr = C2.create(grid, "C2bedsmooth", false); CHKERRQ(ierr);
  ierr = C2.set_attrs(
     "bed_smoother_storage", 
     "polynomial coeff of H^-2, in Schoof (2003) bed roughness parameterization",
		 "m4", ""); CHKERRQ(ierr);
  ierr = C3.create(grid, "C3bedsmooth", false); CHKERRQ(ierr);
  ierr = C3.set_attrs(
     "bed_smoother_storage", 
     "polynomial coeff of H^-3, in Schoof (2003) bed roughness parameterization",
		 "m5", ""); CHKERRQ(ierr);
  ierr = C4.create(grid, "C4bedsmooth", false); CHKERRQ(ierr);
  ierr = C4.set_attrs(
     "bed_smoother_storage", 
     "polynomial coeff of H^-4, in Schoof (2003) bed roughness parameterization",
		 "m6", ""); CHKERRQ(ierr);
  ierr = C5.create(grid, "C5bedsmooth", false); CHKERRQ(ierr);
  ierr = C5.set_attrs(
     "bed_smoother_storage", 
     "polynomial coeff of H^-5, in Schoof (2003) bed roughness parameterization",
		 "m7", ""); CHKERRQ(ierr);

  if (grid.rank == 0) {
    // complete allocation on rank 0 if needed; perhaps not!
  }

  return 0;
}


PetscErrorCode PISMBedSmoother::deallocate() {
  PetscErrorCode ierr;

  ierr = VecDestroy(g2); CHKERRQ(ierr);
  ierr = VecDestroy(g2natural); CHKERRQ(ierr);
  ierr = VecScatterDestroy(scatter); CHKERRQ(ierr);

  ierr = VecDestroy(topgp0); CHKERRQ(ierr);
  ierr = VecDestroy(topgsmoothp0); CHKERRQ(ierr);
  ierr = VecDestroy(C2p0); CHKERRQ(ierr);
  ierr = VecDestroy(C3p0); CHKERRQ(ierr);
  ierr = VecDestroy(C4p0); CHKERRQ(ierr);
  ierr = VecDestroy(C5p0); CHKERRQ(ierr);

  // no need to destroy topgsmooth,C2,C3,C4,C5, because their destructors do
  // that job

  return 0;
}


PetscErrorCode PISMBedSmoother::transfer_to_proc0(IceModelVec2S *source, Vec result) {
  PetscErrorCode ierr;

  ierr = source->copy_to(g2);

  ierr = DAGlobalToNaturalBegin(grid.da2, g2, INSERT_VALUES, g2natural);
    CHKERRQ(ierr);
  ierr =   DAGlobalToNaturalEnd(grid.da2, g2, INSERT_VALUES, g2natural);
    CHKERRQ(ierr);

  ierr = VecScatterBegin(scatter, g2natural, result, INSERT_VALUES, SCATTER_FORWARD);
    CHKERRQ(ierr);
  ierr =   VecScatterEnd(scatter, g2natural, result, INSERT_VALUES, SCATTER_FORWARD);
    CHKERRQ(ierr);

  return 0;
}


PetscErrorCode PISMBedSmoother::transfer_from_proc0(Vec source, IceModelVec2S *result) {
  PetscErrorCode ierr;

  ierr = VecScatterBegin(scatter, source, g2natural, INSERT_VALUES, SCATTER_REVERSE);
    CHKERRQ(ierr);
  ierr =   VecScatterEnd(scatter, source, g2natural, INSERT_VALUES, SCATTER_REVERSE);
    CHKERRQ(ierr);

  ierr = DANaturalToGlobalBegin(grid.da2, g2natural, INSERT_VALUES, g2);
    CHKERRQ(ierr);
  ierr =   DANaturalToGlobalEnd(grid.da2, g2natural, INSERT_VALUES, g2);
    CHKERRQ(ierr);

  ierr = result->copy_from(g2); CHKERRQ(ierr);

  return 0;
}


PetscErrorCode PISMBedSmoother::preprocess_bed(
                 IceModelVec2S *topg, 
                 PetscReal lambdax, PetscReal lambday) {
  PetscErrorCode ierr;

  ierr = transfer_to_proc0(topg, topgp0); CHKERRQ(ierr);

  if (grid.rank == 0) {
    // actually do smoothing on proc 0

    const PetscReal
      n   = 3,  // FIXME: hard-wired!
      k   = (n + 2) / n,
      cc2 = k * (2 * n + 2) / (2 * n),
      cc3 = cc2 * (3 * n + 2) / (3 * n),
      cc4 = cc3 * (4 * n + 2) / (4 * n),
      cc5 = cc4 * (5 * n + 2) / (5 * n);

    // WORK IN HERE
    
    C2.scale(cc2);
    C3.scale(cc3);
    C4.scale(cc4);
    C5.scale(cc5);
    
    // transfer back to all procs in topgsmooth,C2,C3,C4,C5
  }

  return 0;
}


/*!
Implements the strategy for computing \f$\theta(H,x,y)\f$ from previously-
stored coefficients, described on page \ref bedrough.

Specifically, \f$\theta = \omega^{-n}\f$ where \f$\omega\f$ is a local average
of a rational function of thickness, approximated here by a Taylor polynomial:
  \f[ \omega = \fint \left(1 - \frac{\tilde b(x_1,x_2,\xi_1,\xi_2)}{H}\right)^{-(n+2)/n}\,d\xi_1\,d\xi_2 \approx 1 + C_2 H^{-2} + \dots + C_5 H^{-5} \f]
and \f$\tilde b\f$ is the local bed topography, a function with mean zero.
The coefficients \f$C_2,\dots,C_5\f$ are precomputed by preprocess_bed().
 */
PetscErrorCode PISMBedSmoother::get_theta(IceModelVec2S thk, IceModelVec2S *theta) {
  PetscErrorCode ierr;

  const PetscReal n = 3;  // FIXME make configurable

  PetscScalar **mytheta;
  ierr = theta->get_array(mytheta);    CHKERRQ(ierr);
  ierr = thk.begin_access(); CHKERRQ(ierr);
  ierr = C2.begin_access(); CHKERRQ(ierr);
  ierr = C3.begin_access(); CHKERRQ(ierr);
  ierr = C4.begin_access(); CHKERRQ(ierr);
  ierr = C5.begin_access(); CHKERRQ(ierr);

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (thk(i,j) > 100.0) { // FIXME make configurable
        const PetscReal
          Hinv = 1.0 / thk(i,j),
          omega = 1.0 + Hinv * Hinv * 
                   ( C2(i,j) + Hinv * ( C3(i,j) + Hinv * 
                      ( C4(i,j) + Hinv * C5(i,j) ) ) );
        mytheta[i][j] = pow(omega,-n);
      } else {
        mytheta[i][j] = 1.0;  // FIXME = min_theta; make configurable
      }
      // now guarantee in [0,1]
      if (mytheta[i][j] > 1.0)  mytheta[i][j] = 1.0;
      if (mytheta[i][j] < 0.0)  mytheta[i][j] = 0.0;
    }
  }
  
  ierr = C5.end_access(); CHKERRQ(ierr);
  ierr = C4.end_access(); CHKERRQ(ierr);
  ierr = C3.end_access(); CHKERRQ(ierr);
  ierr = C2.end_access(); CHKERRQ(ierr);
  ierr = thk.end_access(); CHKERRQ(ierr);
  ierr = theta->end_access(); CHKERRQ(ierr);

  return 0;
}

