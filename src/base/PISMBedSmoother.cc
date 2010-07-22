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
  ierr = VecDuplicate(topgp0,&maxtlp0); CHKERRQ(ierr);
  ierr = VecDuplicate(topgp0,&C2p0); CHKERRQ(ierr);
  ierr = VecDuplicate(topgp0,&C3p0); CHKERRQ(ierr);
  ierr = VecDuplicate(topgp0,&C4p0); CHKERRQ(ierr);
  ierr = VecDuplicate(topgp0,&C5p0); CHKERRQ(ierr);

  // allocate Vecs that live on all procs
  ierr = topgsmooth.create(grid, "topgsmooth", false); CHKERRQ(ierr);
  ierr = topgsmooth.set_attrs(
     "bed_smoother_tool", 
     "smoothed bed elevation, in bed roughness parameterization",
     "m", ""); CHKERRQ(ierr);
  ierr = maxtl.create(grid, "maxtl", false); CHKERRQ(ierr);
  ierr = maxtl.set_attrs(
     "bed_smoother_tool", 
     "maximum elevation in local topography patch, in bed roughness parameterization",
     "m", ""); CHKERRQ(ierr);
  ierr = C2.create(grid, "C2bedsmooth", false); CHKERRQ(ierr);
  ierr = C2.set_attrs(
     "bed_smoother_tool", 
     "polynomial coeff of H^-2, in bed roughness parameterization",
     "m2", ""); CHKERRQ(ierr);
  ierr = C3.create(grid, "C3bedsmooth", false); CHKERRQ(ierr);
  ierr = C3.set_attrs(
     "bed_smoother_tool", 
     "polynomial coeff of H^-3, in bed roughness parameterization",
     "m3", ""); CHKERRQ(ierr);
  ierr = C4.create(grid, "C4bedsmooth", false); CHKERRQ(ierr);
  ierr = C4.set_attrs(
     "bed_smoother_tool", 
     "polynomial coeff of H^-4, in bed roughness parameterization",
     "m4", ""); CHKERRQ(ierr);
  ierr = C5.create(grid, "C5bedsmooth", false); CHKERRQ(ierr);
  ierr = C5.set_attrs(
     "bed_smoother_tool", 
     "polynomial coeff of H^-5, in bed roughness parameterization",
     "m5", ""); CHKERRQ(ierr);

  return 0;
}


PetscErrorCode PISMBedSmoother::deallocate() {
  PetscErrorCode ierr;

  ierr = VecDestroy(g2); CHKERRQ(ierr);
  ierr = VecDestroy(g2natural); CHKERRQ(ierr);
  ierr = VecScatterDestroy(scatter); CHKERRQ(ierr);

  ierr = VecDestroy(topgp0); CHKERRQ(ierr);
  ierr = VecDestroy(topgsmoothp0); CHKERRQ(ierr);
  ierr = VecDestroy(maxtlp0); CHKERRQ(ierr);
  ierr = VecDestroy(C2p0); CHKERRQ(ierr);
  ierr = VecDestroy(C3p0); CHKERRQ(ierr);
  ierr = VecDestroy(C4p0); CHKERRQ(ierr);
  ierr = VecDestroy(C5p0); CHKERRQ(ierr);

  // no need to destroy topgsmooth,maxtl,C2,C3,C4,C5; their destructors do it

  return 0;
}


PetscErrorCode PISMBedSmoother::transfer_to_proc0(IceModelVec2S source, Vec *result) {
  PetscErrorCode ierr;

  ierr = source.copy_to(g2);

  ierr = DAGlobalToNaturalBegin(grid.da2, g2, INSERT_VALUES, g2natural);
    CHKERRQ(ierr);
  ierr =   DAGlobalToNaturalEnd(grid.da2, g2, INSERT_VALUES, g2natural);
    CHKERRQ(ierr);

  ierr = VecScatterBegin(scatter, g2natural, *result, INSERT_VALUES, SCATTER_FORWARD);
    CHKERRQ(ierr);
  ierr =   VecScatterEnd(scatter, g2natural, *result, INSERT_VALUES, SCATTER_FORWARD);
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


/*!
Input lambda gives physical half-width (in m) of square over which to do the
average.  Only square smoothing domains are allowed.
 */
PetscErrorCode PISMBedSmoother::preprocess_bed(
                 IceModelVec2S topg, PetscReal lambda, PetscReal n) {
  PetscErrorCode ierr;

  if (lambda <= 0) {
    // smoothing completely inactive.  we will transfer the original bed topg
    // to public member topgsmooth ...
    ierr = topgsmooth.copy_from(topg); CHKERRQ(ierr);
    // and tell get_theta() to return theta=1
    Nx = -1;
    Ny = -1;
    return 0;
  }
    
  // determine Nx, Ny, which are always at least one if lambda > 0
  Nx = static_cast<PetscInt>(ceil(lambda / grid.dx));
  Ny = static_cast<PetscInt>(ceil(lambda / grid.dy));
  if (Nx < 1)  Nx = 1;
  if (Ny < 1)  Ny = 1;
  if ((Nx >= grid.Mx) || (Ny >= grid.My)) {
    SETERRQ(1,"PISM ERROR: lambda in bed smoother too large; domain of smoothing exceeds domain");
  }

  // manage the preprocessing by calling procedures which only act on proc 0
  ierr = transfer_to_proc0(topg, &topgp0); CHKERRQ(ierr);
  ierr = smooth_the_bed(); CHKERRQ(ierr);
  ierr = transfer_from_proc0(topgsmoothp0, &topgsmooth); CHKERRQ(ierr);
  ierr = compute_coefficients(n); CHKERRQ(ierr);
  ierr = transfer_from_proc0(maxtlp0, &maxtl); CHKERRQ(ierr);
  ierr = transfer_from_proc0(C2p0, &C2); CHKERRQ(ierr);
  ierr = transfer_from_proc0(C3p0, &C3); CHKERRQ(ierr);
  ierr = transfer_from_proc0(C4p0, &C4); CHKERRQ(ierr);
  ierr = transfer_from_proc0(C5p0, &C5); CHKERRQ(ierr);
  
  return 0;
}


//! Computes the smoothed bed by a simple average over a rectangle of grid points.
PetscErrorCode PISMBedSmoother::smooth_the_bed() {

  if (grid.rank == 0) {
    PetscErrorCode ierr;
    PetscScalar **b0, **bs;
    ierr = VecGetArray2d(topgp0,       grid.Mx, grid.My, 0, 0, &b0); CHKERRQ(ierr);
    ierr = VecGetArray2d(topgsmoothp0, grid.Mx, grid.My, 0, 0, &bs); CHKERRQ(ierr);

    for (PetscInt i=0; i < grid.Mx; i++) {
      for (PetscInt j=0; j < grid.My; j++) {
        // average only over those points which are in the grid; do not wrap
        //   periodically
        PetscReal sum   = 0.0;
        PetscInt  count = 0;
        for (PetscInt r = -Nx; r <= Nx; r++) {
          for (PetscInt s = -Ny; s <= Ny; s++) {
            if ((i+r >= 0) && (i+r < grid.Mx) && (j+s >= 0) && (j+s < grid.My)) {
              sum += b0[i+r][j+s];
              count++;
            }
          }
        }
        bs[i][j] = sum / static_cast<PetscReal>(count);
      }
    }

    ierr = VecRestoreArray2d(topgsmoothp0, grid.Mx, grid.My, 0, 0, &bs); CHKERRQ(ierr);
    ierr = VecRestoreArray2d(topgp0,       grid.Mx, grid.My, 0, 0, &b0); CHKERRQ(ierr);
  }

  return 0;
}


PetscErrorCode PISMBedSmoother::compute_coefficients(PetscReal n) {

  if (grid.rank == 0) {
    PetscErrorCode ierr;
    PetscScalar **b0, **bs, **c2, **c3, **c4, **c5, **mt;
    ierr = VecGetArray2d(topgp0,       grid.Mx, grid.My, 0, 0, &b0); CHKERRQ(ierr);
    ierr = VecGetArray2d(topgsmoothp0, grid.Mx, grid.My, 0, 0, &bs); CHKERRQ(ierr);
    ierr = VecGetArray2d(maxtlp0,      grid.Mx, grid.My, 0, 0, &mt); CHKERRQ(ierr);
    ierr = VecGetArray2d(C2p0,         grid.Mx, grid.My, 0, 0, &c2); CHKERRQ(ierr);
    ierr = VecGetArray2d(C3p0,         grid.Mx, grid.My, 0, 0, &c3); CHKERRQ(ierr);
    ierr = VecGetArray2d(C4p0,         grid.Mx, grid.My, 0, 0, &c4); CHKERRQ(ierr);
    ierr = VecGetArray2d(C5p0,         grid.Mx, grid.My, 0, 0, &c5); CHKERRQ(ierr);

    for (PetscInt i=0; i < grid.Mx; i++) {
      for (PetscInt j=0; j < grid.My; j++) {
        // average only over those points which are in the grid
        // do not wrap periodically
        PetscReal topgs     = bs[i][j],
                  maxtltemp = 0.0,
                  sum2      = 0.0,
                  sum3      = 0.0,
                  sum4      = 0.0,
                  sum5      = 0.0;
        PetscInt  count     = 0;
        for (PetscInt r = -Nx; r <= Nx; r++) {
          for (PetscInt s = -Ny; s <= Ny; s++) {
            if ((i+r >= 0) && (i+r < grid.Mx) && (j+s >= 0) && (j+s < grid.My)) {
              // tl is elevation of local topography at a pt in patch
              const PetscReal tl  = b0[i+r][j+s] - topgs;  
              maxtltemp = PetscMax(maxtltemp, tl);
              // accumulate 2nd, 3rd, 4th, and 5th powers with only 4 mults
              const PetscReal tl2 = tl * tl, tl4 = tl2 * tl2;
              sum2 += tl2;
              sum3 += tl2 * tl;
              sum4 += tl4;
              sum5 += tl4 * tl;
              count++;
            }
          }
        }
        mt[i][j] = maxtltemp;
        // unprotected division by count but r=0,s=0 case guarantees count>=1
        c2[i][j] = sum2 / static_cast<PetscReal>(count);
        c3[i][j] = sum3 / static_cast<PetscReal>(count);
        c4[i][j] = sum4 / static_cast<PetscReal>(count);
        c5[i][j] = sum5 / static_cast<PetscReal>(count);
      }
    }

    ierr = VecRestoreArray2d(C5p0,         grid.Mx, grid.My, 0, 0, &c5); CHKERRQ(ierr);
    ierr = VecRestoreArray2d(C4p0,         grid.Mx, grid.My, 0, 0, &c4); CHKERRQ(ierr);
    ierr = VecRestoreArray2d(C3p0,         grid.Mx, grid.My, 0, 0, &c3); CHKERRQ(ierr);
    ierr = VecRestoreArray2d(C2p0,         grid.Mx, grid.My, 0, 0, &c2); CHKERRQ(ierr);
    ierr = VecRestoreArray2d(maxtlp0,      grid.Mx, grid.My, 0, 0, &mt); CHKERRQ(ierr);
    ierr = VecRestoreArray2d(topgsmoothp0, grid.Mx, grid.My, 0, 0, &bs); CHKERRQ(ierr);
    ierr = VecRestoreArray2d(topgp0,       grid.Mx, grid.My, 0, 0, &b0); CHKERRQ(ierr);

    // scale the coeffs in Taylor series
    const PetscReal
      k  = (n + 2) / n,
      s2 = k * (2 * n + 2) / (2 * n),
      s3 = s2 * (3 * n + 2) / (3 * n),
      s4 = s3 * (4 * n + 2) / (4 * n),
      s5 = s4 * (5 * n + 2) / (5 * n);
    ierr = VecScale(C2p0,s2); CHKERRQ(ierr);
    ierr = VecScale(C3p0,s3); CHKERRQ(ierr);
    ierr = VecScale(C4p0,s4); CHKERRQ(ierr);
    ierr = VecScale(C5p0,s5); CHKERRQ(ierr);
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
PetscErrorCode PISMBedSmoother::get_theta(
      IceModelVec2S thk, PetscReal n, IceModelVec2S *theta) {
  PetscErrorCode ierr;

  if ((Nx <=0) || (Ny <= 0)) {
    ierr = theta->set(1.0); CHKERRQ(ierr);
    return 0;
  }

  PetscScalar **mytheta;
  ierr = theta->get_array(mytheta); CHKERRQ(ierr);
  ierr = thk.begin_access(); CHKERRQ(ierr);
  ierr = maxtl.begin_access(); CHKERRQ(ierr);
  ierr = C2.begin_access(); CHKERRQ(ierr);
  ierr = C3.begin_access(); CHKERRQ(ierr);
  ierr = C4.begin_access(); CHKERRQ(ierr);
  ierr = C5.begin_access(); CHKERRQ(ierr);

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (thk(i,j) > maxtl(i,j)) { 
        // if thickness exceeds maximum variation in patch of local topography,
        // so ice buries local topography; note maxtl >= 0 always
        const PetscReal
          Hinv = 1.0 / thk(i,j),
          omega = 1.0 + Hinv * Hinv * 
                   ( C2(i,j) + Hinv * ( C3(i,j) + Hinv * 
                      ( C4(i,j) + Hinv * C5(i,j) ) ) );
        if (omega <= 0) {
          SETERRQ2(1,"PISM ERROR: omega is negative for i=%d,j=%d\n"
                     "    in PISMBedSmoother.get_theta() ... ending\n",i,j);
        }
        mytheta[i][j] = pow(omega,-n);
      } else {
        mytheta[i][j] = 0.05;  // FIXME = min_theta; make configurable
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
  ierr = maxtl.end_access(); CHKERRQ(ierr);
  ierr = thk.end_access(); CHKERRQ(ierr);
  ierr = theta->end_access(); CHKERRQ(ierr);

  return 0;
}

