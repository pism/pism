// Copyright (C) 2010, 2011, 2012, 2013 Ed Bueler and Constantine Khroulev
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
#include "Mask.hh"


PISMBedSmoother::PISMBedSmoother(
                   IceGrid &g, const NCConfigVariable &conf, PetscInt MAX_GHOSTS)
    : grid(g), config(conf), maxGHOSTS(MAX_GHOSTS) {

  if (allocate() != 0) {
    PetscPrintf(grid.com, "PISMBedSmoother constructor: allocate() failed\n");
    PISMEnd();
  }

  if (config.get("bed_smoother_range") > 0.0) {
    verbPrintf(2, grid.com, 
               "* Initializing bed smoother object with %.3f km half-width ...\n",
               config.get("bed_smoother_range") / 1000.0);
  }

}


PISMBedSmoother::~PISMBedSmoother() {

  if (deallocate() != 0) {
    PetscPrintf(grid.com, "PISMBedSmoother destructor: deallocate() failed\n");
    PISMEnd();
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
  DM da2;
  ierr = grid.get_dm(1, grid.max_stencil_width, da2); CHKERRQ(ierr);

  ierr = DMCreateGlobalVector(da2, &g2); CHKERRQ(ierr);

  // note we want a global Vec but reordered in the natural ordering so when it
  // is scattered to proc zero it is not all messed up; see above
  ierr = DMDACreateNaturalVector(da2, &g2natural); CHKERRQ(ierr);
  // next get scatter context *and* allocate one of Vecs on proc zero
  ierr = VecScatterCreateToZero(g2natural, &scatter, &topgp0); CHKERRQ(ierr);
      
  // allocate Vecs that live on proc 0
  ierr = VecDuplicate(topgp0,&topgsmoothp0); CHKERRQ(ierr);
  ierr = VecDuplicate(topgp0,&maxtlp0); CHKERRQ(ierr);
  ierr = VecDuplicate(topgp0,&C2p0); CHKERRQ(ierr);
  ierr = VecDuplicate(topgp0,&C3p0); CHKERRQ(ierr);
  ierr = VecDuplicate(topgp0,&C4p0); CHKERRQ(ierr);

  // allocate Vecs that live on all procs; all have to be as "wide" as any of
  //   their prospective uses
  ierr = topgsmooth.create(grid, "topgsmooth", true, maxGHOSTS); CHKERRQ(ierr);
  ierr = topgsmooth.set_attrs(
     "bed_smoother_tool", 
     "smoothed bed elevation, in bed roughness parameterization",
     "m", ""); CHKERRQ(ierr);
  ierr = maxtl.create(grid, "maxtl", true, maxGHOSTS); CHKERRQ(ierr);
  ierr = maxtl.set_attrs(
     "bed_smoother_tool", 
     "maximum elevation in local topography patch, in bed roughness parameterization",
     "m", ""); CHKERRQ(ierr);
  ierr = C2.create(grid, "C2bedsmooth", true, maxGHOSTS); CHKERRQ(ierr);
  ierr = C2.set_attrs(
     "bed_smoother_tool", 
     "polynomial coeff of H^-2, in bed roughness parameterization",
     "m2", ""); CHKERRQ(ierr);
  ierr = C3.create(grid, "C3bedsmooth", true, maxGHOSTS); CHKERRQ(ierr);
  ierr = C3.set_attrs(
     "bed_smoother_tool", 
     "polynomial coeff of H^-3, in bed roughness parameterization",
     "m3", ""); CHKERRQ(ierr);
  ierr = C4.create(grid, "C4bedsmooth", true, maxGHOSTS); CHKERRQ(ierr);
  ierr = C4.set_attrs(
     "bed_smoother_tool", 
     "polynomial coeff of H^-4, in bed roughness parameterization",
     "m4", ""); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode PISMBedSmoother::deallocate() {
  PetscErrorCode ierr;

  ierr = VecDestroy(&g2); CHKERRQ(ierr);
  ierr = VecDestroy(&g2natural); CHKERRQ(ierr);
  ierr = VecScatterDestroy(&scatter); CHKERRQ(ierr);

  ierr = VecDestroy(&topgp0); CHKERRQ(ierr);
  ierr = VecDestroy(&topgsmoothp0); CHKERRQ(ierr);
  ierr = VecDestroy(&maxtlp0); CHKERRQ(ierr);
  ierr = VecDestroy(&C2p0); CHKERRQ(ierr);
  ierr = VecDestroy(&C3p0); CHKERRQ(ierr);
  ierr = VecDestroy(&C4p0); CHKERRQ(ierr);
  // no need to destroy topgsmooth,maxtl,C2,C3,C4; their destructors do it
  return 0;
}


/*!
Input lambda gives physical half-width (in m) of square over which to do the
average.  Only square smoothing domains are allowed with this call, which is the
default case.
 */
PetscErrorCode PISMBedSmoother::preprocess_bed(
                 IceModelVec2S topg, PetscReal n, PetscReal lambda) {
  PetscErrorCode ierr;

  if (lambda <= 0.0) {
    // smoothing completely inactive.  we transfer the original bed topg,
    //   including ghosts, to public member topgsmooth ...
    ierr = topg.update_ghosts(topgsmooth); CHKERRQ(ierr);
    // and we tell get_theta() to return theta=1
    Nx = -1;
    Ny = -1;
    return 0;
  }
    
  // determine Nx, Ny, which are always at least one if lambda > 0
  Nx = static_cast<PetscInt>(ceil(lambda / grid.dx));
  Ny = static_cast<PetscInt>(ceil(lambda / grid.dy));
  if (Nx < 1)  Nx = 1;
  if (Ny < 1)  Ny = 1;
  //PetscPrintf(grid.com,"PISMBedSmoother:  Nx = %d, Ny = %d\n",Nx,Ny);

  ierr = preprocess_bed(topg, n, Nx, Ny); CHKERRQ(ierr);
  return 0;
}


/*!
Inputs Nx,Ny gives half-width in number of grid points, over which to do the
average.
 */
PetscErrorCode PISMBedSmoother::preprocess_bed(
                 IceModelVec2S topg, PetscReal n, PetscInt Nx_in, PetscInt Ny_in) {
  PetscErrorCode ierr;

  if ((Nx_in >= grid.Mx) || (Ny_in >= grid.My)) {
    SETERRQ(grid.com, 1,"PISM ERROR: input Nx, Ny in bed smoother is too large because\n"
              "            domain of smoothing exceeds IceGrid domain\n");
  }
  Nx = Nx_in; Ny = Ny_in;

  if ((Nx < 0) || (Ny < 0)) {
    // smoothing completely inactive.  we transfer the original bed topg,
    //   including ghosts, to public member topgsmooth ...
    ierr = topg.update_ghosts(topgsmooth); CHKERRQ(ierr);
    return 0;
  }

  ierr = topg.put_on_proc0(topgp0, scatter, g2, g2natural); CHKERRQ(ierr);
  ierr = smooth_the_bed_on_proc0(); CHKERRQ(ierr);
  // next call *does indeed* fill ghosts in topgsmooth
  ierr = topgsmooth.get_from_proc0(topgsmoothp0, scatter, g2, g2natural); CHKERRQ(ierr);
  
  ierr = compute_coefficients_on_proc0(n); CHKERRQ(ierr);
  // following calls *do* fill the ghosts
  ierr = maxtl.get_from_proc0(maxtlp0, scatter, g2, g2natural); CHKERRQ(ierr);
  ierr = C2.get_from_proc0(C2p0, scatter, g2, g2natural); CHKERRQ(ierr);
  ierr = C3.get_from_proc0(C3p0, scatter, g2, g2natural); CHKERRQ(ierr);
  ierr = C4.get_from_proc0(C4p0, scatter, g2, g2natural); CHKERRQ(ierr);
  return 0;
}


/*!
Call preprocess_bed() first.
 */
PetscErrorCode PISMBedSmoother::get_smoothing_domain(PetscInt &Nx_out, PetscInt &Ny_out) {
  Nx_out = Nx;
  Ny_out = Ny;
  return 0;
}


//! Computes the smoothed bed by a simple average over a rectangle of grid points.
PetscErrorCode PISMBedSmoother::smooth_the_bed_on_proc0() {

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


PetscErrorCode PISMBedSmoother::compute_coefficients_on_proc0(PetscReal n) {

  if (grid.rank == 0) {
    PetscErrorCode ierr;
    PetscScalar **b0, **bs, **c2, **c3, **c4, **mt;
    ierr = VecGetArray2d(topgp0,       grid.Mx, grid.My, 0, 0, &b0); CHKERRQ(ierr);
    ierr = VecGetArray2d(topgsmoothp0, grid.Mx, grid.My, 0, 0, &bs); CHKERRQ(ierr);
    ierr = VecGetArray2d(maxtlp0,      grid.Mx, grid.My, 0, 0, &mt); CHKERRQ(ierr);
    ierr = VecGetArray2d(C2p0,         grid.Mx, grid.My, 0, 0, &c2); CHKERRQ(ierr);
    ierr = VecGetArray2d(C3p0,         grid.Mx, grid.My, 0, 0, &c3); CHKERRQ(ierr);
    ierr = VecGetArray2d(C4p0,         grid.Mx, grid.My, 0, 0, &c4); CHKERRQ(ierr);

    for (PetscInt i=0; i < grid.Mx; i++) {
      for (PetscInt j=0; j < grid.My; j++) {
        // average only over those points which are in the grid
        // do not wrap periodically
        PetscReal topgs     = bs[i][j],
                  maxtltemp = 0.0,
                  sum2      = 0.0,
                  sum3      = 0.0,
                  sum4      = 0.0;
        PetscInt  count     = 0;
        for (PetscInt r = -Nx; r <= Nx; r++) {
          for (PetscInt s = -Ny; s <= Ny; s++) {
            if ((i+r >= 0) && (i+r < grid.Mx) && (j+s >= 0) && (j+s < grid.My)) {
              // tl is elevation of local topography at a pt in patch
              const PetscReal tl  = b0[i+r][j+s] - topgs;  
              maxtltemp = PetscMax(maxtltemp, tl);
              // accumulate 2nd, 3rd, and 4th powers with only 3 mults
              const PetscReal tl2 = tl * tl;
              sum2 += tl2;
              sum3 += tl2 * tl;
              sum4 += tl2 * tl2;
              count++;
            }
          }
        }
        mt[i][j] = maxtltemp;
        // unprotected division by count but r=0,s=0 case guarantees count>=1
        c2[i][j] = sum2 / static_cast<PetscReal>(count);
        c3[i][j] = sum3 / static_cast<PetscReal>(count);
        c4[i][j] = sum4 / static_cast<PetscReal>(count);
      }
    }

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
      s4 = s3 * (4 * n + 2) / (4 * n);
    ierr = VecScale(C2p0,s2); CHKERRQ(ierr);
    ierr = VecScale(C3p0,s3); CHKERRQ(ierr);
    ierr = VecScale(C4p0,s4); CHKERRQ(ierr);
  }

  return 0;
}


//! Computes a smoothed thickness map.
/*!
The result \c thksmooth is the difference between the given upper surface
elevation (\c usurf) and the stored smoothed bed topography (\c topgsmooth),
except where the given original thickness (\c thk) is zero.  In places where
the original thickness is zero, the result \c thksmooth is also set to zero.

Ghosted values are updated directly and no communication occurs.  In fact,
we \e assume \c usurf, \c thk, and \c thksmooth all have stencil width at least
equal to GHOSTS.  We \e check whether \c topgsmooth, which has stencil width 
maxGHOSTS, has at least GHOSTS stencil width, and throw an error if not.

Call preprocess_bed() first.
 */
PetscErrorCode PISMBedSmoother::get_smoothed_thk(IceModelVec2S usurf,
                                                 IceModelVec2S thk,
                                                 IceModelVec2Int mask,
                                                 PetscInt GHOSTS,
                                                 IceModelVec2S *thksmooth) { 
  PetscErrorCode ierr;  

  MaskQuery M(mask);

  if (GHOSTS > maxGHOSTS) {
    SETERRQ2(grid.com, 1,"PISM ERROR:  PISMBedSmoother fields do not have stencil\n"
               "  width sufficient to fill thksmooth with GHOSTS=%d;\n"
               "  construct PISMBedSmoother with MAX_GHOSTS>=%d\n",
               GHOSTS,GHOSTS);
  }

  PetscScalar **thks;  
  ierr = mask.begin_access(); CHKERRQ(ierr);
  ierr = topgsmooth.begin_access(); CHKERRQ(ierr);
  ierr = maxtl.begin_access(); CHKERRQ(ierr);
  ierr = usurf.begin_access(); CHKERRQ(ierr);
  ierr = thk.begin_access(); CHKERRQ(ierr);
  ierr = thksmooth->get_array(thks); CHKERRQ(ierr);
  for (PetscInt i = grid.xs - GHOSTS; i < grid.xs+grid.xm + GHOSTS; ++i) {
    for (PetscInt j = grid.ys - GHOSTS; j < grid.ys+grid.ym + GHOSTS; ++j) {
      if (thk(i,j) < 0.0) {
        SETERRQ2(grid.com, 2,
          "PISM ERROR:  PISMBedSmoother detects negative original thickness\n"
          "  at location (i,j) = (%d,%d) ... ending\n",i,j);
      } else if (thk(i,j) == 0.0) {
        thks[i][j] = 0.0;
      } else if (maxtl(i,j) >= thk(i,j)) {
        thks[i][j] = thk(i,j);
      } else {
        if (M.grounded(i,j)) {
          // if grounded, compute smoothed thickness as the difference of ice
          // surface elevation and smoothed bed elevation
          const PetscScalar thks_try = usurf(i,j) - topgsmooth(i,j);
          thks[i][j] = (thks_try > 0.0) ? thks_try : 0.0;
        } else {
          // if floating, use original thickness (note: surface elevation was
          // computed using this thickness and the sea level elevation)
          thks[i][j] = thk(i,j);
        }
      }
    }
  }
  ierr = topgsmooth.end_access(); CHKERRQ(ierr);
  ierr = maxtl.end_access(); CHKERRQ(ierr);
  ierr = usurf.end_access(); CHKERRQ(ierr);
  ierr = thk.end_access(); CHKERRQ(ierr);
  ierr = thksmooth->end_access(); CHKERRQ(ierr);
  ierr = mask.end_access(); CHKERRQ(ierr);

  return 0;
}


/*!
Implements the strategy for computing \f$\theta(h,x,y)\f$ from previously-
stored coefficients, described on page \ref bedrough and in [\ref 
Schoofbasaltopg2003].

Specifically, \f$\theta = \omega^{-n}\f$ where \f$\omega\f$ is a local average
of a rational function of surface elevation, approximated here by a Taylor polynomial:
  \f[ \omega = \fint \left(1 - \frac{\tilde b(x_1,x_2,\xi_1,\xi_2)}{H}
                           \right)^{-(n+2)/n}\,d\xi_1\,d\xi_2
             \approx 1 + C_2 H^{-2} + C_3 H^{-3} + C_4 H^{-4} \f]
where \f$h =\f$ usurf, \f$H = h -\f$ topgsmooth and \f$\tilde b\f$ is the local
bed topography, a function with mean zero.  The coefficients \f$C_2,C_3,C_4\f$,
which depend on \f$x,y\f$, are precomputed by \c preprocess_bed().

Ghosted values are updated directly and no communication occurs.  In fact,
we \e assume \c usurf and \c theta have stencil width at least
equal to GHOSTS.  We \e check whether \c topgsmooth, which has stencil width 
maxGHOSTS, has at least GHOSTS stencil width, and throw an error if not.

Call preprocess_bed() first.
 */
PetscErrorCode PISMBedSmoother::get_theta(
      IceModelVec2S usurf, PetscReal n,
      PetscInt GHOSTS, IceModelVec2S *theta) {
  PetscErrorCode ierr;

  if (GHOSTS > maxGHOSTS) {
    SETERRQ2(grid.com, 1,
"PISM ERROR:  PISMBedSmoother::topgsmooth,maxtl,C2,C3,C4 do not have stencil\n"
"  width sufficient to fill theta with GHOSTS=%d;  construct PISMBedSmoother\n"
"  with MAX_GHOSTS>=%d\n", GHOSTS,GHOSTS);
  }
  
  if ((Nx < 0) || (Ny < 0)) {
    ierr = theta->set(1.0); CHKERRQ(ierr);
    return 0;
  }

  PetscScalar **mytheta;
  ierr = theta->get_array(mytheta); CHKERRQ(ierr);
  ierr = usurf.begin_access(); CHKERRQ(ierr);
  ierr = topgsmooth.begin_access(); CHKERRQ(ierr);
  ierr = maxtl.begin_access(); CHKERRQ(ierr);
  ierr = C2.begin_access(); CHKERRQ(ierr);
  ierr = C3.begin_access(); CHKERRQ(ierr);
  ierr = C4.begin_access(); CHKERRQ(ierr);
  for (PetscInt i = grid.xs - GHOSTS; i < grid.xs+grid.xm + GHOSTS; ++i) {
    for (PetscInt j = grid.ys - GHOSTS; j < grid.ys+grid.ym + GHOSTS; ++j) {
      const PetscScalar H = usurf(i,j) - topgsmooth(i,j);
      if (H > maxtl(i,j)) { 
        // thickness exceeds maximum variation in patch of local topography,
        // so ice buries local topography; note maxtl >= 0 always
        const PetscReal Hinv = 1.0 / PetscMax(H, 1.0);
        PetscReal omega = 1.0 + Hinv*Hinv * ( C2(i,j) + Hinv * ( C3(i,j) + Hinv*C4(i,j) ) );
        if (omega <= 0) {  // this check *should not* be necessary: p4(s) > 0
          SETERRQ2(grid.com, 1,"PISM ERROR: omega is negative for i=%d,j=%d\n"
                     "    in PISMBedSmoother.get_theta() ... ending\n",i,j);
        }

        if (omega < 0.001)      // this check *should not* be necessary
          omega = 0.001;

        mytheta[i][j] = pow(omega,-n);
        // now guarantee in [0,1]; this check *should not* be necessary, by convexity of p4
        if (mytheta[i][j] > 1.0)  mytheta[i][j] = 1.0;
        if (mytheta[i][j] < 0.0)  mytheta[i][j] = 0.0;
      } else {
        mytheta[i][j] = 0.00;  // FIXME = min_theta; make configurable
      }
    }
  }  
  ierr = C4.end_access(); CHKERRQ(ierr);
  ierr = C3.end_access(); CHKERRQ(ierr);
  ierr = C2.end_access(); CHKERRQ(ierr);
  ierr = maxtl.end_access(); CHKERRQ(ierr);
  ierr = topgsmooth.end_access(); CHKERRQ(ierr);
  ierr = usurf.end_access(); CHKERRQ(ierr);
  ierr = theta->end_access(); CHKERRQ(ierr);

  return 0;
}

