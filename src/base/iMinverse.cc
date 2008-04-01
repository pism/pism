// Copyright (C) 2008 Ed Bueler
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

#include <cstring>
#include <petscda.h>
#include <netcdf.h>
#include "nc_util.hh"
#include "iceModel.hh"


//! Invert ice velocities saved in NetCDF file to find till properties, based on user options.
/*!
Reads user options <tt>-cbar_to_till foo.nc</tt> and <tt>-csurf_to_till foo.nc</tt>.

Must be called after initialization is more-or-less completed.
Uses geometry, temperature field, and map of effective thickness of basal till
water in its inversion.

Expects to find variable \c cbar or \c csurf in NetCDF file \c foo.nc (according to option).  
If found, reads it.  Then calls, in turn, removeVerticalPlaneShearRateFromC(),
to compute the basal sliding, computeBasalShearFromSSA() to compute the basal shear stress, 
and finally computeTFAFromBasalShearStressUsingPseudoPlastic().
The result is to get the till yield stress and the till friction angle in a 
pseudo-plastic till model.  (This model includes purely-plastic and linearly-viscous cases).
A forward time-stepping model could then be based on the resulting till friction angle field.
 */
PetscErrorCode IceModel::invertVelocitiesFromNetCDF() {
  PetscErrorCode ierr;
  
  PetscTruth cbarTillSet, csurfTillSet;
  char cFile[PETSC_MAX_PATH_LEN];
  ierr = PetscOptionsGetString(PETSC_NULL, "-cbar_to_till", cFile, 
                               PETSC_MAX_PATH_LEN, &cbarTillSet); CHKERRQ(ierr);
  ierr = PetscOptionsGetString(PETSC_NULL, "-csurf_to_till", cFile, 
                               PETSC_MAX_PATH_LEN, &csurfTillSet); CHKERRQ(ierr);
  if ((cbarTillSet == PETSC_TRUE) && (csurfTillSet == PETSC_TRUE)) {
    SETERRQ(1,"both -cbar_to_till and -csurf_to_till set; NOT ALLOWED\n");
  }
  if ((cbarTillSet == PETSC_FALSE) && (csurfTillSet == PETSC_FALSE)) {
    return 0;
  }

  PetscInt fileExists = 0;
  int ncid,  stat;
  if (grid.rank == 0) {
    stat = nc_open(cFile, 0, &ncid); fileExists = (stat == NC_NOERR);
  }
  MPI_Bcast(&fileExists, 1, MPI_INT, 0, grid.com);  
  if (!fileExists) {
    SETERRQ(2,"NetCDF file with velocities (either cbar or csurf) not found.\n");
  }

  PetscInt cbarExists = 0, csurfExists = 0;
  int v_cbar, v_csurf;
  if (grid.rank == 0) {
    stat = nc_inq_varid(ncid, "cbar", &v_cbar); cbarExists = stat == NC_NOERR;
    stat = nc_inq_varid(ncid, "csurf", &v_csurf); csurfExists = stat == NC_NOERR;
  }
  ierr = MPI_Bcast(&cbarExists, 1, MPI_INT, 0, grid.com); CHKERRQ(ierr);
  ierr = MPI_Bcast(&csurfExists, 1, MPI_INT, 0, grid.com); CHKERRQ(ierr);

  if ((cbarTillSet == PETSC_TRUE) && (!cbarExists)) {
    SETERRQ1(3,"-cbar_to_till set but cbar not found in %s\n", cFile);
  }
  if ((csurfTillSet == PETSC_TRUE) && (!csurfExists)) {
    SETERRQ1(4,"-csurf_to_till set but csurf not found in %s\n", cFile);
  }

  // our goal is to create "local interpolation context" from dimensions, 
  //   limits, and lengths; compare code in bootstrapFromFile_netCDF()
  size_t dim[5];
  double bdy[7];
  ierr = nct.get_dims_limits_lengths_2d(ncid, dim, bdy, grid.com); CHKERRQ(ierr);
  dim[3] = 1; 
  dim[4] = 1;
  bdy[5] = 0.0;
  bdy[6] = 0.0;
  MPI_Bcast(dim, 5, MPI_LONG, 0, grid.com);
  MPI_Bcast(bdy, 7, MPI_DOUBLE, 0, grid.com);
  double *z_bif, *zb_bif;
  z_bif = new double[dim[3]];
  zb_bif = new double[dim[4]];
  z_bif[0] = 0.0;
  zb_bif[0] = 0.0;
  LocalInterpCtx lic(ncid, dim, bdy, z_bif, zb_bif, grid);
  delete z_bif;
  delete zb_bif;
  //DEBUG:  ierr = lic.printGrid(grid.com); CHKERRQ(ierr);

  ierr = verbPrintf(2, grid.com, 
     "reading variable %s from file %s; computing till yield stress and till friction\n"
     "  angle by removing SIA vertical shear and computing membrane stresses\n",
     (cbarTillSet == PETSC_TRUE) ? "cbar" : "csurf", cFile); CHKERRQ(ierr);

  Vec cIn;
  ierr = VecDuplicate(vh, &cIn); CHKERRQ(ierr);
  if ((cbarTillSet == PETSC_TRUE) && (cbarExists)) {
    ierr = nct.regrid_local_var("cbar", 2, lic, grid, grid.da2, cIn, g2, false);
             CHKERRQ(ierr);
  } else if ((csurfTillSet == PETSC_TRUE) && (csurfExists)) {
    ierr = nct.regrid_local_var("csurf", 2, lic, grid, grid.da2, cIn, g2, false);
             CHKERRQ(ierr);
  } else {
    SETERRQ(998,"how did I get here?");
  }
  if (grid.rank == 0) {
    stat = nc_close(ncid);
  }

  Vec taubxComputed, taubyComputed;
  ierr = VecDuplicate(vh, &taubxComputed); CHKERRQ(ierr);
  ierr = VecDuplicate(vh, &taubyComputed); CHKERRQ(ierr);

  ierr = removeVerticalPlaneShearRateFromC(csurfTillSet, cIn, vubarSSA, vvbarSSA);
            CHKERRQ(ierr);
  ierr = computeBasalShearFromSSA(vubarSSA, vvbarSSA, taubxComputed, taubyComputed);
            CHKERRQ(ierr);
  ierr = computeTFAFromBasalShearStressUsingPseudoPlastic(vubarSSA, vvbarSSA, 
            taubxComputed, taubyComputed, vtauc, vtillphi); CHKERRQ(ierr);

  ierr = VecDestroy(cIn); CHKERRQ(ierr);
  ierr = VecDestroy(taubxComputed); CHKERRQ(ierr);
  ierr = VecDestroy(taubyComputed); CHKERRQ(ierr);
  return 0;
}


//! Compute the vertical plane shear predicted by the nonsliding SIA, and remove it from a given ice speed to give a basal velocity field.
/*!
This procedure computes a z-independent horizontal ice velocity, conceptually
a basal sliding velocity, from either a specified surface horizontal ice speed 
map or a vertically-averaged horizontal speed map.  The nonsliding
thermomechanically coupled SIA equations are used for this purpose.

We could be given a surface horizontal speed \c csurf or a vertically-integrated 
horizontal speed \c cbar.  The argument \c CisSURF is \c PETSC_TRUE in the 
former case and \c PETSC_FALSE in the latter case.

It is assumed that the input speed gives a velocity pointing down the 
surface gradient.  That is,
  \f[ (u,v) = - c \frac{\nabla h}{|\nabla h|} \f]

If CisSURF is FALSE then this procedure calls only \c surfaceGradientSIA() 
and \c velocitySIAStaggered().  Therefore these get modified:
\c Vec s \c vuvbar[2] and \c vWork2d[0..3] and also 
\c IceModelVec3 s \c Istag3[2] and \c Sigmastag3[2].  The results in
\c vuvbar[2] on the staggered grid are the only quantities used here however.

If CisSURF is TRUE then this procedure calls an additional routine
horizontalVelocitySIARegular().  This modifies 3d \c Vec s \c vu and \c vv.
Also \c vub and \c vvb are set to zero.

The result in \c Vecs \c ub_out and \c vb_out is on the regular grid.  These 2D
\c Vec s must be allocated already.
 */
PetscErrorCode IceModel::removeVerticalPlaneShearRateFromC(
                 const PetscTruth CisSURF, Vec myC, Vec ub_out, Vec vb_out) {
  PetscErrorCode ierr;
  PetscScalar    **csurfSIA;

  // compute the vertically integrated velocity from the nonsliding SIA
  //   (onto the staggered grid)
  ierr = surfaceGradientSIA(); CHKERRQ(ierr); // comm may happen here ...
  // surface gradient temporarily stored in vWork2d[0..3]
  // compute vuvbar[2], which is (ubar,vbar) on staggered grid
  ierr = velocitySIAStaggered(); CHKERRQ(ierr);
  // communicate vuvbar[01] for velocities2DSIAToRegular():
  ierr = DALocalToLocalBegin(grid.da2, vuvbar[0], INSERT_VALUES, vuvbar[0]); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2, vuvbar[0], INSERT_VALUES, vuvbar[0]); CHKERRQ(ierr);
  ierr = DALocalToLocalBegin(grid.da2, vuvbar[1], INSERT_VALUES, vuvbar[1]); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(grid.da2, vuvbar[1], INSERT_VALUES, vuvbar[1]); CHKERRQ(ierr);

  if (CisSURF == PETSC_TRUE) {
    // we need to do some extra work.  velocitySIAStaggered() has been called,
    //   but it only computes the 2d vertically-averaged velocity on the staggered grid
    //   we actually need the surface velocity computed by the nonsliding SIA
    ierr = Istag3[0].beginGhostComm(); CHKERRQ(ierr);
    ierr = Istag3[1].beginGhostComm(); CHKERRQ(ierr);
    ierr = Istag3[0].endGhostComm(); CHKERRQ(ierr);
    ierr = Istag3[1].endGhostComm(); CHKERRQ(ierr);
    ierr = VecSet(vub, 0.0); CHKERRQ(ierr);
    ierr = VecSet(vvb, 0.0); CHKERRQ(ierr);
    ierr = horizontalVelocitySIARegular(); CHKERRQ(ierr);
    ierr = u3.beginGhostComm(); CHKERRQ(ierr);
    ierr = v3.beginGhostComm(); CHKERRQ(ierr);
    ierr = u3.endGhostComm(); CHKERRQ(ierr);
    ierr = v3.endGhostComm(); CHKERRQ(ierr);
    ierr = u3.needAccessToVals(); CHKERRQ(ierr);
    ierr = v3.needAccessToVals(); CHKERRQ(ierr);
    ierr = u3.getSurfaceValuesVec2d(vWork2d[4], vH); CHKERRQ(ierr);
    ierr = v3.getSurfaceValuesVec2d(vWork2d[5], vH); CHKERRQ(ierr);
    ierr = u3.doneAccessToVals(); CHKERRQ(ierr);
    ierr = v3.doneAccessToVals(); CHKERRQ(ierr);
    PetscScalar **usurfSIA, **vsurfSIA;
    ierr = DAVecGetArray(grid.da2, vWork2d[4], &usurfSIA); CHKERRQ(ierr);
    ierr = DAVecGetArray(grid.da2, vWork2d[5], &vsurfSIA); CHKERRQ(ierr);
    ierr = DAVecGetArray(grid.da2, vWork2d[0], &csurfSIA); CHKERRQ(ierr);
    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
        csurfSIA[i][j] = sqrt(PetscSqr(usurfSIA[i][j]) + PetscSqr(vsurfSIA[i][j]));
      }
    }
    ierr = DAVecRestoreArray(grid.da2, vWork2d[4], &usurfSIA); CHKERRQ(ierr);
    ierr = DAVecRestoreArray(grid.da2, vWork2d[5], &vsurfSIA); CHKERRQ(ierr);
    ierr = DAVecRestoreArray(grid.da2, vWork2d[0], &csurfSIA); CHKERRQ(ierr);
  }

  PetscScalar **cc, **ub, **vb, **uvbar[2];
  ierr = DAVecGetArray(grid.da2, vuvbar[0], &uvbar[0]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vuvbar[1], &uvbar[1]); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, ub_out, &ub); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vb_out, &vb); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, myC, &cc); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vWork2d[0], &csurfSIA); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      // note (ubarSIA,vbarSIA) points down the surface gradient; see 
      //   velocities2DSIAToRegular()
      const PetscScalar
          ubarSIA = 0.5 * (uvbar[0][i-1][j] + uvbar[0][i][j]),
          vbarSIA = 0.5 * (uvbar[1][i][j-1] + uvbar[1][i][j]),
          magbarSIA = sqrt(PetscSqr(ubarSIA) + PetscSqr(vbarSIA));
      if (magbarSIA < 1.0 / secpera) {
        // this is the problematic case; presumably the surface gradient is so small
        //   that there is little flow driven by the local driving stress;
        //   so in this case we don't know what direction (ub,vb) should be; 
        //   we just set the sliding velocity to zero
        ub[i][j] = 0.0;
        vb[i][j] = 0.0;
      } else {
        if (CisSURF == PETSC_TRUE) { // cc is ice surface speed
          if (csurfSIA[i][j] >= cc[i][j]) {
            // in this case all the speed in cc can be accounted for by SIA
            //   deformation in vertical planes; therefore no sliding
            ub[i][j] = 0.0;
            vb[i][j] = 0.0;
          } else {
            const PetscScalar factor = (cc[i][j] - csurfSIA[i][j]) / magbarSIA;
            ub[i][j] = factor * ubarSIA;
            vb[i][j] = factor * vbarSIA;
          }
        } else { // cc is vertically-averaged
          if (magbarSIA >= cc[i][j]) {
            // in this case all the speed in cc can be accounted for by SIA
            //   deformation in vertical planes; therefore no sliding
            ub[i][j] = 0.0;
            vb[i][j] = 0.0;
          } else {
            // desired case: remove the speed cc, assuming it is in direction of 
            //   SIA velocity
            const PetscScalar factor = (cc[i][j] - magbarSIA) / magbarSIA;
            ub[i][j] = factor * ubarSIA;
            vb[i][j] = factor * vbarSIA;
          }
        }
      }
    }
  }
  ierr = DAVecRestoreArray(grid.da2, myC, &cc); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, ub_out, &ub); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vb_out, &vb); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vuvbar[0], &uvbar[0]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vuvbar[1], &uvbar[1]); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vWork2d[0], &csurfSIA); CHKERRQ(ierr);
  
  return 0;
}


//! Compute basal shear stress from a given sliding velocity using SSA equations.
/*!
From the input sliding velocity, and using ice thickness, ice surface elevation, and an 
ice temperature field, we use the SSA equations to compute the basal shear
stress.  In particular, the equations we are using are
\latexonly
\def\ddt#1{\ensuremath{\frac{\partial #1}{\partial t}}}
\def\ddx#1{\ensuremath{\frac{\partial #1}{\partial x}}}
\def\ddy#1{\ensuremath{\frac{\partial #1}{\partial y}}}
\begin{align*}
  - 2 \ddx{}\left[\nu H \left(2 \ddx{u} + \ddy{v}\right)\right]
        - \ddy{}\left[\nu H \left(\ddy{u} + \ddx{v}\right)\right]
        + \tau_{(b)x}  &=  - \rho g H \ddx{h}, \\
    - \ddx{}\left[\nu H \left(\ddy{u} + \ddx{v}\right)\right]
      - 2 \ddy{}\left[\nu H \left(\ddx{u} + 2 \ddy{v}\right)\right]
        + \tau_{(b)y}  &=  - \rho g H \ddy{h}, 
\end{align*}
\endlatexonly
Note \f$\nu\f$ is temperature-dependent.

In this case the known field is \f$(u,v)\f$.  We are solving for 
\f$\tau_b = (\tau_{(b)x},\tau_{(b)y})\f$.  The equations are solved everywhere 
that there is positive thickness; at other points this procedure returns 
\f$\tau_b = 0\f$.

This procedure must save temporary versions of everything modified by 
assembleSSAMatrix() and assembleSSARhs() and
computeEffectiveViscosity() because it calls those routines.

The result in \c Vecs \c taubx_out and \c taubx_out are on the regular grid.
They must be allocated before calling this routine.
 */
PetscErrorCode IceModel::computeBasalShearFromSSA(
                 Vec myu, Vec myv, Vec taubx_out, Vec tauby_out) {
  PetscErrorCode ierr;

  // effective viscosity for myu, myv
  const PetscTruth leaveNuAloneSSA_save = leaveNuAloneSSA, 
                   useConstantNuForSSA_save = useConstantNuForSSA;
  leaveNuAloneSSA = PETSC_FALSE;
  useConstantNuForSSA = PETSC_FALSE;
  Vec myvNu[2] = {vWork2d[0], vWork2d[1]}; // already allocated space
  Vec vubarSSAold = vWork2d[2], vvbarSSAold = vWork2d[3];
  ierr = VecCopy(vubarSSA, vubarSSAold); CHKERRQ(ierr);
  ierr = VecCopy(vvbarSSA, vvbarSSAold); CHKERRQ(ierr);
  ierr = VecCopy(myu, vubarSSA); CHKERRQ(ierr);
  ierr = VecCopy(myv, vvbarSSA); CHKERRQ(ierr);
  // eps=0.0 in bdd-below regularization; Schoof type regularization does occur;
  ierr = computeEffectiveViscosity(myvNu, 0.0); CHKERRQ(ierr);
  ierr = VecCopy(vubarSSAold, vubarSSA); CHKERRQ(ierr);
  ierr = VecCopy(vvbarSSAold, vvbarSSA); CHKERRQ(ierr);
  leaveNuAloneSSA = leaveNuAloneSSA_save;
  useConstantNuForSSA = useConstantNuForSSA_save;

  // allocate Mat and Vecs; compare linear system in velocitySSA()
  Mat A;
  Vec x, result, rhs;
  ierr = MatDuplicate(SSAStiffnessMatrix, MAT_COPY_VALUES, &A); CHKERRQ(ierr);
  ierr = VecDuplicate(SSAX, &result); CHKERRQ(ierr);
  ierr = VecDuplicate(SSAX, &rhs); CHKERRQ(ierr);
  ierr = VecDuplicate(SSAX, &x); CHKERRQ(ierr);

  // Build approximation of SSA equations (A x = rhs)
  //   Generally we want basal shear as part of matrix:
  //                L[U] + beta U = (driving),
  //   so  A = L + beta I. In this case we don't want basal shear terms 
  //   as part of the matrix:
  //                tau_b := L[U] - (driving),
  //   so  A = L.
  ierr = assembleSSAMatrix(false, myvNu, A); CHKERRQ(ierr);
 
  // Note rhs contains driving terms  - \rho g H \grad h
  ierr = assembleSSARhs(false, rhs); CHKERRQ(ierr);

  // Set x = [u v]'  (i.e. interleaved).
  PetscScalar **u, **v;
  const PetscInt  M = 2 * grid.My;
  ierr = VecSet(x, 0.0); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, myu, &u); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, myv, &v); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      ierr = VecSetValue(x, i*M + 2*j, u[i][j], INSERT_VALUES); CHKERRQ(ierr);
      ierr = VecSetValue(x, i*M + 2*j+1, v[i][j], INSERT_VALUES); CHKERRQ(ierr);
    }
  }
  ierr = DAVecRestoreArray(grid.da2, myu, &u); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, myv, &v); CHKERRQ(ierr);

  // communicate!
  ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(x); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(x); CHKERRQ(ierr);

  // With new A and x, compute  result = Ax - rhs = L[u] - (driving)
  ierr = MatMult(A,x,result); CHKERRQ(ierr);
  ierr = VecAXPY(result,-1.0,rhs); CHKERRQ(ierr);  // note result = 0 if MASK_SHEET
  
  // transfer result to taub output vectors
  PetscScalar     **tbx, **tby, *res;
  Vec             resultLoc = SSAXLocal;
  ierr = VecScatterBegin(SSAScatterGlobalToLocal, result, resultLoc, 
           INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecScatterEnd(SSAScatterGlobalToLocal, result, resultLoc, 
           INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGetArray(resultLoc, &res); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, taubx_out, &tbx); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, tauby_out, &tby); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      tbx[i][j] = res[i*M + 2*j];
      tby[i][j] = res[i*M + 2*j+1];
    }
  }
  ierr = DAVecGetArray(grid.da2, taubx_out, &tbx); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, tauby_out, &tby); CHKERRQ(ierr);
  ierr = VecRestoreArray(resultLoc, &res); CHKERRQ(ierr);

  // de-allocate
  ierr = MatDestroy(A); CHKERRQ(ierr);
  ierr = VecDestroy(x); CHKERRQ(ierr);
  ierr = VecDestroy(result); CHKERRQ(ierr);
  ierr = VecDestroy(rhs); CHKERRQ(ierr);

  return 0;
}


//! Compute till friction angle in degrees.
/*!
The more complete brief description might be: "Compute till friction angle in 
degrees from saved basal shear stress map and saved basal velocity 
(both vector-valued) using pseudo-plastic till constitutive relation 
and a pore water pressure model."

Here is the context for this procedure:

This procedure should not be called in the purely plastic case, that is, 
if <tt>PlasticBasalType::pseudo_plastic</tt> is \c FALSE.

If <tt>PlasticBasalType::pseudo_plastic</tt> is \c TRUE then
PlasticBasalType::drag() will compute the basal shear stress
as
    \f[ (\tau_{(b)x},\tau_{(b)y})
           = - \frac{\tau_c}{|U|^{1-q} |U_{\mathtt{thresh}}|^q} U\f]
where \f$U=(u,v)\f$, \f$q=\f$ <tt>pseudo_q</tt>, and
\f$U_{\mathtt{thresh}}=\f$ <tt>pseudo_u_threshold</tt>.
Typical values for the latter two constants are \f$q=0.25\f$ 
and 100 m/a for \f$U_{\mathtt{thresh}}\f$, but \f$q\f$ could even be
set to 1 for a linear till model, in which case
\f$\beta = \tau_c/|U_{\mathtt{thresh}}|\f$.

On the other hand, updateYieldStressFromHmelt() computes
the till yield stress \f$\tau_c\f$, a scalar, by 
    \f[ \tau_c = c_0 + \tan\left(\frac{\pi}{180} \phi\right)\, N \f]
where \f$\phi\f$ is the map of till friction angle, in degrees, stored in
\c vtillphi.  The till cohesion \f$c_0\f$ (= \c plastic_till_c_0) is usually zero.
Here \f$N\f$ is the effective pressure on the till, namely \f$N = \rho g H - p_w\f$.
The porewater pressure \f$p_w\f$ is estimated in terms of the stored basal water
\c vHmelt; see getEffectivePressureOnTill().

In this context, the current procedure uses a sliding velocity 
and a basal shear stress (sign choice so that it is stress applied to the ice)
to compute the yield stress and then the till friction angle.  
There are some cases to check:
  - If the velocity is below 1 m/a in magnitude then \f$\tau_c\f$ 
    is set equal to \f$10 |\tau_{(b)}|\f$.
  - It is not assumed that the vectors \f$U\f$ and \f$\tau_{(b)}\f$
    at each grid point are pointing in the same direction 
    on input to this routine.  If the velocity is larger in magnitude than
    1 m/a then the angle between \f$U\f$ and \f$\tau_{(b)}\f$ is checked; 
    if that angle exceeds 45 degrees then \f$\tau_c\f$ is set equal to 
    \f$10 |\tau_{(b)}|\f$.
  - If the last two cases did not apply, compute the yield stress through magnitudes:
    \f[ \tau_c = \frac{|\tau_{(b)}|\,|U_{\mathtt{thresh}}|^q}{|U|^q} \f]
    by PlasticBasalType::taucFromMagnitudes().
  - Compute the till friction angle \e from the till yield stress
    \f[ \phi = \frac{180}{\pi} \arctan\left(\frac{\tau_c - c_0}{N}\right) \f]

All of the above is done only if the ice is grounded and there is positive
ice thickness.  If a point is marked \c MASK_FLOATING or \c MASK_FLOATING_OCEAN0 
then the yield stress is set to zero and the friction angle (quite arbitrarily) to 
30 degrees.  If a point is grounded but the ice thickness is zero the the
yield stress is set to a high value of 1000 kPa = 10 bar and again the friction angle 
is set to 30 degrees.

Values of both \f$\tau_c\f$ and \f$\phi\f$ are returned by this routine
but it is envisioned that the value of \f$\tau_c\f$ may evolve in a run
while \f$\phi\f$ is treated as time-independent.  No communication occurs
in this routine.

This procedure reads three fields from the current model state, namely ice 
thickness \c vH, effective thickness of basal meltwater \c vHmelt, and the 
mask \c vMask.
 */
PetscErrorCode IceModel::computeTFAFromBasalShearStressUsingPseudoPlastic(
                 const Vec myu, const Vec myv, const Vec mytaubx, const Vec mytauby, 
                 Vec tauc_out, Vec tfa_out) {
  PetscErrorCode ierr;
  
  const PetscScalar slowOrWrongDirFactor = 10.0,
                    sufficientSpeed = 1.0 / secpera,
                    sameDirMaxAngle = pi/4.0;
  PetscScalar **tauc, **phi, **u, **v, **taubx, **tauby, 
              **Hmelt, **H, **mask;
  ierr = DAVecGetArray(grid.da2, myu, &u); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, myv, &v); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, mytaubx, &taubx); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, mytauby, &tauby); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, tauc_out, &tauc); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, tfa_out, &phi); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vHmelt, &Hmelt); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vMask, &mask); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (modMask(mask[i][j]) == MASK_FLOATING) {
        tauc[i][j] = 0.0;
        phi[i][j] = 30.0;
      } else if (H[i][j] <= 0.0) {
        tauc[i][j] = 1000.0e3;  // large; 1000 kPa = 10 bar if no ice
        phi[i][j] = 30.0;
      } else { // grounded and ice is present
        // first get tauc
        const PetscScalar
            speed = sqrt(PetscSqr(u[i][j]) + PetscSqr(v[i][j])),
            taubmag = sqrt(PetscSqr(taubx[i][j]) + PetscSqr(tauby[i][j]));
        if (speed < sufficientSpeed) {
          tauc[i][j] = slowOrWrongDirFactor * taubmag;
        } else {
          const PetscScalar 
              ubDOTtaub = u[i][j] * taubx[i][j] + tauby[i][j] * v[i][j],
              angle = acos(ubDOTtaub / (speed * taubmag));
          if (PetscAbs(angle) > sameDirMaxAngle) {
            tauc[i][j] = slowOrWrongDirFactor * taubmag;
          } else {
            // use the formula which inverts PlasticBasalType::drag()
            tauc[i][j] = basal->taucFromMagnitudes(taubmag, speed);
          }
        }
        // now get till friction angle
        if (tauc[i][j] > plastic_till_c_0) {
          const PetscScalar N = getEffectivePressureOnTill(H[i][j], Hmelt[i][j]);
          phi[i][j] = (180.0 / pi) * atan( (tauc[i][j] - plastic_till_c_0) / N );
        } else {
          phi[i][j] = 0.0;
        }
      }
    }
  }
  ierr = DAVecRestoreArray(grid.da2, myu, &u); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, myv, &v); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, mytaubx, &taubx); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, mytauby, &tauby); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, tauc_out, &tauc); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, tfa_out, &phi); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vHmelt, &Hmelt); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vMask, &mask); CHKERRQ(ierr);

  return 0;
}

