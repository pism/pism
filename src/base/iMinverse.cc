// Copyright (C) 2008 Ed Bueler and Constantine Khroulev
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
  
  // check whether either option is set; get filename if so
  PetscTruth cbarTillSet, csurfTillSet;
  char filename[PETSC_MAX_PATH_LEN];
  ierr = PetscOptionsGetString(PETSC_NULL, "-cbar_to_till", filename, 
                               PETSC_MAX_PATH_LEN, &cbarTillSet); CHKERRQ(ierr);
  ierr = PetscOptionsGetString(PETSC_NULL, "-csurf_to_till", filename, 
                               PETSC_MAX_PATH_LEN, &csurfTillSet); CHKERRQ(ierr);
  if ((cbarTillSet == PETSC_TRUE) && (csurfTillSet == PETSC_TRUE)) {
    SETERRQ(1,"both -cbar_to_till and -csurf_to_till set; NOT ALLOWED\n");
  }
  if ((cbarTillSet == PETSC_FALSE) && (csurfTillSet == PETSC_FALSE)) {
    return 0;  // if neither option set
  }

  // check whether file exists, and whether it contains cbar/csurf variable
  bool file_exists = false;
  NCTool nc(&grid);

  ierr = nc.open_for_reading(filename, file_exists); CHKERRQ(ierr);

  if (!file_exists) {
    SETERRQ(2,"NetCDF file with velocities (either cbar or csurf) not found.\n");
  }
  bool cbar_exists = false, csurf_exists = false;
  ierr = nc.find_variable("cbar",  NULL, NULL, cbar_exists); CHKERRQ(ierr);
  ierr = nc.find_variable("csurf", NULL, NULL, csurf_exists); CHKERRQ(ierr);

  if ((cbarTillSet == PETSC_TRUE) && (!cbar_exists)) {
    SETERRQ1(3,"-cbar_to_till set but cbar not found in %s\n", filename);
  }
  if ((csurfTillSet == PETSC_TRUE) && (!csurf_exists)) {
    SETERRQ1(4,"-csurf_to_till set but csurf not found in %s\n", filename);
  }

  char varName[10];
  strcpy(varName, (cbarTillSet == PETSC_TRUE) ? "cbar" : "csurf");
  ierr = verbPrintf(2, grid.com, 
     "apply inverse model to %s in file %s\n"
     "  reading and storing %s ...\n",
     varName, filename, varName); CHKERRQ(ierr);

  // create "local interpolation context" from dimensions, 
  //   limits, and lengths; compare code in bootstrapFromFile_netCDF()
  size_t dim[5];
  double bdy[7];
  ierr = nc.get_dims_limits_lengths_2d(dim, bdy); CHKERRQ(ierr);
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
  LocalInterpCtx lic(dim, bdy, z_bif, zb_bif, grid);
  delete z_bif;
  delete zb_bif;
  //DEBUG:  ierr = lic.printGrid(grid.com); CHKERRQ(ierr);

  // space for read-in var; read it (regridding as you go)
  IceModelVec2 cIn;
  ierr = cIn.create(grid, varName, true); CHKERRQ(ierr);

  ierr = cIn.regrid(filename, lic, true); CHKERRQ(ierr);

  ierr = nc.close(); CHKERRQ(ierr);

  // space for computed basal shear
  IceModelVec2 taubxComputed, taubyComputed;
  ierr = taubxComputed.create(grid, "taubx", true);
  ierr = taubyComputed.create(grid, "tauby", true);

  // do inverse model; result is tauc and tillphi fields
  ierr = verbPrintf(2, grid.com, 
     "  removing SIA vertical shear ...\n"); CHKERRQ(ierr);
  ierr = removeVerticalPlaneShearRateFromC(csurfTillSet, cIn, vubarSSA, vvbarSSA);
            CHKERRQ(ierr);
  ierr = verbPrintf(2, grid.com, 
     "  computing membrane stresses and basal shear stress from SSA ...\n"); CHKERRQ(ierr);
  ierr = computeBasalShearFromSSA(vubarSSA, vvbarSSA, taubxComputed, taubyComputed);
            CHKERRQ(ierr);
  ierr = verbPrintf(2, grid.com, 
     "  computing till friction angle using (pseudo-)plastic model ...\n"); 
     CHKERRQ(ierr);
  ierr = computeTFAFromBasalShearStressUsingPseudoPlastic(vubarSSA, vvbarSSA, 
            taubxComputed, taubyComputed, vtauc, vtillphi); CHKERRQ(ierr);

  // clean up
  ierr =           cIn.destroy(); CHKERRQ(ierr);
  ierr = taubxComputed.destroy(); CHKERRQ(ierr);
  ierr = taubyComputed.destroy(); CHKERRQ(ierr);
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
                 const PetscTruth CisSURF,
		 IceModelVec2 myC, IceModelVec2 ub_out, IceModelVec2 vb_out) {
  PetscErrorCode ierr;
  PetscScalar    **csurfSIA;

  // compute the vertically integrated velocity from the nonsliding SIA
  //   (onto the staggered grid)
  ierr = surfaceGradientSIA(); CHKERRQ(ierr); // comm may happen here ...
  // surface gradient temporarily stored in vWork2d[0..3]
  // compute vuvbar[2], which is (ubar,vbar) on staggered grid
  ierr = velocitySIAStaggered(); CHKERRQ(ierr);
  // communicate vuvbar[01] for velocities2DSIAToRegular():
  ierr = vuvbar[0].beginGhostComm(); CHKERRQ(ierr);
  ierr = vuvbar[1].beginGhostComm(); CHKERRQ(ierr);
  ierr = vuvbar[0].endGhostComm(); CHKERRQ(ierr);
  ierr = vuvbar[1].endGhostComm(); CHKERRQ(ierr);

  if (CisSURF == PETSC_TRUE) {
    // we need to do some extra work.  velocitySIAStaggered() has been called,
    //   but it only computes the 2d vertically-averaged velocity on the staggered grid
    //   we actually need the surface velocity computed by the nonsliding SIA
    ierr = Istag3[0].beginGhostComm(); CHKERRQ(ierr);
    ierr = Istag3[1].beginGhostComm(); CHKERRQ(ierr);
    ierr = Istag3[0].endGhostComm(); CHKERRQ(ierr);
    ierr = Istag3[1].endGhostComm(); CHKERRQ(ierr);
    ierr = vub.set(0.0); CHKERRQ(ierr);
    ierr = vvb.set(0.0); CHKERRQ(ierr);
    ierr = horizontalVelocitySIARegular(); CHKERRQ(ierr);
    ierr = u3.beginGhostComm(); CHKERRQ(ierr);
    ierr = v3.beginGhostComm(); CHKERRQ(ierr);
    ierr = u3.endGhostComm(); CHKERRQ(ierr);
    ierr = v3.endGhostComm(); CHKERRQ(ierr);
    ierr = u3.begin_access(); CHKERRQ(ierr);
    ierr = v3.begin_access(); CHKERRQ(ierr);
    ierr = u3.getSurfaceValues(vWork2d[4], vH); CHKERRQ(ierr);
    ierr = v3.getSurfaceValues(vWork2d[5], vH); CHKERRQ(ierr);
    ierr = u3.end_access(); CHKERRQ(ierr);
    ierr = v3.end_access(); CHKERRQ(ierr);
    PetscScalar **usurfSIA, **vsurfSIA;
    ierr = vWork2d[4].get_array(usurfSIA); CHKERRQ(ierr);
    ierr = vWork2d[5].get_array(vsurfSIA); CHKERRQ(ierr);
    ierr = vWork2d[0].get_array(csurfSIA); CHKERRQ(ierr);
    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
        csurfSIA[i][j] = sqrt(PetscSqr(usurfSIA[i][j]) + PetscSqr(vsurfSIA[i][j]));
      }
    }
    ierr = vWork2d[4].end_access(); CHKERRQ(ierr);
    ierr = vWork2d[5].end_access(); CHKERRQ(ierr);
    ierr = vWork2d[0].end_access(); CHKERRQ(ierr);
  }

  PetscScalar **cc, **ub, **vb, **uvbar[2];
  ierr =  vuvbar[0].get_array(uvbar[0]); CHKERRQ(ierr);
  ierr =  vuvbar[1].get_array(uvbar[1]); CHKERRQ(ierr);
  ierr =     ub_out.get_array(ub);       CHKERRQ(ierr);
  ierr =     vb_out.get_array(vb);       CHKERRQ(ierr);
  ierr =        myC.get_array(cc);       CHKERRQ(ierr);
  ierr = vWork2d[0].get_array(csurfSIA); CHKERRQ(ierr);
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
  ierr =        myC.end_access(); CHKERRQ(ierr);
  ierr =     ub_out.end_access(); CHKERRQ(ierr);
  ierr =     vb_out.end_access(); CHKERRQ(ierr);
  ierr =  vuvbar[0].end_access(); CHKERRQ(ierr);
  ierr =  vuvbar[1].end_access(); CHKERRQ(ierr);
  ierr = vWork2d[0].end_access(); CHKERRQ(ierr);
  
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
PetscErrorCode IceModel::computeBasalShearFromSSA(IceModelVec2 myu, IceModelVec2 myv,
						  IceModelVec2 taubx_out, IceModelVec2 tauby_out) {
  PetscErrorCode ierr;

  // FIXME: there is a change from "Nu" to "NuH" notation; follow it more completely:
  // effective viscosity for myu, myv
  const PetscTruth leaveNuAloneSSA_save = leaveNuAloneSSA, 
                   useConstantNuHForSSA_save = useConstantNuHForSSA;
  leaveNuAloneSSA = PETSC_FALSE;
  useConstantNuHForSSA = PETSC_FALSE;
  IceModelVec2 myvNu[2] = {vWork2d[0], vWork2d[1]}; // already allocated space
  IceModelVec2 vubarSSAold = vWork2d[2], vvbarSSAold = vWork2d[3];
  ierr = vubarSSA.copy_to(vubarSSAold); CHKERRQ(ierr);
  ierr = vvbarSSA.copy_to(vvbarSSAold); CHKERRQ(ierr);
  ierr = myu.copy_to(vubarSSA); CHKERRQ(ierr);
  ierr = myv.copy_to(vvbarSSA); CHKERRQ(ierr);
  // eps=0.0 in bdd-below regularization; Schoof type regularization does occur;
  ierr = computeEffectiveViscosity(myvNu, 0.0); CHKERRQ(ierr);
  ierr = vubarSSAold.copy_to(vubarSSA); CHKERRQ(ierr);
  ierr = vvbarSSAold.copy_to(vvbarSSA); CHKERRQ(ierr);
  leaveNuAloneSSA = leaveNuAloneSSA_save;
  useConstantNuHForSSA = useConstantNuHForSSA_save;

  // allocate Mat and Vecs; compare linear system in velocitySSA()
  Mat A;
  Vec x, result, rhs;
//  ierr = MatDuplicate(SSAStiffnessMatrix, MAT_COPY_VALUES, &A); CHKERRQ(ierr);
  const PetscInt M = 2 * grid.Mx * grid.My;
  ierr = MatCreateMPIAIJ(grid.com, PETSC_DECIDE, PETSC_DECIDE, M, M,
                         13, PETSC_NULL, 13, PETSC_NULL,
                         &A); CHKERRQ(ierr);
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
  const PetscInt  twoMy = 2 * grid.My;
  ierr = VecSet(x, 0.0); CHKERRQ(ierr);
  ierr = myu.get_array(u); CHKERRQ(ierr);
  ierr = myv.get_array(v); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      ierr = VecSetValue(x, i*twoMy + 2*j, u[i][j], INSERT_VALUES); CHKERRQ(ierr);
      ierr = VecSetValue(x, i*twoMy + 2*j+1, v[i][j], INSERT_VALUES); CHKERRQ(ierr);
    }
  }
  ierr = myu.end_access(); CHKERRQ(ierr);
  ierr = myv.end_access(); CHKERRQ(ierr);

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
  ierr = taubx_out.get_array(tbx); CHKERRQ(ierr);
  ierr = tauby_out.get_array(tby); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      tbx[i][j] = res[i*twoMy + 2*j];
      tby[i][j] = res[i*twoMy + 2*j+1];
    }
  }
  ierr = taubx_out.end_access(); CHKERRQ(ierr);
  ierr = tauby_out.end_access(); CHKERRQ(ierr);
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
                 IceModelVec2 myu, IceModelVec2 myv,
		 IceModelVec2 mytaubx, IceModelVec2 mytauby, 
                 IceModelVec2 tauc_out, IceModelVec2 tfa_out) {
  PetscErrorCode ierr;
  
  // FIXME: this should really ask the PlasticBasalType instance basal-> whether
  //   it is pseudo or not
  if (doPseudoPlasticTill == PETSC_FALSE) {
    ierr = verbPrintf(1, grid.com, 
       "WARNING: computeTFAFromBasalShearStress() should only be called with q > 0.0\n");
       CHKERRQ(ierr);
  }
  
  const PetscScalar slowOrWrongDirFactor = 10.0,
                    sufficientSpeed = 1.0 / secpera,
                    sameDirMaxAngle = pi/4.0;
  PetscScalar **tauc, **phi, **u, **v, **taubx, **tauby, 
              **Hmelt, **H, **mask;
  ierr =      myu.get_array(u);     CHKERRQ(ierr);
  ierr =      myv.get_array(v);     CHKERRQ(ierr);
  ierr =  mytaubx.get_array(taubx); CHKERRQ(ierr);
  ierr =  mytauby.get_array(tauby); CHKERRQ(ierr);
  ierr = tauc_out.get_array(tauc);  CHKERRQ(ierr);
  ierr =  tfa_out.get_array(phi);   CHKERRQ(ierr);
  ierr =   vHmelt.get_array(Hmelt); CHKERRQ(ierr);
  ierr =       vH.get_array(H);     CHKERRQ(ierr);
  ierr =    vMask.get_array(mask);  CHKERRQ(ierr);
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
  ierr =      myu.end_access(); CHKERRQ(ierr);
  ierr =      myv.end_access(); CHKERRQ(ierr);
  ierr =  mytaubx.end_access(); CHKERRQ(ierr);
  ierr =  mytauby.end_access(); CHKERRQ(ierr);
  ierr = tauc_out.end_access(); CHKERRQ(ierr);
  ierr =  tfa_out.end_access(); CHKERRQ(ierr);
  ierr =   vHmelt.end_access(); CHKERRQ(ierr);
  ierr =       vH.end_access(); CHKERRQ(ierr);
  ierr =    vMask.end_access(); CHKERRQ(ierr);

  return 0;
}

