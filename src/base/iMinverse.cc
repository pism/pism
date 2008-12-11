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


//! Invert given ice surface velocities to find till friction angle using a pseudo-plastic model.
/*!
Reads user option <tt>-surf_vel_to_tfa foo.nc</tt>.  If this option is not found,
this method returns without doing anything.

If the option is found, this method expects to find variables
\c us and \c vs, namely x and y components of surface velocity,
in NetCDF file \c foo.nc.  If found, reads them.

Als0 looks for option <tt>-fofv bar.nc</tt>.  If this optional option is given,
then expects variable \c fofv in bar.nc.  This procedure uses \c fofv in 
removing the vertical plane shear rate (i.e. the SIA computed velocity field).

This method calls the three major steps in turn:
-# removeVerticalPlaneShearRate() to compute the basal sliding velocity,
-# computeBasalShearFromSSA() to compute the basal shear stress, and finally
-# computeTFAFromBasalShearStressUsingPseudoPlastic().

This method must be called after initialization is essentially completed.
Uses geometry, temperature field, and map of effective thickness of basal till
water in the inversion.

The result is a field of values for till yield stress and the till friction angle in a 
pseudo-plastic till model.  This model includes purely-plastic and linearly-viscous cases.

A forward time-stepping model could then be based on the resulting till friction angle field.
 */
PetscErrorCode IceModel::invertSurfaceVelocities(const PetscTruth writeWhenDone) {
  PetscErrorCode ierr;

  IceModelVec2 usIn, vsIn, taubxComputed, taubyComputed, fofv;
  
  PetscTruth svTOtfaSet, fofvSet;
  char filename[PETSC_MAX_PATH_LEN], fofvname[PETSC_MAX_PATH_LEN];
  ierr = PetscOptionsGetString(PETSC_NULL, "-surf_vel_to_tfa", filename, 
                               PETSC_MAX_PATH_LEN, &svTOtfaSet); CHKERRQ(ierr);
  if (svTOtfaSet == PETSC_FALSE) {
    return 0;  // leave if you are not wanted ...
  }

  ierr = PetscOptionsGetString(PETSC_NULL, "-fofv", fofvname, 
                               PETSC_MAX_PATH_LEN, &fofvSet); CHKERRQ(ierr);

  // create inverse model state variables with metadata
  ierr = usIn.create(grid, "us", true); CHKERRQ(ierr);
  ierr = usIn.set_attrs(
     NULL, "x component of ice surface velocity","m s-1", NULL); CHKERRQ(ierr);
  ierr = vsIn.create(grid, "vs", true); CHKERRQ(ierr);
  ierr = vsIn.set_attrs(
     NULL, "y component of ice surface velocity","m s-1", NULL); CHKERRQ(ierr);
  ierr = taubxComputed.create(grid, "taubxOUT", false);  // global
  ierr = taubxComputed.set_attrs(
     NULL, "x component of basal shear stress", "Pa", NULL); CHKERRQ(ierr);
  ierr = taubyComputed.create(grid, "taubyOUT", false);  // global
  ierr = taubyComputed.set_attrs(
     NULL, "y component of basal shear stress", "Pa", NULL); CHKERRQ(ierr);
  ierr = fofv.create(grid, "fofv", false);  // global
  ierr = fofv.set_attrs(
     NULL, "fraction of final velocity which comes from SIA in hybrid model", 
     "1", NULL); CHKERRQ(ierr);

  // read in
  ierr = verbPrintf(2, grid.com, 
     "option -surf_vel_to_tfa seen ... doing ad hoc inverse model;\n"
     "  reading us, vs from file %s ...\n", filename); CHKERRQ(ierr);
  ierr = getIMV2FromFile(filename, usIn); CHKERRQ(ierr);
  ierr = getIMV2FromFile(filename, vsIn); CHKERRQ(ierr);
  if (fofvSet == PETSC_TRUE) {
    ierr = verbPrintf(2, grid.com, 
       "option -fofv seen ... reading f(|v|) from file %s ...\n",
       fofvname); CHKERRQ(ierr);
    ierr = getIMV2FromFile(fofvname, fofv); CHKERRQ(ierr);
  } else {
    // compute an f(|v|) field from current ubarSSA, vbarSSA
    ierr = verbPrintf(2, grid.com, 
       "option -fofv NOT seen ...\n"
       "  (computing f(|v|) from CURRENT values (from -if file?) of vubarSSA,vvbarSSA ...)\n");
       CHKERRQ(ierr);
    ierr = computeFofV(vubarSSA, vvbarSSA, fofv); CHKERRQ(ierr);
  }

  // do three steps of inverse model; result is tauc and tillphi fields
  ierr = verbPrintf(2, grid.com, 
           "  removing SIA vertical shear ...\n"); CHKERRQ(ierr);
  ierr = removeVerticalPlaneShearRate(
           usIn, vsIn, fofv, vubarSSA, vvbarSSA);    CHKERRQ(ierr);
  ierr = verbPrintf(2, grid.com, 
           "  computing membrane stresses and basal shear stress from SSA;\n"
           "  (applying SSA differential operator to us,vs) ...\n");
           CHKERRQ(ierr);
  ierr = computeBasalShearFromSSA(
           vubarSSA, vvbarSSA, taubxComputed, taubyComputed);   CHKERRQ(ierr);
  ierr = verbPrintf(2, grid.com, 
           "  computing till friction angle using (pseudo-)plastic model ...\n"); 
           CHKERRQ(ierr);
  ierr = computeTFAFromBasalShearStressUsingPseudoPlastic(
           vubarSSA, vvbarSSA, taubxComputed, taubyComputed, vtauc, vtillphi);
           CHKERRQ(ierr);

  // write out stored inverse info for user's edification; mostly debug
  if (writeWhenDone == PETSC_TRUE) {
    char filename[PETSC_MAX_PATH_LEN];
    PetscTruth o_set;
    ierr = PetscOptionsGetString(PETSC_NULL, "-o", filename, 
                                 PETSC_MAX_PATH_LEN, &o_set); CHKERRQ(ierr);
    if (!o_set) {
      ierr = PetscPrintf(grid.com,
		"PISM ERROR: inverse model can only write to file if -o option given\n");
		CHKERRQ(ierr);
      PetscEnd();
    }
    ierr = usIn.write(filename, NC_DOUBLE); CHKERRQ(ierr);
    ierr = vsIn.write(filename, NC_DOUBLE); CHKERRQ(ierr);
    ierr = taubxComputed.write(filename, NC_DOUBLE); CHKERRQ(ierr);
    ierr = taubyComputed.write(filename, NC_DOUBLE); CHKERRQ(ierr);
    ierr = fofv.write(filename, NC_DOUBLE); CHKERRQ(ierr);
  }
  
  // clean up
  ierr =          usIn.destroy(); CHKERRQ(ierr);
  ierr =          vsIn.destroy(); CHKERRQ(ierr);
  ierr = taubxComputed.destroy(); CHKERRQ(ierr);
  ierr = taubyComputed.destroy(); CHKERRQ(ierr);
  ierr =          fofv.destroy(); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceModel::getIMV2FromFile(const char *fname, IceModelVec2 v) {
  PetscErrorCode ierr;

  // check whether file exists
  bool file_exists = false;
  NCTool nc(&grid);
  ierr = nc.open_for_reading(fname, file_exists); CHKERRQ(ierr);
  if (!file_exists) {
    ierr = PetscPrintf(grid.com,
	"PISM ERROR: Can't open file '%s'.\n", fname); CHKERRQ(ierr);
    PetscEnd();
  }

  // create "local interpolation context" from dimensions, 
  //   limits, and lengths; compare code in bootstrapFromFile()
  grid_info g;
  ierr = nc.get_grid_info_2d(g); CHKERRQ(ierr);
  ierr = nc.close(); CHKERRQ(ierr);
  LocalInterpCtx lic(g, NULL, NULL, grid); // 2D only
  //DEBUG:  ierr = lic.printGrid(grid.com); CHKERRQ(ierr);

  // read in (by regridding); checks for existence of file and variable in file
  ierr = v.regrid(fname, lic, true); CHKERRQ(ierr);  // at this point it *is* critical

  ierr = nc.close(); CHKERRQ(ierr);
  return 0;
}


// compare broadcastSSAVelocity(); this should be put in a better spot
//   and/or used in broadcastSSAVelocity(); 
// this variable is "f(|\bv|)" in preprint Bueler & Brown 2008;
//   see equations (21) and (22) in that preprint; http://arxiv.org/abs/0810.3449
PetscErrorCode IceModel::computeFofV(IceModelVec2 uSSA, IceModelVec2 vSSA, 
                                    IceModelVec2 &fofv) {
  PetscErrorCode ierr;
  PetscScalar **ubarssa, **vbarssa, **f;

  const PetscScalar 
     inC_fofv = 1.0e-4 * PetscSqr(secpera),
     outC_fofv = 2.0 / pi;

  ierr = uSSA.get_array(ubarssa); CHKERRQ(ierr);
  ierr = vSSA.get_array(vbarssa); CHKERRQ(ierr);
  ierr = fofv.get_array(f); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      const PetscScalar c2 = PetscSqr(ubarssa[i][j]) + PetscSqr(vbarssa[i][j]);
      f[i][j] = 1.0 - outC_fofv * atan(inC_fofv * c2);
    }
  }
  ierr = fofv.end_access(); CHKERRQ(ierr);
  ierr = uSSA.end_access(); CHKERRQ(ierr);
  ierr = vSSA.end_access(); CHKERRQ(ierr);
  return 0;
}


//! Compute the vertical plane shear predicted by nonsliding SIA.  Remove it from ice surface velocity to give ice basal velocity.
/*!
This procedure computes a z-independent horizontal ice velocity, conceptually
a basal sliding velocity, from a specified surface horizontal ice velocity
(map).  The nonsliding thermomechanically coupled SIA equations are used for this purpose.

In this implementation we extract the SSA velocity by assuming \f$f(|\mathbf{v}|)\f$ is
known in formula (21) of \lo\cite{BBssasliding}\elo.

This procedure calls surfaceGradientSIA(), velocitySIAStaggered(), and
horizontalVelocitySIARegular().  Therefore these IceModel members get modified:

- IceModelVec2 vuvbar[2]
- IceModelVec2 vWork2d[0..5]
- IceModelVec2 vub  (set to zero)
- IceModelVec2 vvb  (set to zero)
- IceModelVec3 Istag3[2] 
- IceModelVec3 Sigmastag3[2].
- IceModelVec3 vu
- IceModelVec3 vv
 
Only \c vu and \c vv are directly used here, however; they are evaluated for 
their surface velocity.

The results in ub_out and vb_out are on the regular grid.  These IceModelVec2
must be allocated (created) before calling this procedure.
 */
PetscErrorCode IceModel::removeVerticalPlaneShearRate(
		 IceModelVec2 us_in, IceModelVec2 vs_in, IceModelVec2 vfofv,
		 IceModelVec2 &ub_out, IceModelVec2 &vb_out) {
  PetscErrorCode ierr;

  // surface values computed by SIA
  IceModelVec2   vusSIA = vWork2d[4],  vvsSIA = vWork2d[5];

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

  // we have the 2d vertically-averaged velocity on the staggered grid, but
  //   we actually need the surface velocity computed by the nonsliding SIA
  //   on the regular grid, so now we call horizontalVelocitySIARegular();
  // re communication pattern: compare code in velocity()
  ierr = Istag3[0].beginGhostComm(); CHKERRQ(ierr);
  ierr = Istag3[1].beginGhostComm(); CHKERRQ(ierr);
  ierr = Istag3[0].endGhostComm(); CHKERRQ(ierr);
  ierr = Istag3[1].endGhostComm(); CHKERRQ(ierr);
  ierr = vub.set(0.0); CHKERRQ(ierr);  // no sliding SIA
  ierr = vvb.set(0.0); CHKERRQ(ierr);
  ierr = horizontalVelocitySIARegular(); CHKERRQ(ierr);
  ierr = u3.beginGhostComm(); CHKERRQ(ierr);
  ierr = v3.beginGhostComm(); CHKERRQ(ierr);
  ierr = u3.endGhostComm(); CHKERRQ(ierr);
  ierr = v3.endGhostComm(); CHKERRQ(ierr);
  ierr = u3.begin_access(); CHKERRQ(ierr);
  ierr = v3.begin_access(); CHKERRQ(ierr);
  ierr = u3.getSurfaceValues(vusSIA, vH); CHKERRQ(ierr);
  ierr = v3.getSurfaceValues(vvsSIA, vH); CHKERRQ(ierr);
  ierr = u3.end_access(); CHKERRQ(ierr);
  ierr = v3.end_access(); CHKERRQ(ierr);

  // remove SIA computed velocity from surface velocity to get basal velocity
  PetscScalar **us, **vs, **ub, **vb, **usSIA, **vsSIA, **fofv;
  ierr =  us_in.get_array(us);    CHKERRQ(ierr);
  ierr =  vs_in.get_array(vs);    CHKERRQ(ierr);
  ierr = vusSIA.get_array(usSIA); CHKERRQ(ierr);
  ierr = vvsSIA.get_array(vsSIA); CHKERRQ(ierr);
  ierr = ub_out.get_array(ub);    CHKERRQ(ierr);
  ierr = vb_out.get_array(vb);    CHKERRQ(ierr);
  ierr =  vfofv.get_array(fofv);    CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      // (ubarSIA,vbarSIA) is on reg grid and points down surface gradient
      const PetscScalar f = fofv[i][j];
//      const PetscScalar 
//          csSIA = sqrt(PetscSqr(usSIA[i][j]) + PetscSqr(vsSIA[i][j])),
//          cs    = sqrt(PetscSqr(us[i][j]) + PetscSqr(vs[i][j]))
//      if (csSIA >= cs) {
//        // in this case all the surface speed can be accounted for by SIA
//        //   deformation in vertical planes; therefore no sliding
//        ub[i][j] = 0.0;
//        vb[i][j] = 0.0;
//      } else { // main case: remove according to fofv[i][j];

        // FIXME:  I am skeptical this will work!
        // recall eqn (21) of Bueler&Brown(2008) is
        //    U = f(|u_SSA|) u_SIA + (1-f(|u_SSA|)) u_SSA
        // where u is from SIA and v is from SSA; if f=f(|u_SSA|) is
        // assumed known then
        //    u_b = u_SSA = (Us - f u_SIA) / (1-f)
        ub[i][j] = (us[i][j] - f * usSIA[i][j]) / (1.0 - f);
        vb[i][j] = (vs[i][j] - f * vsSIA[i][j]) / (1.0 - f);

//        ub[i][j] = us[i][j] - usSIA[i][j]; // straightforward  removal might work
//        vb[i][j] = vs[i][j] - vsSIA[i][j];
//      }
    }
  }
  ierr =  us_in.end_access(); CHKERRQ(ierr);
  ierr =  vs_in.end_access(); CHKERRQ(ierr);
  ierr = ub_out.end_access(); CHKERRQ(ierr);
  ierr = vb_out.end_access(); CHKERRQ(ierr);
  ierr = vusSIA.end_access(); CHKERRQ(ierr);
  ierr = vvsSIA.end_access(); CHKERRQ(ierr);
  ierr =  vfofv.end_access(); CHKERRQ(ierr);
  
  return 0;
}


//! Compute basal shear stress from a given sliding velocity using SSA equations.
/*!
We assume the basal sliding velocity is known.  Using ice thickness, 
ice surface elevation, and an ice temperature field, we use 
SSA equations to compute the basal shear stress.

There is no smoothing in this initial implementation.

In particular:
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

This procedure calls

- assembleSSAMatrix(),
- assembleSSARhs(), and
- computeEffectiveViscosity().

It saves temporary versions of various flags modified by those routines.

The result in IceModelVec2 taubx_out, tauby_out are on the regular grid.
They must be allocated (created) before calling this routine.
 */
PetscErrorCode IceModel::computeBasalShearFromSSA(
                          IceModelVec2 ub_in, IceModelVec2 vb_in,
			  IceModelVec2 &taubx_out, IceModelVec2 &tauby_out) {
  PetscErrorCode ierr;

  // effective viscosity for ub_in, vb_in: save flags before changing them
  const PetscTruth leaveNuHAloneSSA_save = leaveNuHAloneSSA, 
                   useConstantNuHForSSA_save = useConstantNuHForSSA;
  leaveNuHAloneSSA = PETSC_FALSE;
  useConstantNuHForSSA = PETSC_FALSE;
  IceModelVec2 myvNuH[2] = {vWork2d[0], vWork2d[1]}; // already allocated space
  IceModelVec2 vubarSSAold = vWork2d[2], 
               vvbarSSAold = vWork2d[3];
  ierr = vubarSSA.copy_to(vubarSSAold); CHKERRQ(ierr);
  ierr = vvbarSSA.copy_to(vvbarSSAold); CHKERRQ(ierr);
  ierr = ub_in.copy_to(vubarSSA); CHKERRQ(ierr);
  ierr = vb_in.copy_to(vvbarSSA); CHKERRQ(ierr);
  // eps=0.0 in bdd-below regularization; Schoof type regularization does occur;
  ierr = computeEffectiveViscosity(myvNuH, 0.0); CHKERRQ(ierr);
  ierr = vubarSSAold.copy_to(vubarSSA); CHKERRQ(ierr);
  ierr = vvbarSSAold.copy_to(vvbarSSA); CHKERRQ(ierr);
  leaveNuHAloneSSA = leaveNuHAloneSSA_save;
  useConstantNuHForSSA = useConstantNuHForSSA_save;

  // allocate Mat and Vecs; compare setup in IceModel::createVecs() and
  //   IceModel::velocitySSA()
  Mat A;
  Vec x, result, rhs;
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
  ierr = assembleSSAMatrix(false, myvNuH, A); CHKERRQ(ierr);
 
  // Note rhs contains driving terms  - \rho g H \grad h
  ierr = assembleSSARhs(false, rhs); CHKERRQ(ierr);

  // Set x = [u v]'  (interleaved).
  PetscScalar **ub, **vb;
  const PetscInt  twoMy = 2 * grid.My;
  ierr = VecSet(x, 0.0); CHKERRQ(ierr);
  ierr = ub_in.get_array(ub); CHKERRQ(ierr);
  ierr = vb_in.get_array(vb); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      ierr = VecSetValue(x, i*twoMy + 2*j, ub[i][j], INSERT_VALUES); CHKERRQ(ierr);
      ierr = VecSetValue(x, i*twoMy + 2*j+1, vb[i][j], INSERT_VALUES); CHKERRQ(ierr);
    }
  }
  ierr = ub_in.end_access(); CHKERRQ(ierr);
  ierr = vb_in.end_access(); CHKERRQ(ierr);

  // assemble matrix and vec before multiplication; communicate!
  ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(x); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(x); CHKERRQ(ierr);

  // With new A and x, compute  result = Ax - rhs = L[u] - (driving)
  // HERE IS INVERSE MODEL!!
  
  ierr = MatMult(A,x,result); CHKERRQ(ierr);  
  
  ierr = VecAXPY(result,-1.0,rhs); CHKERRQ(ierr);  // note result = 0 if MASK_SHEET
  
  // transfer result to taub output vectors; compare code in
  //   IceModel::moveVelocityToDAVectors()
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

Values of both \f$\tau_c\f$ and \f$\phi\f$ are returned by this routine
but it is envisioned that the value of \f$\tau_c\f$ may evolve in a run
while \f$\phi\f$ is treated as time-independent.

Here is the context for this procedure:

This procedure should not be called in the purely plastic case, that is, 
if <tt>PlasticBasalType::pseudo_plastic</tt> is \c FALSE.

In the pseudo plastic case PlasticBasalType::drag() will compute 
the basal shear stress as
    \f[ \tau_{(b)}
           = - \frac{\tau_c}{|U|^{1-q} |U_{\mathtt{thresh}}|^q} U\f]
where \f$\tau_{(b)}=(\tau_{(b)x},\tau_{(b)y})\f$, \f$U=(u,v)\f$,
\f$q=\f$ <tt>pseudo_q</tt>, and 
\f$U_{\mathtt{thresh}}=\f$ <tt>pseudo_u_threshold</tt>.
Typical values for the constants are \f$q=0.25\f$ 
and 100 m/a for \f$U_{\mathtt{thresh}}\f$.

The linearly viscous till case \f$q=1\f$ is allowed, in which case
\f$\beta = \tau_c/|U_{\mathtt{thresh}}|\f$.

On the other hand, updateYieldStressFromHmelt() computes
the till yield stress \f$\tau_c\f$, a scalar, by 
    \f[ \tau_c = c_0 + \tan\left(\frac{\pi}{180} \phi\right)\, N \f]
where \f$\phi\f$ is the map of till friction angle, in degrees, stored in
\c vtillphi.  The till cohesion \f$c_0\f$ (= \c plastic_till_c_0) may be zero.
Here \f$N\f$ is the effective pressure on the till, namely \f$N = \rho g H - p_w\f$,
where \f$p_w\f$ is the porewater pressure.  And \f$p_w\f$ is estimated in 
terms of the stored basal water \c vHmelt; see getEffectivePressureOnTill().

In this context, the current procedure uses a sliding velocity 
and a basal shear stress (sign choice so that the stress is applied to the ice)
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

No communication occurs in this routine.

This procedure reads three fields from the current model state, namely ice 
thickness \c vH, effective thickness of basal meltwater \c vHmelt, and the 
mask \c vMask.
 */
PetscErrorCode IceModel::computeTFAFromBasalShearStressUsingPseudoPlastic(
                 IceModelVec2 ub_in, IceModelVec2 vb_in,
		 IceModelVec2 taubx_in, IceModelVec2 tauby_in, 
                 IceModelVec2 &tauc_out, IceModelVec2 &tfa_out) {
  PetscErrorCode ierr;
  
  if (doPseudoPlasticTill == PETSC_FALSE) {
    ierr = verbPrintf(1, grid.com, 
       "WARNING: computeTFAFromBasalShearStress() should only be called with q > 0.0\n"
       "  in pseudo-plastic model;  here is PlasticBasalType::printInfo() output:\n");
       CHKERRQ(ierr);
    ierr = basal->printInfo(1,grid.com); CHKERRQ(ierr);
  }
  
  const PetscScalar slowOrWrongDirFactor = 10.0,
                    sufficientSpeed = 1.0 / secpera,
                    sameDirMaxAngle = pi/4.0,
                    phiDefault = 30.0,
                    largeYieldStress = 1000.0e3;  // large; 1000 kPa = 10 bar
  PetscScalar **tauc, **phi, **ub, **vb, **taubx, **tauby, 
              **Hmelt, **H, **mask;
  ierr =    ub_in.get_array(ub);     CHKERRQ(ierr);
  ierr =    vb_in.get_array(vb);     CHKERRQ(ierr);
  ierr = taubx_in.get_array(taubx); CHKERRQ(ierr);
  ierr = tauby_in.get_array(tauby); CHKERRQ(ierr);
  ierr = tauc_out.get_array(tauc);  CHKERRQ(ierr);
  ierr =  tfa_out.get_array(phi);   CHKERRQ(ierr);
  ierr =   vHmelt.get_array(Hmelt); CHKERRQ(ierr);
  ierr =       vH.get_array(H);     CHKERRQ(ierr);
  ierr =    vMask.get_array(mask);  CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (modMask(mask[i][j]) == MASK_FLOATING) {
        tauc[i][j] = 0.0;
        phi[i][j] = phiDefault;
      } else if (H[i][j] <= 0.0) { // if no ice use large resistance
        tauc[i][j] = largeYieldStress;
        phi[i][j] = phiDefault;
      } else { // case: grounded and ice is present
        // get tauc = yield stress
        const PetscScalar
            speed = sqrt(PetscSqr(ub[i][j]) + PetscSqr(vb[i][j])),
            taubmag = sqrt(PetscSqr(taubx[i][j]) + PetscSqr(tauby[i][j]));
        if (speed < sufficientSpeed) {
          tauc[i][j] = slowOrWrongDirFactor * taubmag;
        } else {
          const PetscScalar 
              ubDOTtaub = ub[i][j] * taubx[i][j] + tauby[i][j] * vb[i][j],
              angle = acos(ubDOTtaub / (speed * taubmag));
          if (PetscAbs(angle) > sameDirMaxAngle) {
            tauc[i][j] = slowOrWrongDirFactor * taubmag;
          } else {
            // use the formula which inverts PlasticBasalType::drag()
            tauc[i][j] = basal->taucFromMagnitudes(taubmag, speed);
          }
        }
        // get phi = till friction angle
        if (tauc[i][j] > plastic_till_c_0) {
          const PetscScalar N = getEffectivePressureOnTill(H[i][j], Hmelt[i][j]);
          phi[i][j] = (180.0 / pi) * atan( (tauc[i][j] - plastic_till_c_0) / N );
        } else {
          phi[i][j] = 0.0;
        }
      }
    }
  }
  ierr =    ub_in.end_access(); CHKERRQ(ierr);
  ierr =    vb_in.end_access(); CHKERRQ(ierr);
  ierr = taubx_in.end_access(); CHKERRQ(ierr);
  ierr = tauby_in.end_access(); CHKERRQ(ierr);
  ierr = tauc_out.end_access(); CHKERRQ(ierr);
  ierr =  tfa_out.end_access(); CHKERRQ(ierr);
  ierr =   vHmelt.end_access(); CHKERRQ(ierr);
  ierr =       vH.end_access(); CHKERRQ(ierr);
  ierr =    vMask.end_access(); CHKERRQ(ierr);

  return 0;
}
