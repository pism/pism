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

#include <petscda.h>
#include "nc_util.hh"
#include "LocalInterpCtx.hh"
#include "iceModelVec.hh"
#include "iceModel.hh"

/*
example usage:

[START from P1.nc in examples/pst/pst.sh]

pisms -pst -P1 -if P1.nc -y 10 -f3d -o P1plus10.nc  # writes uvelsurf,vvelsurf in addition to 3D

[NOW USE NCO to convert units on uvelsurf,vvelsurf; saved in m a-1, needed in m s-1]

pisms -pst -P1 -if P1plus10.nc -y 1 -surf_vel_to_tfa P1plus10.nc -write_inverse_fields foo.nc

*/

//! Invert given ice surface velocities to find till friction angle using a pseudo-plastic model.
/*!
Reads user option <tt>-surf_vel_to_tfa foo.nc</tt>.  If this option is \e not found,
this method returns without doing anything.

If the option is found, then expects to find variables
\c uvelsurf and \c vvelsurf, namely x and y components of observed surface velocity,
in NetCDF file \c foo.nc.  If found, reads them.

Note that these "observed" values must have full area coverage with no missing 
values, so in practice they will come from a modeling or kriging procedure 
which fills holes in coverage.

Then calls computeSIASurfaceVelocity() to get the value of the surface
velocity which would be computed from SIA flow model for the current
geometry and temperature.

If <tt>-super</tt> is also set, we call the method computeFofVforInverse()
which computes the factor \f$f(|\mathbf{v}|)\f$ used in combining the SIA 
and SSA velocities \lo\cite[equation (21)]{BBssasliding}\elo.  If
<tt>-super</tt> is not set, we do not call computeFofVforInverse(), and
instead just set \f$f(|\mathbf{v}|)=0\f$, so the observed surface velocity
is treated as all described by the SSA model.

This method then calls four steps in turn:

- removeSIApart() to compute the basal sliding velocity \f$\mathbf{u}_b\f$ by removing the SIA
modeled shear from the surface velocity,

- computeBasalShearFromSSA() to compute the basal shear stress \f$\tau_b\f$ which solves 
the SSA stress balance in which the basal sliding velocity \f$\mathbf{u}_b\f$ is used as the
depth independent velocity,

- computeYieldStressFromBasalShearUsingPseudoPlastic() to compute the yield stress \f$\tau_c\f$
for the given basal shear stress \f$\tau_b\f$ and basal sliding velocity \f$\mathbf{u}_b\f$,

- computeTFAFromYieldStress() to compute the till friction angle \f$\phi\f$ from the 
yield stress \f$\tau_c\f$.

The result is a field of values for till friction angle in a pseudo-plastic till model.  
PISM can do the forward time-stepping model using the resulting till friction angle
field.

This model includes purely-plastic and linearly-viscous cases.  There is no 
smoothing in this initial implementation.  This method must be called after 
initialization is essentially completed.  It uses geometry, temperature field, 
and map of effective thickness of basal till water in the inversion.
 */
PetscErrorCode IceModel::invertSurfaceVelocities() {
  PetscErrorCode ierr;
  
  PetscTruth svTOtfaSet,invfieldsSet;
  char filename[PETSC_MAX_PATH_LEN], invfieldsname[PETSC_MAX_PATH_LEN];
  ierr = PetscOptionsGetString(PETSC_NULL, "-surf_vel_to_tfa", filename, 
                               PETSC_MAX_PATH_LEN, &svTOtfaSet); CHKERRQ(ierr);
  if (svTOtfaSet == PETSC_FALSE) {
    return 0;  // leave if you are not wanted ...
  }
  ierr = PetscOptionsGetString(PETSC_NULL, "-write_inverse_fields", invfieldsname, 
                               PETSC_MAX_PATH_LEN, &invfieldsSet); CHKERRQ(ierr);
  // note computeYieldStressFromBasalShearUsingPseudoPlastic() reads additional
  //   option "-inverse_reg_eps"

  IceModelVec2 usIn, vsIn, usSIA, vsSIA, taubxComputed, taubyComputed, fofv, taucComputed;

  // create inverse model state variables with metadata
  ierr = usIn.create(grid, "uvelsurf", true); CHKERRQ(ierr);
  ierr = usIn.set_attrs(
     NULL, "x component of velocity of ice at ice surface","m s-1", NULL); CHKERRQ(ierr);
  ierr = vsIn.create(grid, "vvelsurf", true); CHKERRQ(ierr);
  ierr = vsIn.set_attrs(
     NULL, "y component of velocity of ice at ice surface","m s-1", NULL); CHKERRQ(ierr);
  ierr = usSIA.create(grid, "uvelsurfSIA", false); CHKERRQ(ierr);  // global
  ierr = usSIA.set_attrs(
     NULL, "x component of ice surface velocity predicted by non-sliding SIA",
     "m s-1", NULL); CHKERRQ(ierr);
  ierr = vsSIA.create(grid, "vvelsurfSIA", false); CHKERRQ(ierr);  // global
  ierr = vsSIA.set_attrs(
     NULL, "y component of ice surface velocity predicted by non-sliding SIA",
     "m s-1", NULL); CHKERRQ(ierr);
  ierr = taubxComputed.create(grid, "taubxOUT", false);  // global
  ierr = taubxComputed.set_attrs(
     NULL, "inverse-model-computed x component of basal shear stress", 
     "Pa", NULL); CHKERRQ(ierr);
  ierr = taubyComputed.create(grid, "taubyOUT", false);  // global
  ierr = taubyComputed.set_attrs(
     NULL, "inverse-model-computed y component of basal shear stress", 
     "Pa", NULL); CHKERRQ(ierr);
  ierr = fofv.create(grid, "fofv", false);  // global
  ierr = fofv.set_attrs(
     NULL, "inverse-model-computed fraction of final velocity which comes from SIA in hybrid model", 
     "1", NULL); CHKERRQ(ierr);
  ierr = taucComputed.create(grid, "taucOUT", false);  // global
  ierr = taucComputed.set_attrs(
     NULL, "inverse-model-computed yield stress for basal till (plastic or pseudo-plastic model)", 
     "Pa", NULL); CHKERRQ(ierr);

  // read in surface velocity
  ierr = verbPrintf(2, grid.com, 
     "option -surf_vel_to_tfa seen;  doing ad hoc inverse model;\n"
     "  reading observed surface velocities us, vs from file %s ...\n",
     filename); CHKERRQ(ierr);
  bool file_exists = false; // check whether file exists
  NCTool nc(&grid);
  ierr = nc.open_for_reading(filename, file_exists); CHKERRQ(ierr);
  if (!file_exists) {
    ierr = PetscPrintf(grid.com,"PISM ERROR: Can't open file '%s'.\n",
    filename); CHKERRQ(ierr);    PetscEnd();
  }
  grid_info g;
  ierr = nc.get_grid_info_2d(g); CHKERRQ(ierr);
  ierr = nc.close(); CHKERRQ(ierr);
  LocalInterpCtx lic(g, NULL, NULL, grid); // 2D only
  // read in by regridding; checks for existence of file and variable in file
  ierr = usIn.regrid(filename, lic, true); CHKERRQ(ierr);// it *is* critical
  ierr = vsIn.regrid(filename, lic, true); CHKERRQ(ierr);

  // compute the surface velocity which the SIA predicts
  ierr = verbPrintf(2, grid.com, 
    "computing SIA surface velocity (for removal from observed before inverting w SSA ...\n");
    CHKERRQ(ierr);
  ierr = computeSIASurfaceVelocity(usSIA, vsSIA); CHKERRQ(ierr);

#if 1
  // compute f(|v|) factor, or set to constant if -super not used
  if (doSuperpose == PETSC_TRUE) {
    ierr = verbPrintf(2, grid.com, 
       "  flag doSuperpose (option -super) seen;  computing f(|v|) from observed\n"
       "    surface velocities by solving transcendental equations ...\n",
       fofvname); CHKERRQ(ierr);
    ierr = computeFofVforInverse(usIn, vsIn, usSIA, vsSIA, fofv); CHKERRQ(ierr);
  } else {
    // compute an f(|v|) field from current ubarSSA, vbarSSA
    ierr = verbPrintf(2, grid.com, 
       "  option -super NOT seen;  setting f(|v|) to 0.0, so NONE of SIA velocity\n"
       "    is removed from observed surface velocity ...\n");
       CHKERRQ(ierr);
    ierr = fofv.set(0.0); CHKERRQ(ierr);
  }

  // do four steps of inverse model; result is tillphi field
  ierr = verbPrintf(2, grid.com, 
           "  removing SIA-computed (vertical shear) part of surface velocity ...\n");
           CHKERRQ(ierr);
  ierr = removeSIApart(
           usIn, vsIn, usSIA, vsSIA, fofv, vubarSSA, vvbarSSA);    CHKERRQ(ierr);

  ierr = verbPrintf(2, grid.com, 
           "  computing membrane stresses and basal shear stress from SSA;\n"
           "    (applying SSA differential operator ...)\n");
           CHKERRQ(ierr);
  ierr = computeBasalShearFromSSA(
           vubarSSA, vvbarSSA, taubxComputed, taubyComputed);   CHKERRQ(ierr);

  ierr = verbPrintf(2, grid.com, 
           "  computing till yield stress tau_c using (pseudo-)plastic model ...\n"); 
           CHKERRQ(ierr);
  ierr = computeYieldStressFromBasalShearUsingPseudoPlastic(
           vubarSSA, vvbarSSA, taubxComputed, taubyComputed, taucComputed); CHKERRQ(ierr);

  ierr = verbPrintf(2, grid.com, 
           "  computing till friction angle phi from tau_c by Mohr-Coulomb criterion ...\n"); 
           CHKERRQ(ierr);
  ierr = computeTFAFromYieldStress(taucComputed,vtillphi); CHKERRQ(ierr);
#endif

  // write out stored inverse info for user's edification; mostly for debug
  if (invfieldsSet == PETSC_TRUE) {
    // FIXME: message that we are writing to file here
    // prepare the file
    ierr = nc.open_for_writing(invfieldsname,PETSC_TRUE); CHKERRQ(ierr);
    ierr = nc.append_time(grid.year * secpera); CHKERRQ(ierr);
    ierr = nc.write_history("-write_inverse_fields option read"); CHKERRQ(ierr);
    ierr = nc.write_polar_stereographic(psParams.svlfp, psParams.lopo, psParams.sp); CHKERRQ(ierr);
    ierr = nc.write_global_attrs(PETSC_FALSE, "CF-1.0"); CHKERRQ(ierr);
    ierr = nc.close(); CHKERRQ(ierr);
    // write the fields
    ierr =          usIn.write(invfieldsname, NC_DOUBLE); CHKERRQ(ierr);
    ierr =          vsIn.write(invfieldsname, NC_DOUBLE); CHKERRQ(ierr);
    ierr =         usSIA.write(invfieldsname, NC_DOUBLE); CHKERRQ(ierr);
    ierr =         vsSIA.write(invfieldsname, NC_DOUBLE); CHKERRQ(ierr);
    ierr = taubxComputed.write(invfieldsname, NC_DOUBLE); CHKERRQ(ierr);
    ierr = taubyComputed.write(invfieldsname, NC_DOUBLE); CHKERRQ(ierr);
    ierr =          fofv.write(invfieldsname, NC_DOUBLE); CHKERRQ(ierr);
    ierr =  taucComputed.write(invfieldsname, NC_DOUBLE); CHKERRQ(ierr);
  }
  
  // clean up
  ierr =          usIn.destroy(); CHKERRQ(ierr);
  ierr =          vsIn.destroy(); CHKERRQ(ierr);
  ierr = taubxComputed.destroy(); CHKERRQ(ierr);
  ierr = taubyComputed.destroy(); CHKERRQ(ierr);
  ierr =          fofv.destroy(); CHKERRQ(ierr);
  ierr =  taucComputed.destroy(); CHKERRQ(ierr);
  return 0;
}


//! Compute the surface velocity from only the vertical plane shear predicted by the nonsliding SIA.
/*!
This is one of several routines called by invertSurfaceVelocities() for
inverse model-based initialization.  See comments for that method.

The nonsliding thermomechanically coupled SIA equations are used.

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
 */
PetscErrorCode IceModel::computeSIASurfaceVelocity(
		 IceModelVec2 &usSIA_out, IceModelVec2 &vsSIA_out) {
  PetscErrorCode ierr;

  ierr = surfaceGradientSIA(); CHKERRQ(ierr); // comm may happen here ...
  // surface gradient temporarily stored in vWork2d[0..3] ...

  // compute the vertically integrated velocity from the nonsliding SIA
  //   (namely vuvbar[2], which is (ubar,vbar) on staggered grid)
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
  ierr = vub.set(0.0); CHKERRQ(ierr);  // no sliding in SIA
  ierr = vvb.set(0.0); CHKERRQ(ierr);
  ierr = horizontalVelocitySIARegular(); CHKERRQ(ierr);

  // now grab surface value of 3D velocity fields
  ierr = u3.beginGhostComm(); CHKERRQ(ierr);
  ierr = v3.beginGhostComm(); CHKERRQ(ierr);
  ierr = u3.endGhostComm(); CHKERRQ(ierr);
  ierr = v3.endGhostComm(); CHKERRQ(ierr);
  ierr = u3.begin_access(); CHKERRQ(ierr);
  ierr = v3.begin_access(); CHKERRQ(ierr);
  ierr = u3.getSurfaceValues(usSIA_out, vH); CHKERRQ(ierr);
  ierr = v3.getSurfaceValues(vsSIA_out, vH); CHKERRQ(ierr);
  ierr = u3.end_access(); CHKERRQ(ierr);
  ierr = v3.end_access(); CHKERRQ(ierr);
  
  return 0;
}


//! Computes G(x) and its derivative needed by computeFofVforInverse().
/*!
See comments for computeFofVforInverse() for formula for \f$G(x)\f$.
 */
PetscErrorCode IceModel::getGforInverse(const PetscScalar x, 
                   const PetscScalar UsuSIAdiffsqr, const PetscScalar UsuSIAdiffdotuSIA,
                   const PetscScalar uSIAsqr, PetscScalar &G, PetscScalar &Gprime) {

  const PetscScalar  v0sqr  = PetscSqr(100.0 / secpera),
                     outC   = 2.0 / pi,
                     F      = outC * atan(x / v0sqr),
                     Fprime = outC * (1.0 / (1.0 + PetscSqr(x / v0sqr))) / v0sqr;
  G      = (x - uSIAsqr) * F * F - 2.0 * F * UsuSIAdiffdotuSIA - UsuSIAdiffsqr;
  Gprime = F * F + 2.0 * Fprime * ( (x - uSIAsqr) * F - UsuSIAdiffdotuSIA );
  return 0;
}


//! Given observed surface velocity and SIA computed surface velocity, find f(|v|) factor.
/*!
This is one of several routines called by invertSurfaceVelocities() for
inverse model-based initialization.  See comments for that method.

The scalar field which gets computed is called \f$f(|\mathbf{v}|)\f$ in 
\lo\cite[equations (21) and (22)]{BBssasliding}\elo.  In fact, if 
\f$\mathbf{U}_s\f$ is the observed surface velocity then
at each point in the map plane we seek \f$f(|\mathbf{v}|)\f$ in the equation
  \f[ \mathbf{U}_s = f(|\mathbf{v}|) \mathbf{u}_s
                   + \left(1 - f(|\mathbf{v}|)\right) \mathbf{v} \f]
where \f$\mathbf{u}_s\f$ is the surface value of the nonsliding SIA.  Thus 
\f$\mathbf{u}_s\f$ is directly computable from the temperature and geometry,
and \f$\mathbf{v}\f$ is the depth-independent velocity which goes in the 
SSA stress balance.  Recall
  \f[ f(|\mathbf{v}|) 
         = 1 - \frac{2}{\pi} \arctan\left(\frac{|\mathbf{v}|^2}{v_0^2}\right) \f]
where \f$v_0=100\f$ m/a.

Our approach here is to find \f$x = |\mathbf{v}|^2\f$ by solving a transcendental
equation by numerical root-finding.  Define \f$F(x) = (2/\pi) \arctan\left(x\,v_0^{-2}\right)\f$
so \f$f(|\mathbf{v}|) = 1 - F(x)\f$ and the scalar \f$x\f$ solves the scalar equation
  \f[ |\mathbf{U}_s - (1-F(x))\mathbf{u}_s|^2 = F(x)^2 x. \f]
We treat this as a root-finding problem \f$G(x)=0\f$ for \f$x\f$, where, with 
a little more rewriting,
  \f[ G(x) = (x - |\mathbf{u}_s|^2) F(x)^2 
             - 2 F(x) (\mathbf{U}_s - \mathbf{u}_s) \cdot \mathbf{u}_s
             - |\mathbf{U}_s - \mathbf{u}_s|^2. \f]
Note \f$F(0)=0\f$ and \f$F(x)\to 1\f$ as \f$x\to\infty\f$.  Thus \f$G(0)\le 0\f$ 
and \f$G(x) \sim x\f$ as \f$x\to\infty\f$.  Because \f$G(x)\f$ is continuous, there is
a nonnegative root.

We use the derivative
  \f[ G'(x) = F(x)^2 
              + 2 F'(x) \left[(x - |\mathbf{u}_s|^2) F(x) 
                              - (\mathbf{U}_s - \mathbf{u}_s) \cdot \mathbf{u}_s\right]\f]
in Newton's method
  \f[x_{n+1} = x_n - G(x_n)/G'(x_n).\f]

The method getGforInverse() actually computes \f$G(x)\f$ and 
\f$G'(x)\f$.  The initial guess is \f$x_0 = |\mathbf{u}_s|^2 
+ 4 (\mathbf{U}_s - \mathbf{u}_s) \cdot \mathbf{U}_s \f$ because that value
solves the equation that results from replacing the \f$F(x)\f$ term in the equation
\f$G(x)=0\f$ above by 1/2.  That is, \f$x_0\f$ would be the root of \f$G(x)\f$ if
\f$F(x)\f$ were equal to 1/2.

We iterate Newton's method to \f$10^{-12}\f$ relative tolerance and warn on nonconvergence.
 */
PetscErrorCode IceModel::computeFofVforInverse(
                            IceModelVec2 us_in, IceModelVec2 vs_in, 
                            IceModelVec2 usSIA_in, IceModelVec2 vsSIA_in, 
                            IceModelVec2 &fofv_out) {
  PetscErrorCode ierr;
  PetscScalar **us, **vs, **usSIA, **vsSIA, **f;

  ierr =    us_in.get_array(us);    CHKERRQ(ierr);
  ierr =    vs_in.get_array(vs);    CHKERRQ(ierr);
  ierr = usSIA_in.get_array(usSIA); CHKERRQ(ierr);
  ierr = vsSIA_in.get_array(vsSIA); CHKERRQ(ierr);
  ierr = fofv_out.get_array(f);     CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      // Newton's method:
      const PetscScalar 
         Usmus_x = us[i][j] - usSIA[i][j],
         Usmus_y = vs[i][j] - vsSIA[i][j];
      const PetscScalar
         us2        = PetscSqr(usSIA[i][j]) + PetscSqr(vsSIA[i][j]),
         Usmus2     = PetscSqr(Usmus_x)     + PetscSqr(Usmus_y),
         Usmusdotus = Usmus_x * usSIA[i][j] + Usmus_y * vsSIA[i][j],
         UsmusdotUs = Usmus_x * us[i][j]    + Usmus_y * vs[i][j];
      PetscScalar xold = us2 + 4.0 * UsmusdotUs;  // initial guess solves nearby eqn
      PetscScalar xnew, G, Gprime;
      PetscInt k;
      for (k=0; k<20; ++k) {
         ierr = getGforInverse(xold, Usmus2, Usmusdotus, us2, G, Gprime); CHKERRQ(ierr);
         xnew = xold - G / Gprime;  // Newton; practicing unprotected steps here ...
         const PetscScalar rel = PetscAbs(xnew - xold) / PetscAbs(xold);
         if (rel < 1.0e-12)  break;
         xold = xnew;
      }
      if (k >= 10) {
        ierr = verbPrintf(1, grid.com, 
           "WARNING: failure of Newton to converge in IceModel::computeFofVforInverse();\n"
           "  info: i=%d, j=%d, k=%d, Usmus2=%f, Usmusdotus=%f, us2=%f, xold=%f, xnew=%f\n",
           i,j,k,Usmus2,Usmusdotus,us2,xold,xnew); CHKERRQ(ierr);
      }
      // at this point xnew contains  x=|v|^2
      const PetscScalar  v0sqr  = PetscSqr(100.0 / secpera),
                         outC   = 2.0 / pi,
                         F      = outC * atan(xnew / v0sqr);
      if ( (F < 0.0) || (F > 1.0) ) {
        ierr = verbPrintf(1, grid.com, 
           "ERROR: failure to compute reasonable f(|v|) in [0,1] in\n"
           "  IceModel::computeFofVforInverse(); info: i=%d, j=%d, f[i][j]=%f\n",
           i,j,1.0 - F); CHKERRQ(ierr);
        PetscEnd();
      }
      f[i][j] = 1.0 - F;
    }
  }
  ierr = fofv_out.end_access(); CHKERRQ(ierr);
  ierr =    us_in.end_access(); CHKERRQ(ierr);
  ierr =    vs_in.end_access(); CHKERRQ(ierr);
  ierr = usSIA_in.end_access(); CHKERRQ(ierr);
  ierr = vsSIA_in.end_access(); CHKERRQ(ierr);

  return 0;
}


//! Remove the SIA surface velocity from ice surface velocity to give depth-independent, membrane stress-balanced (SSA) velocity.
/*!
This is one of several routines called by invertSurfaceVelocities() for
inverse model-based initialization.  See comments for that method.

This procedure computes a z-independent horizontal ice velocity, conceptually
a basal sliding velocity, from a specified surface horizontal ice velocity
(map).  The portion already found from SIA equations is removed.

The results in IceModelVec2 s ub_out and vb_out
must be allocated (created) before calling this procedure.
 */
PetscErrorCode IceModel::removeSIApart(
                 IceModelVec2 us_in, IceModelVec2 vs_in, 
                 IceModelVec2 usSIA_in, IceModelVec2 vsSIA_in, 
                 IceModelVec2 fofv_in,
		 IceModelVec2 &ub_out, IceModelVec2 &vb_out) {
  PetscErrorCode ierr;

  // max amount before div by (1-f) will be bypassed:
  const PetscScalar maxfofv = 0.98;

  PetscScalar **us, **vs, **ub, **vb, **usSIA, **vsSIA, **fofv;
  ierr =    us_in.get_array(us);    CHKERRQ(ierr);
  ierr =    vs_in.get_array(vs);    CHKERRQ(ierr);
  ierr = usSIA_in.get_array(usSIA); CHKERRQ(ierr);
  ierr = vsSIA_in.get_array(vsSIA); CHKERRQ(ierr);
  ierr =  fofv_in.get_array(fofv);  CHKERRQ(ierr);
  ierr =   ub_out.get_array(ub);    CHKERRQ(ierr);
  ierr =   vb_out.get_array(vb);    CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      const PetscScalar f = fofv[i][j];
      if (f < maxfofv) {
        ub[i][j] = (us[i][j] - f * usSIA[i][j]) / (1.0 - f);
        vb[i][j] = (vs[i][j] - f * vsSIA[i][j]) / (1.0 - f);
      } else {
        // if f is approximately 1, just treat all observed surface
        //   velocity as SIA, so no sliding
        ub[i][j] = 0.0;
        vb[i][j] = 0.0;
      }
    }
  }
  ierr =    us_in.end_access(); CHKERRQ(ierr);
  ierr =    vs_in.end_access(); CHKERRQ(ierr);
  ierr = usSIA_in.end_access(); CHKERRQ(ierr);
  ierr = vsSIA_in.end_access(); CHKERRQ(ierr);
  ierr =  fofv_in.end_access(); CHKERRQ(ierr);
  ierr =   ub_out.end_access(); CHKERRQ(ierr);
  ierr =   vb_out.end_access(); CHKERRQ(ierr);
  
  return 0;
}


//! Compute till friction angle in degrees given the (scalar) till yield stress.
/*!
The more complete brief description might be: "Compute till friction angle 
\f$\phi\f$ in degrees from pre-computed till yield stress \f$\tau_c\f$ 
and modeled pore water pressure."

This is one of several routines called by invertSurfaceVelocities() for
inverse model-based initialization.  See comments for that method.

It is envisioned that in forward modeling the value of \f$\tau_c\f$ evolves
while \f$\phi\f$ is treated as time-independent.

The scalar till yield stress \f$\tau_c\f$ is related to the till friction
angle by the Mohr-Coulomb criterion \lo\cite{Paterson}\elo,
    \f[ \tau_c = c_0 + \tan\left(\frac{\pi}{180} \phi\right)\,\left(\rho g H - p_w\right) \f]
Here \f$\phi\f$ is the map of till friction angle, in degrees, stored in
\c vtillphi.  The till cohesion \f$c_0\f$ (= \c plastic_till_c_0) may be zero.

The effective pressure on the till, namely \f$N = \rho g H - p_w\f$,
where \f$p_w\f$ is the porewater pressure, is computed by 
getEffectivePressureOnTill().  In summary, \f$p_w\f$ is estimated in 
terms of the stored basal water \c vHmelt which is the time-integrated basal
melt rate.

This procedure computes
    \f[ \phi = \frac{180}{\pi} \arctan\left(\frac{\tau_c - c_0}{N}\right) \f]
If \f$\tau_c - c_0 \le 0\f$ then the till friction angle is set to 
0 degrees.  If the effective pressure \f$N\f$ is smaller than 1/10 of the 
overburden then the till friction angle is set (arbitrarily) to twenty degrees.

No communication occurs in this routine.  This procedure reads the ice thickness
\c vH and the effective thickness of basal meltwater \c vHmelt from the 
current model state.

Compare updateYieldStressFromHmelt() which uses the Mohr-Coulomb criterion
to compute yield stress from till friction angle in a forward model.
 */
PetscErrorCode IceModel::computeTFAFromYieldStress(
                 IceModelVec2 tauc_in, IceModelVec2 &tfa_out) {
  PetscErrorCode ierr;
  
  const PetscScalar phiDefault = 20.0; // degrees
  PetscScalar **tauc, **phi, **Hmelt, **H;
  ierr = tauc_in.get_array(tauc);  CHKERRQ(ierr);
  ierr = tfa_out.get_array(phi);   CHKERRQ(ierr);
  ierr =  vHmelt.get_array(Hmelt); CHKERRQ(ierr);
  ierr =      vH.get_array(H);     CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (tauc[i][j] > plastic_till_c_0) {
        const PetscScalar N = getEffectivePressureOnTill(H[i][j], Hmelt[i][j]);
        if (N > 0.1 * ice->rho * grav * H[i][j]) {
          phi[i][j] = (180.0 / pi) * atan( (tauc[i][j] - plastic_till_c_0) / N );
        } else {
          phi[i][j] = phiDefault;
        }
      } else {
        phi[i][j] = 0.0;
      }
    }
  }
  ierr = tauc_in.end_access(); CHKERRQ(ierr);
  ierr = tfa_out.end_access(); CHKERRQ(ierr);
  ierr =  vHmelt.end_access(); CHKERRQ(ierr);
  ierr =      vH.end_access(); CHKERRQ(ierr);

  return 0;
}

