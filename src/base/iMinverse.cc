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


//! Invert given ice surface velocities to find till friction angle using a pseudo-plastic model.
/*!
Reads user option <tt>-surf_vel_to_tfa foo.nc</tt>.  If this option is \e not found,
this method returns without doing anything.

If the option is found, then expects to find variables
\c us and \c vs, namely x and y components of observed surface velocity,
in NetCDF file \c foo.nc.  If found, reads them.  (Note these
"observed" values must have full area coverage with no missing values, so
in practice these "observed" values will have a modeling or kriging procedure
which fills holes in coverage.)

First calls computeSIASurfaceVelocity() to get the value of the surface
velocity which would be computed from SIA flow model for the current
geometry and temperature.

If <tt>-super</tt> is also set, we call the method computeFofVforInverse()
which computes the factor \f$f(|\mathbf{v}|)\f$ used in combining the SIA 
and SSA velocities \lo\cite[equation (21)]{BBssasliding}\elo.  If
<tt>-super</tt> is not set, we do not call computeFofVforInverse(), and
instead just set \f$f(|\mathbf{v}|)=0\f$.  (Thus the observed surface velocity
is treated as all described by the SSA model.)

This method then calls the three steps in turn:

- removeSIApart() to compute the basal sliding velocity,

- computeBasalShearFromSSA() to compute the basal shear stress according to 
the SSA stress balance, and

- computeTFAFromBasalShearStressUsingPseudoPlastic().

There is no smoothing in this initial implementation.

The current method must be called after initialization is essentially completed.
It uses geometry, temperature field, and map of effective thickness of basal till
water in the inversion.

The result is a field of values for till yield stress and the till friction angle in a 
pseudo-plastic till model.  This model includes purely-plastic and linearly-viscous cases.

PISM can then do the forward time-stepping model from the resulting till friction angle
field.
 */
PetscErrorCode IceModel::invertSurfaceVelocities(const PetscTruth writeWhenDone) {
  PetscErrorCode ierr;
  
  PetscTruth svTOtfaSet, superSet;
  char filename[PETSC_MAX_PATH_LEN], fofvname[PETSC_MAX_PATH_LEN];
  ierr = PetscOptionsGetString(PETSC_NULL, "-surf_vel_to_tfa", filename, 
                               PETSC_MAX_PATH_LEN, &svTOtfaSet); CHKERRQ(ierr);
  if (svTOtfaSet == PETSC_FALSE) {
    return 0;  // leave if you are not wanted ...
  }
  ierr = PetscOptionsHasName(PETSC_NULL, "-super", &superSet); CHKERRQ(ierr);

  IceModelVec2 usIn, vsIn, usSIA, vsSIA, taubxComputed, taubyComputed, fofv;

  // create inverse model state variables with metadata
  ierr = usIn.create(grid, "us", true); CHKERRQ(ierr);
  ierr = usIn.set_attrs(
     NULL, "x component of ice surface velocity","m s-1", NULL); CHKERRQ(ierr);
  ierr = vsIn.create(grid, "vs", true); CHKERRQ(ierr);
  ierr = vsIn.set_attrs(
     NULL, "y component of ice surface velocity","m s-1", NULL); CHKERRQ(ierr);
  ierr = usSIA.create(grid, "usSIA", false); CHKERRQ(ierr);  // global
  ierr = usSIA.set_attrs(
     NULL, "x component of ice surface velocity predicted by non-sliding SIA",
     "m s-1", NULL); CHKERRQ(ierr);
  ierr = vsSIA.create(grid, "vsSIA", false); CHKERRQ(ierr);  // global
  ierr = vsSIA.set_attrs(
     NULL, "y component of ice surface velocity predicted by non-sliding SIA",
     "m s-1", NULL); CHKERRQ(ierr);
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

  // compute f(|v|) factor, or set to constant if -super not used
  if (superSet == PETSC_TRUE) {
    ierr = verbPrintf(2, grid.com, 
       "option -super seen;  computing f(|v|) from observed surface velocities\n"
       "  by solving transcendental equations ...\n",
       fofvname); CHKERRQ(ierr);
    ierr = computeFofVforInverse(usIn, vsIn, usSIA, vsSIA, fofv); CHKERRQ(ierr);
  } else {
    // compute an f(|v|) field from current ubarSSA, vbarSSA
    ierr = verbPrintf(2, grid.com, 
       "option -super NOT seen;  setting f(|v|) to 0.0, so NONE of SIA velocity\n"
       "  is removed from observed surface velocity ...\n");
       CHKERRQ(ierr);
    ierr = fofv.set(0.0); CHKERRQ(ierr);
  }

  // do three steps of inverse model; result is tauc and tillphi fields
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
           "  computing till friction angle using (pseudo-)plastic model ...\n"); 
           CHKERRQ(ierr);
  ierr = computeTFAFromBasalShearStressUsingPseudoPlastic(
           vubarSSA, vvbarSSA, taubxComputed, taubyComputed, vtauc, vtillphi);
           CHKERRQ(ierr);

  // write out stored inverse info for user's edification; mostly for debug
  if (writeWhenDone == PETSC_TRUE) {
    char filename[PETSC_MAX_PATH_LEN];
    PetscTruth o_set;
    ierr = PetscOptionsGetString(PETSC_NULL, "-o", filename, 
                                 PETSC_MAX_PATH_LEN, &o_set); CHKERRQ(ierr);
    if (!o_set) {
      ierr = PetscPrintf(grid.com,
		"WARNING: inverse model can only write to file if -o option given,\n"
		"  so no writing of inverse model fields; sorry ...\n");
		CHKERRQ(ierr);
    } else {
      ierr =          usIn.write(filename, NC_DOUBLE); CHKERRQ(ierr);
      ierr =          vsIn.write(filename, NC_DOUBLE); CHKERRQ(ierr);
      ierr = taubxComputed.write(filename, NC_DOUBLE); CHKERRQ(ierr);
      ierr = taubyComputed.write(filename, NC_DOUBLE); CHKERRQ(ierr);
      ierr =          fofv.write(filename, NC_DOUBLE); CHKERRQ(ierr);
    }
  }
  
  // clean up
  ierr =          usIn.destroy(); CHKERRQ(ierr);
  ierr =          vsIn.destroy(); CHKERRQ(ierr);
  ierr = taubxComputed.destroy(); CHKERRQ(ierr);
  ierr = taubyComputed.destroy(); CHKERRQ(ierr);
  ierr =          fofv.destroy(); CHKERRQ(ierr);
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
equation by numerical root-finding.  Define
  \f[ F(x) = \frac{2}{\pi} \arctan\left(\frac{x}{v_0^2}\right) \f]
so \f$f(|\mathbf{v}|) = 1 - F(x)\f$ and \f$x\f$ solves
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
\f$G'(x)\f$.  The initial guess
  \f[x_0 = |\mathbf{u}_s|^2 + 4 (\mathbf{U}_s - \mathbf{u}_s) \cdot \mathbf{U}_s \f]
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


//! Compute basal shear stress from a given sliding velocity using SSA equations.
/*!
This is one of several routines called by invertSurfaceVelocities() for
inverse model-based initialization.  See comments for that method.

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

- assembleSSARhs(),

- computeEffectiveViscosity().

It saves temporary versions of various flags modified by those routines.

The result in IceModelVec2 taubx_out, tauby_out must be allocated (created) 
before calling this routine.
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
This is one of several routines called by invertSurfaceVelocities() for
inverse model-based initialization.  See comments for that method.

The more complete brief description might be: "Compute till friction angle in 
degrees from saved basal shear stress map and saved basal velocity 
(both vector-valued) using pseudo-plastic till constitutive relation 
and a pore water pressure model."

FIXME:  This procedure should be minimizing a regularized L^2 norm
(i.e. doing least squares).

Values of both \f$\tau_c\f$ and \f$\phi\f$ are returned by this routine
but it is envisioned that the value of \f$\tau_c\f$ may evolve in a run
while \f$\phi\f$ is treated as time-independent.

Here is the context for this procedure:

This procedure should not be called in the purely plastic case, that is, 
if <tt>PlasticBasalType::pseudo_plastic</tt> is \c FALSE.

In the pseudo plastic case PlasticBasalType::drag() will compute 
the basal shear stress as
    \f[ \tau_{(b)}
           = - \frac{\tau_c}{|U|^{1-q} |U_{\mathtt{thresh}}|^q} U    \f]
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
