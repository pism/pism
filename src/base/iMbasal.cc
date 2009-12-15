// Copyright (C) 2004--2009 Jed Brown, Ed Bueler and Constantine Khroulev
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

#include <cmath>
#include <petscda.h>
#include "iceModel.hh"

/*
This file collects procedures related to SSA-as-sliding law in grounded
areas.  IceModel::basalVelocitySIA() is in iMsia.cc (and is not recommended,
generally).
*/

/*** for ice stream regions (MASK_DRAGGING): ***/
PetscScalar IceModel::basalDragx(PetscScalar **tauc,
                                 PetscScalar **u, PetscScalar **v,
                                 PetscInt i, PetscInt j) const {
  return basal->drag(tauc[i][j], u[i][j], v[i][j]);
}

PetscScalar IceModel::basalDragy(PetscScalar **tauc,
                                 PetscScalar **u, PetscScalar **v,
                                 PetscInt i, PetscInt j) const {
  return basal->drag(tauc[i][j], u[i][j], v[i][j]);
}


//! Initialize the pseudo-plastic till mechanical model.
/*! 
See IceBasalResistancePlasticLaw and updateYieldStressFromHmelt() and getEffectivePressureOnTill()
for model equations.  

Calls either invertSurfaceVelocities(), for one way to get a map of till friction angle
\c vtillphi, or computePhiFromBedElevation() for another way, or leaves \c vtillphi
unchanged.  First two of these are according to options \c -surf_vel_to_phi
and \c -topg_to_phi, respectively.

Also initializes a SIA-type sliding law, but use of that model is not recommended
and is turned off by default.
 */
PetscErrorCode IceModel::initBasalTillModel() {
  PetscErrorCode ierr;

  PetscScalar pseudo_plastic_q = config.get("pseudo_plastic_q"),
    pseudo_plastic_uthreshold = config.get("pseudo_plastic_uthreshold") / secpera,
    plastic_regularization = config.get("plastic_regularization") / secpera;

  bool do_pseudo_plastic_till = config.get_flag("do_pseudo_plastic_till"),
    use_ssa_velocity = config.get_flag("use_ssa_velocity");

  if (basal == NULL)
    basal = new IceBasalResistancePlasticLaw(plastic_regularization, do_pseudo_plastic_till, 
                                 pseudo_plastic_q, pseudo_plastic_uthreshold);

  if (basalSIA == NULL)
    basalSIA = new BasalTypeSIA();  // initialize it; USE NOT RECOMMENDED!
  
  if (use_ssa_velocity) {
    ierr = basal->printInfo(3,grid.com); CHKERRQ(ierr);
  }

  ierr = vtauc.set(config.get("default_tauc")); CHKERRQ(ierr);

  // initialize till friction angle (vtillphi) from options
  PetscTruth  topgphiSet,svphiSet;
  char filename[PETSC_MAX_PATH_LEN];
  ierr = check_option("-topg_to_phi", topgphiSet); CHKERRQ(ierr);
  ierr = PetscOptionsGetString(PETSC_NULL, "-surf_vel_to_phi", filename, 
                               PETSC_MAX_PATH_LEN, &svphiSet); CHKERRQ(ierr);
  if ((svphiSet == PETSC_TRUE) && (topgphiSet == PETSC_TRUE)) {
    SETERRQ(1,"conflicting options for initializing till friction angle; ENDING ...\n");
  }
  if (topgphiSet == PETSC_TRUE) {
    ierr = verbPrintf(2, grid.com, 
      "option -topg_to_phi seen; creating till friction angle map from bed elev ...\n");
      CHKERRQ(ierr);
    // note option -topg_to_phi will be read again to get comma separated array of parameters
    ierr = computePhiFromBedElevation(); CHKERRQ(ierr);
  }
  if (svphiSet == PETSC_TRUE) {
    ierr = verbPrintf(2, grid.com, 
      "option -surf_vel_to_phi seen; doing ad hoc inverse model ...\n"); CHKERRQ(ierr);
    ierr = invertSurfaceVelocities(filename); CHKERRQ(ierr);
  }
  // if neither -surf_vel_to_phi OR -topg_to_phi then pass through; vtillphi is set from
  //   default constant, or -i value, or -boot_from (?)
  return 0;
}


//! Computes the till friction angle phi as a piecewise linear function of bed elevation, according to user options.
/*!
Computes the till friction angle \f$\phi(x,y)\f$ at a location, namely
\c IceModel::vtillphi, as the following increasing, piecewise-linear function of 
the bed elevation \f$b(x,y)\f$.  Let 
	\f[ M = (\phi_{\text{max}} - \phi_{\text{min}}) / (b_{\text{max}} - b_{\text{min}}) \f]
be the slope of the nontrivial part.  Then
	\f[ \phi(x,y) = \begin{cases}
	        \phi_{\text{min}}, & b(x,y) \le b_{\text{min}}, \\
	        \phi_{\text{min}} + (b(x,y) - b_{\text{min}}) \,M,
	                          &  b_{\text{min}} < b(x,y) < b_{\text{max}}, \\
	        \phi_{\text{max}}, & b_{\text{max}} \le b(x,y), \end{cases} \f]
The exception is if the point is marked as floating, in which case the till friction angle
is set to the value \c phi_ocean.

The default values are vaguely suitable for Antarctica, perhaps:
- \c phi_min = 5.0 degrees,
- \c phi_max = 15.0 degrees,
- \c topg_min = -1000.0 m,
- \c topg_max = 1000.0 m,
- \c phi_ocean = 10.0 degrees.
 */
PetscErrorCode IceModel::computePhiFromBedElevation() {

  PetscErrorCode ierr;

  PetscInt    Nparam=5;
  PetscReal   inarray[5] = {5.0, 15.0, -1000.0, 1000.0, 10.0};

  // read comma-separated array of zero to five values
  PetscTruth  topgphiSet;
  ierr = PetscOptionsGetRealArray(PETSC_NULL, "-topg_to_phi", inarray, &Nparam, &topgphiSet);
     CHKERRQ(ierr);
  if (topgphiSet != PETSC_TRUE) {
    SETERRQ(1,"HOW DID I GET HERE? ... ending...\n");
  }
  if (Nparam > 5) {
    ierr = verbPrintf(1, grid.com, 
      "WARNING: option -topg_to_phi read more than 5 parameters ... effect may be bad ...\n");
      CHKERRQ(ierr);
  }
  PetscReal   phi_min = inarray[0],
              phi_max = inarray[1],
              topg_min = inarray[2],
              topg_max = inarray[3],
              phi_ocean = inarray[4];

  ierr = verbPrintf(2, grid.com, 
      "  till friction angle (phi) is piecewise-linear function of bed elev (topg):\n"
      "            /  %5.2f                                 for   topg < %.f\n"
      "      phi = |  %5.2f + (topg - %.f) * (%.2f / %.f)   for   %.f < topg < %.f\n"
      "            \\  %5.2f                                 for   %.f < topg\n",
      phi_min, topg_min,
      phi_min, topg_min, phi_max-phi_min, topg_max - topg_min, topg_min, topg_max,
      phi_max, topg_max);
      CHKERRQ(ierr);

  PetscReal slope = (phi_max - phi_min) / (topg_max - topg_min);
  PetscScalar **tillphi, **bed;
  ierr = vMask.begin_access(); CHKERRQ(ierr);
  ierr = vbed.get_array(bed); CHKERRQ(ierr);
  ierr = vtillphi.get_array(tillphi); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (!vMask.is_floating(i,j)) {
        if (bed[i][j] <= topg_min) {
          tillphi[i][j] = phi_min;
        } else if (bed[i][j] >= topg_max) {
          tillphi[i][j] = phi_max;
        } else {
          tillphi[i][j] = phi_min + (bed[i][j] - topg_min) * slope;
        }
      } else {
        tillphi[i][j] = phi_ocean;
      }
    }
  }
  ierr = vMask.end_access(); CHKERRQ(ierr);
  ierr = vbed.end_access(); CHKERRQ(ierr);
  ierr = vtillphi.end_access(); CHKERRQ(ierr);

  return 0;
}


//! Compute effective pressure on till using effective thickness of stored till water.
/*!
Uses ice thickness to compute overburden pressure.

Provides very simple model of pore water pressure:  Pore water pressure is assumed
to be a fixed fraction of the overburden pressure.

Note \c bwat is thickness of basal water.  It should be zero at points where
base of ice is frozen.

Need \f$0 \le\f$ \c bwat \f$\le\f$ \c max_hmelt before calling this.  There is
no error checking.
 */
PetscScalar IceModel::getEffectivePressureOnTill(PetscScalar thk, PetscScalar bwat,
						 PetscScalar till_pw_fraction,
						 PetscScalar max_hmelt) const {
  const PetscScalar  overburdenP = ice->rho * standard_gravity * thk;
  return overburdenP * (1.0 - till_pw_fraction * (bwat / max_hmelt));
}


//! Update the till yield stress for the pseudo-plastic till model.
/*!
Expanded brief description: Update the till yield stress \e and the mask, for 
the pseudo-plastic till model, based on pressure and stored till water.

We implement formula (2.4) in [\ref SchoofStream] .  That formula is
    \f[   \tau_c = \mu (\rho g H - p_w)\f]
We modify it by:
    - adding a small till cohesion \f$c_0\f$ (see [\ref Paterson] table 8.1);
    - replacing \f$p_w \to \lambda p_w\f$ where \f$\lambda =\f$ 
      Hmelt/DEFAULT_MAX_HMELT; thus \f$0 \le \lambda \le 1\f$ always while 
      \f$\lambda = 0\f$ when the bed is frozen; and
    - computing porewater pressure \f$p_w\f$ as a fixed fraction \f$\varphi\f$ 
      of the overburden pressure \f$\rho g H\f$.
The effective pressure \f$\rho g H - p_w\f$ is actually computed by 
getEffectivePressureOnTill().

With these replacements our formula looks like
    \f[   \tau_c = c_0 + \mu \left(1 - \lambda \varphi\right) \rho g H \f]
Note also that \f$\mu = \tan(\theta)\f$ where \f$\theta\f$ is a "friction angle."
The parameters \f$c_0\f$, \f$\varphi\f$, \f$\theta\f$ can be set by options 
\c -till_cohesion, \c -till_pw_fraction, and \c -till_friction_angle, respectively.
 */
PetscErrorCode IceModel::updateYieldStressFromHmelt() {
  PetscErrorCode  ierr;
  //      (compare the porewater pressure computed by formula (4) in 
  //      C. Ritz et al 2001 J. G. R. vol 106 no D23 pp 31943--31964;
  //      the modification of this porewater pressure as in Lingle&Brown 1987 is not 
  //      implementable because the "elevation of the bed at the grounding line"
  //      is at an unknowable location as we are not doing a flow line model!)

  bool do_plastic_till = config.get_flag("do_plastic_till");
  // only makes sense when do_plastic_till == TRUE
  if (do_plastic_till == PETSC_FALSE) {
    SETERRQ(1,"do_plastic_till == PETSC_FALSE but updateYieldStressFromHmelt() called");
  }

  if (holdTillYieldStress == PETSC_FALSE) { // usual case: use Hmelt to determine tauc
    PetscScalar **tauc, **H, **Hmelt, **tillphi; 

    PetscScalar till_pw_fraction = config.get("till_pw_fraction"),
      till_c_0 = config.get("till_c_0") * 1e3, // convert from kPa to Pa
      till_mu  = tan((pi/180.0)*config.get("default_till_phi")),
      max_hmelt = config.get("max_hmelt");

    ierr =    vMask.begin_access();    CHKERRQ(ierr);
    ierr =    vtauc.get_array(tauc);    CHKERRQ(ierr);
    ierr =       vH.get_array(H);       CHKERRQ(ierr);
    ierr =   vHmelt.get_array(Hmelt);   CHKERRQ(ierr);
    ierr = vtillphi.get_array(tillphi); CHKERRQ(ierr);

    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
        if (vMask.is_floating(i,j)) {
          tauc[i][j] = 0.0;  
        } else if (H[i][j] == 0.0) {
          tauc[i][j] = 1000.0e3;  // large yield stress of 1000 kPa = 10 bar if no ice
        } else { // grounded and there is some ice
          const PetscScalar N = getEffectivePressureOnTill(H[i][j], Hmelt[i][j],
							   till_pw_fraction, max_hmelt);
          if (useConstantTillPhi == PETSC_TRUE) {
            tauc[i][j] = till_c_0 + till_mu * N;
          } else {
            const PetscScalar mymu = tan((pi/180.0) * tillphi[i][j]);
            tauc[i][j] = till_c_0 + mymu * N;
          }
        }
      }
    }
    ierr =    vMask.end_access(); CHKERRQ(ierr);
    ierr =    vtauc.end_access(); CHKERRQ(ierr);
    ierr =       vH.end_access(); CHKERRQ(ierr);
    ierr = vtillphi.end_access(); CHKERRQ(ierr);
    ierr =   vHmelt.end_access(); CHKERRQ(ierr);
  }

  return 0;
}


//! Apply explicit time step for pure diffusion to basal layer of melt water.
/*!
See preprint \ref BBssasliding .

Uses vWork2d[0] to temporarily store new values for Hmelt.
 */
PetscErrorCode IceModel::diffuseHmelt() {
  PetscErrorCode  ierr;
  
  // diffusion constant K in u_t = K \nabla^2 u is chosen so that fundmental
  //   solution has standard deviation \sigma = 20 km at time t = 1000 yrs;
  //   2 \sigma^2 = 4 K t
  const PetscScalar K = 2.0e4 * 2.0e4 / (2.0 * 1000.0 * secpera),
                    Rx = K * dtTempAge / (grid.dx * grid.dx),
                    Ry = K * dtTempAge / (grid.dy * grid.dy);

  // NOTE: restriction that
  //    1 - 2 R_x - 2 R_y \ge 0
  // is a maximum principle restriction; therefore new Hmelt will be between
  // zero and max_hmelt if old Hmelt has that property
  const PetscScalar oneM4R = 1.0 - 2.0 * Rx - 2.0 * Ry;
  if (oneM4R <= 0.0) {
    SETERRQ(1,
       "diffuseHmelt() has 1 - 2Rx - 2Ry <= 0 so explicit method for diffusion unstable\n"
       "  (timestep restriction believed so rare that is not part of adaptive scheme)");
  }

  // communicate ghosted values so neighbors are valid
  ierr = vHmelt.beginGhostComm(); CHKERRQ(ierr);
  ierr = vHmelt.endGhostComm(); CHKERRQ(ierr);

  PetscScalar **Hmelt, **Hmeltnew; 
  ierr = vHmelt.get_array(Hmelt); CHKERRQ(ierr);
  ierr = vWork2d[0].get_array(Hmeltnew); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      Hmeltnew[i][j] = oneM4R * Hmelt[i][j]
                       + Rx * (Hmelt[i+1][j] + Hmelt[i-1][j])
                       + Ry * (Hmelt[i][j+1] + Hmelt[i][j-1]);
    }
  }
  ierr = vHmelt.end_access(); CHKERRQ(ierr);
  ierr = vWork2d[0].end_access(); CHKERRQ(ierr);

  // finally copy new into vHmelt (and communicate ghosted values at the same time)
  ierr = vWork2d[0].beginGhostComm(vHmelt); CHKERRQ(ierr);
  ierr = vWork2d[0].endGhostComm(vHmelt); CHKERRQ(ierr);

  return 0;
}

